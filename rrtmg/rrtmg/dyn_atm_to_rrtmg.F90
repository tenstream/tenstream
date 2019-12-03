!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

!> \page Routines to call tenstream with optical properties from RRTM
!! The routines in this file take care of concatenating the model input
!! with a background profile which is needed to do radiative transfer up
!! until TopOfAtmosphere
!!

module m_dyn_atm_to_rrtmg
#include "petsc/finclude/petsc.h"
  use petsc

  use iso_fortran_env, only: REAL32, REAL64
  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb

  use m_data_parameters, only : iintegers, mpiint, ireals, default_str_len, &
    zero, one, pi, i1, i2, i9, init_mpi_data_parameters

  use m_helper_functions, only: CHKWARN, CHKERR, reverse, &
    imp_allreduce_min, imp_allreduce_max, &
    imp_bcast, read_ascii_file_2d, &
    gradient, get_arg, itoa, ftoa, &
    meanvec, meanval, assert_arr_is_monotonous

  use m_search, only: search_sorted_bisection

  use m_tenstream_interpolation, only : interp_1d

  implicit none

  !logical,parameter :: ldebug=.True.
  logical,parameter :: ldebug=.False.

  ! specific gas constant for dry air [J kg−1 K−1] and standard gravity on earth
  real(ireals), parameter :: Ra  =287.058_ireals, grav  =9.80665_ireals
  real(REAL32), parameter :: Ra32=287.058_REAL32, grav32=9.80665_REAL32
  real(REAL64), parameter :: Ra64=287.058_REAL64, grav64=9.80665_REAL64

  interface
    real function PLKINT(WVLLO, WVLHI, T)
      real :: WVLLO, WVLHI, T
    end function
  end interface

  interface hydrostat_dz
    module procedure hydrostat_dz_real32, hydrostat_dz_real64
  end interface
  interface hydrostat_dp
    module procedure hydrostat_dp_real32, hydrostat_dp_real64
  end interface

  type t_bg_atm
    real(ireals),allocatable, dimension(:) :: plev, tlev, zt, h2o_lev, &
      o3_lev, co2_lev, ch4_lev, n2o_lev, o2_lev, &
      play, zm, dz, tlay, &
      h2o_lay, o3_lay, co2_lay, ch4_lay, n2o_lay, o2_lay
  end type


    type t_tenstr_atm
      type(t_bg_atm),allocatable :: bg_atm
      ! level quantities dim(nlay+1 of merged grid, ncol)
      real(ireals),allocatable :: plev   (:,:) ! interface pressure     [hPa]
      real(ireals),allocatable :: tlev   (:,:) ! interface temperature  [K]
      real(ireals),allocatable :: zt     (:,:) ! interface heights      [m]

      ! layer quantities dim(nlay of merged grid, ncol)
      real(ireals),allocatable :: zm     (:,:) ! layer mean height        [m]
      real(ireals),allocatable :: dz     (:,:) ! vertical layer thickness [m]
      real(ireals),allocatable :: tlay   (:,:) ! layer mean temperature   [K]
      real(ireals),allocatable :: h2o_lay(:,:) ! watervapor volume mixing ratio [e.g. 1e-3]
      real(ireals),allocatable :: o3_lay (:,:) ! ozone volume mixing ratio      [e.g. .1e-6]
      real(ireals),allocatable :: co2_lay(:,:) ! CO2 volume mixing ratio        [e.g. 407e-6]
      real(ireals),allocatable :: ch4_lay(:,:) ! methane volume mixing ratio    [e.g. 2e-6]
      real(ireals),allocatable :: n2o_lay(:,:) ! n2o volume mixing ratio        [e.g. .32]
      real(ireals),allocatable :: o2_lay (:,:) ! oxygen volume mixing ratio     [e.g. .2]
      real(ireals),allocatable :: lwc    (:,:) ! liq water content              [g/kg]
      real(ireals),allocatable :: reliq  (:,:) ! effective radius               [micron]
      real(ireals),allocatable :: iwc    (:,:) ! ice water content              [g/kg]
      real(ireals),allocatable :: reice  (:,:) ! ice effective radius           [micron]
      real(ireals),allocatable :: cfrac  (:,:) ! cloud fraction

      real(ireals),allocatable :: opt_tau(:,:,:) ! optional optical properties: tau, w0, g dim (Nlay_dynamics, Ncol, Nbands(solar or thermal))
      real(ireals),allocatable :: opt_w0 (:,:,:) ! will be added to the rrtmg optical properties
      real(ireals),allocatable :: opt_g  (:,:,:) ! if only tau is allocated, assume it is absorption only


      real(ireals),allocatable :: tskin  (:) ! skin temperature   [K] dim(ncol)

      logical :: lTOA_to_srfc

      ! index of lowermost layer in atm: search for level where height is bigger and pressure is lower:
      ! number of vertical levels of background profile, i.e. 1:atm_ke will be put on top of the dynamics grid
      ! indexing into the atmosphere grid:
      ! 1:d_ke / 1:d_ke1 gives dynamics portion for layer/level variables
      ! d_ke+1: / d_ke1+1: gives background portion for layer/level variables
      ! size(x)-atm_ke: gives background portion
      integer(iintegers), allocatable :: atm_ke
      integer(iintegers) :: d_ke, d_ke1
    end type

    type t_rrtmg_dyn_atm_log_events
      PetscLogEvent :: setup_tenstr_atm
    end type
    type(t_rrtmg_dyn_atm_log_events), allocatable :: logs
  contains

    subroutine setup_tenstr_atm(comm, lTOA_to_srfc, atm_filename, d_plev, d_tlev, atm, &
        d_tlay, d_h2ovmr, d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr,  d_o2vmr, &
        d_lwc, d_reliq, d_iwc, d_reice, d_cloud_fraction, d_surface_height, d_skin_temperature)
      integer(mpiint), intent(in) :: comm
      logical, intent(in) :: lTOA_to_srfc    ! True if provided variables go from TOA to srfc or False if starting at surface

      ! Filename of background atmosphere file. ASCII file with columns:
      ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
      character(default_str_len), intent(in) :: atm_filename

      ! dim(nlay_dynamics+1, ncol)
      real(ireals),intent(in) :: d_plev(:,:) ! pressure on layer interfaces [hPa]
      real(ireals),intent(in) :: d_tlev(:,:) ! Temperature on layer interfaces [K]

      type(t_tenstr_atm), intent(inout) :: atm

      ! all have dim(nlay_dynamics, ncol)
      real(ireals),intent(in),optional :: d_tlay   (:,:)      ! layer mean temperature         [K]
      real(ireals),intent(in),optional :: d_h2ovmr (:,:)      ! watervapor volume mixing ratio [e.g. 1e-3]
      real(ireals),intent(in),optional :: d_o3vmr  (:,:)      ! ozone volume mixing ratio      [e.g. .1e-6]
      real(ireals),intent(in),optional :: d_co2vmr (:,:)      ! CO2 volume mixing ratio        [e.g. 407e-6]
      real(ireals),intent(in),optional :: d_ch4vmr (:,:)      ! methane volume mixing ratio    [e.g. 2e-6]
      real(ireals),intent(in),optional :: d_n2ovmr (:,:)      ! n2o volume mixing ratio        [e.g. .32]
      real(ireals),intent(in),optional :: d_o2vmr  (:,:)      ! oxygen volume mixing ratio     [e.g. .2]
      real(ireals),intent(in),optional :: d_lwc    (:,:)      ! liq water content              [g/kg]
      real(ireals),intent(in),optional :: d_reliq  (:,:)      ! effective radius               [micron]
      real(ireals),intent(in),optional :: d_iwc    (:,:)      ! ice water content              [g/kg]
      real(ireals),intent(in),optional :: d_reice  (:,:)      ! ice effective radius           [micron]
      real(ireals),intent(in),optional :: d_cloud_fraction(:,:) ! cloud fraction
      real(ireals),intent(in),optional :: d_surface_height(:)   ! surface height above sea     [m]
      real(ireals),intent(in),optional :: d_skin_temperature(:) ! skin tempertaure             [K]

      integer(iintegers) :: icol
      PetscClassId, parameter :: cid=0
      integer(mpiint) :: ierr

      if(lTOA_to_srfc) then
        call CHKERR(1_mpiint, 'currently not possible to supply dynamics input starting at the TOP,'// &
          'input should be starting at the surface')
      endif

      call init_mpi_data_parameters(comm)
      if(.not.allocated(logs)) then
        allocate(logs)
        call PetscLogEventRegister('setup_tenstr_atm', cid, logs%setup_tenstr_atm, ierr); call CHKERR(ierr)
      endif

      call PetscLogEventBegin(logs%setup_tenstr_atm, ierr); call CHKERR(ierr)

      if(.not.allocated(atm%bg_atm)) then
        call load_atmfile(comm, atm_filename, atm%bg_atm)
        call sanitize_input(.True., atm%bg_atm%plev, atm%bg_atm%tlev, ierr, atm%bg_atm%tlay)
        call CHKERR(ierr, 'bad input in bg_atmosphere file')
      endif

      call check_shape_2d(d_tlev          ,size(d_plev, 1, kind=iintegers)   , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_tlay          ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_h2ovmr        ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_o3vmr         ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_co2vmr        ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_ch4vmr        ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_n2ovmr        ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_o2vmr         ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_lwc           ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_reliq         ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_iwc           ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_reice         ,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_cloud_fraction,size(d_plev, 1, kind=iintegers)-1 , size(d_plev, 2, kind=iintegers))
      call check_shape_1d(d_surface_height, ncol=size(d_plev, 2, kind=iintegers))

      do icol=lbound(d_plev,2),ubound(d_plev,2)
        if(present(d_tlay)) then
          call sanitize_input(lTOA_to_srfc, d_plev(:,icol), d_tlev(:,icol), ierr, d_tlay(:,icol))
        else
          call sanitize_input(lTOA_to_srfc, d_plev(:,icol), d_tlev(:,icol), ierr)
        endif
        call CHKERR(ierr, 'bad input from dynamics grid, column: '//itoa(icol))
      enddo

      call merge_dyn_rad_grid(comm, atm, &
        d_plev, d_tlev, d_tlay, d_h2ovmr, &
        d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr, &
        d_o2vmr, d_lwc, d_reliq, d_iwc, d_reice, &
        d_cfrac=d_cloud_fraction, &
        d_surface_height=d_surface_height)

      call check_shape_2d(d_tlev          , atm%d_ke1,size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_tlay          , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_h2ovmr        , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_o3vmr         , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_co2vmr        , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_ch4vmr        , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_n2ovmr        , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_o2vmr         , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_lwc           , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_reliq         , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_iwc           , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_reice         , atm%d_ke, size(d_plev, 2, kind=iintegers))
      call check_shape_2d(d_cloud_fraction, atm%d_ke, size(d_plev, 2, kind=iintegers))

      call check_shape_1d(d_skin_temperature, size(d_plev, 2, kind=iintegers))
      if(present(d_skin_temperature)) then
        if(.not.allocated(atm%tskin)) allocate(atm%tskin(size(d_plev, 2, kind=iintegers)))
        atm%tskin = d_skin_temperature
      endif

      call PetscLogEventEnd(logs%setup_tenstr_atm, ierr); call CHKERR(ierr)
      contains
        subroutine check_shape_1d(d_arr, ncol)
          real(ireals), intent(in), optional :: d_arr(:)
          integer(iintegers), intent(in), optional :: ncol

          if(present(d_arr)) then
            if(present(ncol)) call CHKERR(int(size(d_arr,1)-ncol, mpiint), &
              'bad nr cols got'//itoa(size(d_arr,1))//' expect '//itoa(ncol))
          endif
        end subroutine
        subroutine check_shape_2d(d_arr, k, ncol)
          real(ireals), intent(in), optional :: d_arr(:,:)
          integer(iintegers), intent(in), optional :: k, ncol

          if(present(d_arr)) then
            if(present(k))    call CHKERR(int(size(d_arr,1)-k   , mpiint), &
              'bad vert size! got '//itoa(size(d_arr,1))//' but expected '//itoa(k))
            if(present(ncol)) call CHKERR(int(size(d_arr,2)-ncol, mpiint), &
              'bad nr cols got '//itoa(size(d_arr,2))//' but expected '//itoa(ncol))
          endif
        end subroutine
    end subroutine

    subroutine print_tenstr_atm(atm, icol)
      type(t_tenstr_atm),intent(in) :: atm
      integer(iintegers), optional :: icol
      integer(iintegers) :: k, j

      j = get_arg(i1, icol)

      print *,'atm%bg_atm ', allocated(atm%bg_atm)
      call alloc_info(atm%plev   ,'atm%plev   ')
      call alloc_info(atm%tlev   ,'atm%tlev   ')
      call alloc_info(atm%zt     ,'atm%zt     ')

      call alloc_info(atm%zm     ,'atm%zm     ')
      call alloc_info(atm%dz     ,'atm%dz     ')
      call alloc_info(atm%tlay   ,'atm%tlay   ')
      call alloc_info(atm%h2o_lay,'atm%h2o_lay')
      call alloc_info(atm%o3_lay ,'atm%o3_lay ')
      call alloc_info(atm%co2_lay,'atm%co2_lay')
      call alloc_info(atm%ch4_lay,'atm%ch4_lay')
      call alloc_info(atm%n2o_lay,'atm%n2o_lay')
      call alloc_info(atm%o2_lay ,'atm%o2_lay ')
      call alloc_info(atm%lwc    ,'atm%lwc    ')
      call alloc_info(atm%reliq  ,'atm%reliq  ')
      call alloc_info(atm%iwc    ,'atm%iwc    ')
      call alloc_info(atm%reice  ,'atm%reice  ')
      print *,'atm%atm_ke ', allocated(atm%atm_ke)

      do k=size(atm%plev,1),1,-1
        print *,k,'zt',atm%zt(k,j), 'plev', atm%plev(k,j), 'Tlev', atm%tlev(k,j)
      enddo

      do k=size(atm%tlay,1),1,-1
        print *,k,'dz',atm%dz(k,j), 'Tlay', atm%tlay(k,j), &
          'H2O', atm%h2o_lay(k,j), 'CO2', atm%co2_lay(k,j),'O3', atm%o3_lay(k,j),'N2O', atm%n2o_lay(k,j),&
          'O2', atm%o2_lay(k,j)
      enddo
    contains
      subroutine alloc_info(x, varname)
        real(ireals), allocatable, dimension(:,:), intent(in) :: x
        character(len=*), intent(in) :: varname
        if(allocated(x)) then
          print *,varname, allocated(x), shape(x)
        else
          print *,varname, allocated(x)
        endif
      end subroutine
    end subroutine


    subroutine destroy_tenstr_atm(atm)
      type(t_tenstr_atm),intent(inout) :: atm
      if(allocated(atm%bg_atm )) deallocate(atm%bg_atm )
      if(allocated(atm%plev   )) deallocate(atm%plev   )
      if(allocated(atm%tlev   )) deallocate(atm%tlev   )
      if(allocated(atm%zt     )) deallocate(atm%zt     )

      if(allocated(atm%zm     )) deallocate(atm%zm     )
      if(allocated(atm%dz     )) deallocate(atm%dz     )
      if(allocated(atm%tlay   )) deallocate(atm%tlay   )
      if(allocated(atm%tskin  )) deallocate(atm%tskin  )
      if(allocated(atm%h2o_lay)) deallocate(atm%h2o_lay)
      if(allocated(atm%o3_lay )) deallocate(atm%o3_lay )
      if(allocated(atm%co2_lay)) deallocate(atm%co2_lay)
      if(allocated(atm%ch4_lay)) deallocate(atm%ch4_lay)
      if(allocated(atm%n2o_lay)) deallocate(atm%n2o_lay)
      if(allocated(atm%o2_lay )) deallocate(atm%o2_lay )
      if(allocated(atm%lwc    )) deallocate(atm%lwc    )
      if(allocated(atm%reliq  )) deallocate(atm%reliq  )
      if(allocated(atm%iwc    )) deallocate(atm%iwc    )
      if(allocated(atm%reice  )) deallocate(atm%reice  )
      if(allocated(atm%cfrac  )) deallocate(atm%cfrac  )
      if(allocated(atm%atm_ke )) deallocate(atm%atm_ke )
    end subroutine

    !> Reorder the atmosphere object so that it either goes from TOA to surface or the other way around
    !! atm%bg_atm always goes from TOA to srfc
    subroutine reorder_atm_vertically(atm, lTOA_to_srfc)
      type(t_tenstr_atm),intent(inout) :: atm
      logical, intent(in) :: lTOA_to_srfc

      if(atm%lTOA_to_srfc.eqv.lTOA_to_srfc) return ! nothing to do

      if(allocated(atm%plev   )) atm%plev    = reverse(atm%plev   )
      if(allocated(atm%tlev   )) atm%tlev    = reverse(atm%tlev   )
      if(allocated(atm%zt     )) atm%zt      = reverse(atm%zt     )

      if(allocated(atm%zm     )) atm%zm      = reverse(atm%zm     )
      if(allocated(atm%dz     )) atm%dz      = reverse(atm%dz     )
      if(allocated(atm%tlay   )) atm%tlay    = reverse(atm%tlay   )
      if(allocated(atm%h2o_lay)) atm%h2o_lay = reverse(atm%h2o_lay)
      if(allocated(atm%o3_lay )) atm%o3_lay  = reverse(atm%o3_lay )
      if(allocated(atm%co2_lay)) atm%co2_lay = reverse(atm%co2_lay)
      if(allocated(atm%ch4_lay)) atm%ch4_lay = reverse(atm%ch4_lay)
      if(allocated(atm%n2o_lay)) atm%n2o_lay = reverse(atm%n2o_lay)
      if(allocated(atm%o2_lay )) atm%o2_lay  = reverse(atm%o2_lay )
      if(allocated(atm%lwc    )) atm%lwc     = reverse(atm%lwc    )
      if(allocated(atm%reliq  )) atm%reliq   = reverse(atm%reliq  )
      if(allocated(atm%iwc    )) atm%iwc     = reverse(atm%iwc    )
      if(allocated(atm%reice  )) atm%reice   = reverse(atm%reice  )
      if(allocated(atm%cfrac  )) atm%cfrac   = reverse(atm%cfrac  )

      atm%lTOA_to_srfc = lTOA_to_srfc
    end subroutine

    !> Concatenate the dynamical grid with the background profile
    subroutine merge_dyn_rad_grid(comm, atm,   &
        d_plev, d_tlev, d_tlay, d_h2ovmr,   &
        d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr, &
        d_o2vmr, d_lwc, d_reliq, d_iwc, d_reice, &
        d_cfrac, &
        d_surface_height)

      integer(mpiint), intent(in) :: comm
      type(t_tenstr_atm),intent(inout) :: atm

      real(ireals),intent(in) :: d_plev (:,:), d_tlev(:,:) ! dim(nlay_dynamics+1, ncol)

      real(ireals),intent(in),optional :: d_tlay   (:,:) ! all have
      real(ireals),intent(in),optional :: d_h2ovmr (:,:) ! dim(nlay_dynamics, ncol)
      real(ireals),intent(in),optional :: d_o3vmr  (:,:) !
      real(ireals),intent(in),optional :: d_co2vmr (:,:) !
      real(ireals),intent(in),optional :: d_ch4vmr (:,:) !
      real(ireals),intent(in),optional :: d_n2ovmr (:,:) !
      real(ireals),intent(in),optional :: d_o2vmr  (:,:) !
      real(ireals),intent(in),optional :: d_lwc    (:,:) !
      real(ireals),intent(in),optional :: d_reliq  (:,:) !
      real(ireals),intent(in),optional :: d_iwc    (:,:) !
      real(ireals),intent(in),optional :: d_reice  (:,:) !
      real(ireals),intent(in),optional :: d_cfrac  (:,:) !
      real(ireals),intent(in),optional :: d_surface_height(:)

      !real(ireals) :: d_plev (ubound(in_d_plev,1), ubound(in_d_plev,2))

      integer(iintegers) :: ke, ke1 ! number of vertical levels of merged grid
      integer(iintegers) :: is, ie, icol, l,m

      real(ireals),allocatable :: d_hhl(:,:), d_dz(:)
      real(ireals) :: global_maxheight, global_minplev
      logical :: lupdate_bg_entries
      real(ireals) :: hsrfc, dz, interface_pressures(2), tmp_lay_temp
      real(ireals) :: minval_plev, maxval_hhl

      hsrfc = zero ! default value

      if(.not.allocated(atm%bg_atm)) call CHKERR(1_mpiint,'bg_atm has to be allocated before merging with dynamics grid variables')
      associate( bg_atm => atm%bg_atm )

        is = lbound(d_plev,2); ie = ubound(d_plev,2)

        if(allocated(atm%atm_ke)) then
          call CHKERR(int(atm%d_ke1-ubound(d_plev,1),mpiint), &
            'Seems you changed the vertical dimension of the input between calls.'// &
            'atm_d_ke1 '//itoa(atm%d_ke1)//' d_plev '//itoa(ubound(d_plev,1))// &
            ' You have to destroy the atmosphere first and recreate')
        endif

        ! find out how many layers we have to put on top of the dynamics grid
        if(.not.allocated(atm%atm_ke)) then
          atm%d_ke1 = ubound(d_plev,1); atm%d_ke = atm%d_ke1-1

          ! First get top height of dynamics grid
          allocate(d_hhl(atm%d_ke1, ie))
          allocate(d_dz(atm%d_ke))

          maxval_hhl = zero
          do icol = is, ie
            if(present(d_surface_height)) then
              hsrfc = d_surface_height(icol)
            endif
            if(present(d_tlay)) then
              call hydrostat_lev(d_plev(:,icol),d_tlay(:,icol), hsrfc, d_hhl(:, icol), d_dz)
            else
              call hydrostat_lev(d_plev(:,icol),(d_tlev(1:atm%d_ke,icol)+d_tlev(2:atm%d_ke1,icol))/2, &
                hsrfc, d_hhl(:, icol), d_dz)
            endif
            maxval_hhl = max(maxval_hhl, d_hhl(size(d_hhl,dim=1),icol))
          enddo

          ! first put a tiny increment on the top of dynamics pressure value,
          ! to handle a cornercase where dynamics and background profile are the same
          ! because sometimes it happens that the layer between dynamics grid and profile
          ! is so small that dz is very small and consequently not in LUT's
          minval_plev = minval(d_plev(size(d_plev, dim=1), :))
          m = floor(search_sorted_bisection(bg_atm%plev, minval_plev))
          interface_pressures(2) = bg_atm%plev(m)
          tmp_lay_temp = (bg_atm%tlev(m) + bg_atm%tlev(m+1)) * .5_ireals
          minval_plev = d_plev(size(d_plev, dim=1), is)
          do icol = is, ie
            interface_pressures(1) = d_plev(size(d_plev, dim=1),icol)
            dz = hydrostat_dz(abs(interface_pressures(1)-interface_pressures(2)), &
              meanval(interface_pressures), tmp_lay_temp )
            if(dz.lt.one) then
              !call CHKWARN(1_mpiint, 'bg atmosphere and dynamics grid pressure are very close.' // &
              !  'Note that I`ll drop one layer here.')
              minval_plev = min(minval_plev, (bg_atm%plev(m) + bg_atm%plev(max(i1,m-i1)))*.5_ireals)
            endif
          enddo

          ! index of lowermost layer in atm: search for level where height is bigger and
          ! pressure is lower
          call imp_allreduce_max(comm, maxval_hhl, global_maxheight)
          call imp_allreduce_min(comm, minval_plev, global_minplev)

          if(global_maxheight.ge.bg_atm%zt(1)) &
            call CHKERR(1_mpiint, 'background profile TOA height is smaller than dynamics grid height')
          if(global_minplev.le.bg_atm%plev(1)) &
            call CHKERR(1_mpiint, 'background profile TOA pressure is larger than dynamics grid pressure')

          l = floor(search_sorted_bisection(bg_atm%zt, global_maxheight))
          m = floor(search_sorted_bisection(bg_atm%plev, global_minplev))
          allocate(atm%atm_ke)
          atm%atm_ke = min(l,m)
          ke  = atm%atm_ke + atm%d_ke
          ke1 = atm%atm_ke + atm%d_ke1

          ! then from there on couple background atm data on top of that
          if(.not.allocated(atm%plev   )) allocate(atm%plev   (ke1, ie))
          if(.not.allocated(atm%tlev   )) allocate(atm%tlev   (ke1, ie))
          if(.not.allocated(atm%zt     )) allocate(atm%zt     (ke1, ie))
          if(.not.allocated(atm%dz     )) allocate(atm%dz     (ke,  ie))
          if(.not.allocated(atm%tlay   )) allocate(atm%tlay   (ke,  ie))
          if(.not.allocated(atm%h2o_lay)) allocate(atm%h2o_lay(ke,  ie))
          if(.not.allocated(atm%o3_lay )) allocate(atm%o3_lay (ke,  ie))
          if(.not.allocated(atm%co2_lay)) allocate(atm%co2_lay(ke,  ie))
          if(.not.allocated(atm%ch4_lay)) allocate(atm%ch4_lay(ke,  ie))
          if(.not.allocated(atm%n2o_lay)) allocate(atm%n2o_lay(ke,  ie))
          if(.not.allocated(atm%o2_lay )) allocate(atm%o2_lay (ke,  ie))
          if(.not.allocated(atm%lwc    )) allocate(atm%lwc    (ke,  ie))
          if(.not.allocated(atm%reliq  )) allocate(atm%reliq  (ke,  ie))
          if(.not.allocated(atm%iwc    )) allocate(atm%iwc    (ke,  ie))
          if(.not.allocated(atm%reice  )) allocate(atm%reice  (ke,  ie))

          lupdate_bg_entries = .True.
        else
          lupdate_bg_entries = .False.
        endif

        call alloc_if_present(d_cfrac, atm%cfrac, [size(atm%lwc,dim=1,kind=iintegers), size(atm%lwc,dim=2,kind=iintegers)])

        associate(atm_ke => atm%atm_ke)
          ke  = atm_ke + atm%d_ke
          ke1 = atm_ke + atm%d_ke1

          do icol=is,ie

            ! First merge pressure levels .. pressure is always given..
            if(lupdate_bg_entries) then
              atm%plev(ke1-atm_ke+1:ke1, icol) = reverse(bg_atm%plev(1:atm_ke))
            endif
            atm%plev(1:atm%d_ke1, icol) = d_plev(:, icol)
            if(atm%plev(ke1-atm_ke+1, icol) .gt. atm%plev(atm%d_ke1, icol)) then
              print *,'background profile pressure is .ge. than uppermost pressure &
                & level of dynamics grid -- this suggests the dynamics grid is way &
                & off hydrostatic balance... please check', icol, &
                atm%plev(ke1-atm_ke+1, icol), atm%plev(atm%d_ke1, icol), atm%plev(:, icol)
              call CHKWARN(1_mpiint, 'error in rrtm_lw merging grids')
            endif

            ! And also Tlev has to be present always
            if(lupdate_bg_entries) then
              atm%tlev(ke1-atm_ke+1:ke1, icol) = reverse(bg_atm%tlev(1:atm_ke))
            endif
            atm%tlev(1:atm%d_ke1, icol) = d_tlev(:,icol)

            ! compute dz
            if(present(d_surface_height)) hsrfc = d_surface_height(icol)
            if(lupdate_bg_entries) then
              call hydrostat_lev(atm%plev(:,icol),(atm%tlev(1:ke,icol)+atm%tlev(2:ke1,icol))/2, &
                hsrfc, atm%zt(:, icol), atm%dz(:, icol))
            else
              call hydrostat_lev(atm%plev(1:atm%d_ke1,icol),(d_tlev(1:atm%d_ke,icol)+d_tlev(2:atm%d_ke1,icol))/2, &
                hsrfc, atm%zt(1:atm%d_ke1, icol), atm%dz(1:atm%d_ke, icol))
            endif


            if(present(d_tlay)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%tlay, bg_atm%tlev, atm%tlay(:, icol), d_tlay(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%tlay, bg_atm%tlev, atm%tlay(:, icol), &
                (d_tlev(1:atm%d_ke,icol)+d_tlev(2:atm%d_ke1,icol))/2)
            endif

            if(present(d_lwc)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                zero*bg_atm%tlay, zero*bg_atm%tlev, atm%lwc(:,icol), d_lwc(:,icol))
            else
              atm%lwc(:,icol) = zero
            endif
            if(present(d_reliq)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                zero*bg_atm%tlay, zero*bg_atm%tlev, atm%reliq(:,icol), d_reliq(:,icol))
            else
              atm%reliq = zero
            endif

            if(present(d_iwc)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                zero*bg_atm%tlay, zero*bg_atm%tlev, atm%iwc(:,icol), d_iwc(:,icol))
            else
              atm%iwc(:,icol) = zero
            endif
            if(present(d_reice)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                zero*bg_atm%tlay, zero*bg_atm%tlev, atm%reice(:,icol), d_reice(:,icol))
            else
              atm%reice = zero
            endif

            if(present(d_cfrac)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                zero*bg_atm%tlay, zero*bg_atm%tlev, atm%cfrac(:,icol), d_cfrac(:,icol))
            endif

            if(present(d_h2ovmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%h2o_lay, bg_atm%h2o_lev, atm%h2o_lay(:,icol), d_h2ovmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
              bg_atm%h2o_lay, bg_atm%h2o_lev, atm%h2o_lay(:,icol))
            endif
            if(present(d_o3vmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%o3_lay, bg_atm%o3_lev, atm%o3_lay(:,icol), d_o3vmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%o3_lay, bg_atm%o3_lev, atm%o3_lay(:,icol))
            endif

            if(present(d_co2vmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%co2_lay, bg_atm%co2_lev, atm%co2_lay(:,icol), d_co2vmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%co2_lay, bg_atm%co2_lev, atm%co2_lay(:,icol))
            endif
            if(present(d_ch4vmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%ch4_lay, bg_atm%ch4_lev, atm%ch4_lay(:,icol), d_ch4vmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%ch4_lay, bg_atm%ch4_lev, atm%ch4_lay(:,icol))
            endif
            if(present(d_n2ovmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%n2o_lay, bg_atm%n2o_lev, atm%n2o_lay(:,icol), d_n2ovmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%n2o_lay, bg_atm%n2o_lev, atm%n2o_lay(:,icol))
            endif
            if(present(d_o2vmr)) then
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%o2_lay, bg_atm%o2_lev, atm%o2_lay(:,icol), d_o2vmr(:,icol))
            else
              call merge_grid_var(lupdate_bg_entries, bg_atm%zt, atm%zt(:,icol), atm_ke, &
                bg_atm%o2_lay, bg_atm%o2_lev, atm%o2_lay(:,icol))
            endif
          enddo
        end associate
      end associate
      atm%lTOA_to_srfc = .False.
      contains
        subroutine alloc_if_present(opt_input_arr, atm_arr, dims)
          real(ireals), optional, intent(in) :: opt_input_arr(:,:)
          real(ireals), allocatable, intent(inout) :: atm_arr(:,:)
          integer(iintegers), intent(in) :: dims(2)
          if(present(opt_input_arr)) then
            if(.not.allocated(atm_arr)) allocate(atm_arr(dims(1),dims(2)))
          endif
        end subroutine
    end subroutine

    ! merge the dynamics grid and the background profile together at lvl atm_ke
    ! NOTE! Only use with variables on layer
    subroutine merge_grid_var(lupdate_bg_entries, a_hhl, d_hhl, atm_ke, a_lay, a_lev, col_var, d_var)
      logical, intent(in) :: lupdate_bg_entries
      integer(iintegers),intent(in) :: atm_ke
      real(ireals),intent(in) :: a_hhl(:), d_hhl(:), a_lay(:), a_lev(:) ! a_arr is from atm%, d_arr corresponds to dynamics grids
      real(ireals),intent(out) :: col_var(:)
      real(ireals),intent(in),optional :: d_var(:)
      integer(iintegers) :: k, kt ! kt is reverse index
      real(ireals) :: h

      ! Top of atmosphere layers are always given by background profile
      if(lupdate_bg_entries) then
        do k=1,atm_ke
          kt = size(col_var)-k+1
          col_var(kt) = a_lay(k)
        enddo
      endif

      if(present(d_var)) then ! dynamics grid variable is provided, use that
        do k=1,size(d_var)
          col_var(k) = d_var(k)
        enddo
      else ! we may still use atmospheric grid file instead...
        if(lupdate_bg_entries) then
          do k=1,size(col_var)-atm_ke
            h = (d_hhl(k+1) + d_hhl(k)) / 2
            col_var(k) = interp_1d(search_sorted_bisection(a_hhl, h), a_lev)
          enddo
        endif
      endif
    end subroutine

    subroutine sanitize_input(lTOA_to_srfc, plev, tlev, ierr, tlay)
      logical, intent(in) :: lTOA_to_srfc
      real(ireals),intent(in),dimension(:) :: plev, tlev
      integer(mpiint), intent(out) :: ierr
      real(ireals),intent(in),dimension(:),optional :: tlay
      integer(iintegers) :: k
      logical :: lerr

      ierr = 0

      lerr = .not.assert_arr_is_monotonous(plev, lincreasing=lTOA_to_srfc, lstrict=.True.)
      if(lerr) then
        print *,'Pressure is not strictly monotous decreasing, however, '// &
          'we need this for hydrostatic integration. '// &
          'If you cannot guarantee that, please ask me for non-hydrostatic support.', plev
        ierr = ierr+1
      endif

      if(lTOA_to_srfc) then
        lerr = plev(size(plev)) .gt. 1050
      else
        lerr = plev(1) .gt. 1050
      endif
      if(lerr) then
        print *,'Pressure above 1050 hPa -- are you sure this is earth?', maxval(plev)
        ierr = ierr+1
      endif

      if(lTOA_to_srfc) then
        lerr = plev(1) .lt. zero
      else
        lerr = plev(size(plev)) .lt. zero
      endif
      if(lerr) then
        print *,'Pressure negative -- are you sure this is physically correct?', minval(plev)
        ierr = ierr+1
      endif

      lerr = minval(tlev) .lt. 159
      if(lerr) then
        print *,'Temperature is very low -- are you sure RRTMG can handle that?', minval(tlev)
        do k=lbound(tlev,1), ubound(tlev,1)
          print *,'lev',k,'T',tlev(k)
        enddo
        ierr = ierr+1
      endif

      lerr = maxval(tlev) .gt. 400
      if(lerr) then
        print *,'Temperature is very high -- are you sure RRTMG can handle that?', maxval(tlev)
        do k=lbound(tlev,1), ubound(tlev,1)
          print *,'lev',k,'T',tlev(k)
        enddo
        ierr = ierr+1
      endif

      if(present(tlay) .and. ldebug) then
        do k=lbound(tlay,1), ubound(tlay,1)
          lerr = (tlay(k)-tlev(k).ge.zero) .eqv. (tlay(k)-tlev(k+1).gt.zero) ! different sign says its in between
          if(lerr) then
            print *,'Layer Temperature not between level temps?', k, tlev(k), '|', tlay(k), '|', tlev(k+1)
            ierr = ierr+1
          endif
        enddo
      endif

      if(ierr.gt.0) then
        print *,'Found wonky input to pprts_rrtm_lw -- please check! -- you should probably abort now.'
      endif
    end subroutine

    subroutine load_atmfile(comm, atm_filename, atm)
      integer(mpiint), intent(in) :: comm
      character(default_str_len), intent(in) :: atm_filename
      type(t_bg_atm),allocatable,intent(inout) :: atm

      integer(mpiint) :: myid, ierr
      integer(iintegers) :: k, nlev
      real(ireals),allocatable :: prof(:,:) ! # z(km)  p(mb)  T(K) air(cm-3) o3(cm-3) o2(cm-3)  h2o(cm-3) co2(cm-3) no2(cm-3)

      if(allocated(atm)) return

      call mpi_comm_rank(comm, myid, ierr)
      allocate(atm)

      if(myid.eq.0) then
        call read_ascii_file_2d(atm_filename, prof, ierr)
        if(ierr.ne.0) then
          print *,'************* Error occured reading the atmosphere file:', atm_filename, '::', ierr
          call CHKERR(ierr)
        endif

        nlev = ubound(prof,1)

        allocate(atm%plev   (nlev))
        allocate(atm%zt     (nlev))
        allocate(atm%tlev   (nlev))
        allocate(atm%h2o_lev(nlev))
        allocate(atm%o3_lev (nlev))
        allocate(atm%co2_lev(nlev))
        allocate(atm%ch4_lev(nlev))
        allocate(atm%n2o_lev(nlev))
        allocate(atm%o2_lev (nlev))

        atm%zt   = prof(:,1)*1e3
        atm%plev = prof(:,2)
        atm%tlev = prof(:,3)
        atm%h2o_lev = prof(:,7) / prof(:,4)
        atm%o3_lev  = prof(:,5) / prof(:,4)
        atm%co2_lev = prof(:,8) / prof(:,4)
        atm%ch4_lev = atm%co2_lev / 1e2
        atm%n2o_lev = prof(:,9) / prof(:,4)
        atm%o2_lev  = prof(:,6) / prof(:,4)

        if(ldebug .and. myid.eq.0) then
          do k=1, nlev
            print *,k,'zt', atm%zt(k), 'plev', atm%plev(k), 'T', atm%tlev(k), 'CO2', atm%co2_lev(k), &
              'H2O', atm%h2o_lev(k), 'O3', atm%o3_lev(k),'N2O' , atm%n2o_lev(k), 'O2', atm%o2_lev(k)
          enddo

        endif
      endif
      call imp_bcast(comm, atm%plev   , 0_mpiint)
      call imp_bcast(comm, atm%zt     , 0_mpiint)
      call imp_bcast(comm, atm%tlev   , 0_mpiint)
      call imp_bcast(comm, atm%h2o_lev, 0_mpiint)
      call imp_bcast(comm, atm%o3_lev , 0_mpiint)
      call imp_bcast(comm, atm%co2_lev, 0_mpiint)
      call imp_bcast(comm, atm%ch4_lev, 0_mpiint)
      call imp_bcast(comm, atm%n2o_lev, 0_mpiint)
      call imp_bcast(comm, atm%o2_lev , 0_mpiint)

      nlev = size(atm%plev)

      allocate(atm%play   (nlev-1))
      allocate(atm%zm     (nlev-1))
      allocate(atm%dz     (nlev-1))
      allocate(atm%tlay   (nlev-1))
      allocate(atm%h2o_lay(nlev-1))
      allocate(atm%o3_lay (nlev-1))
      allocate(atm%co2_lay(nlev-1))
      allocate(atm%ch4_lay(nlev-1))
      allocate(atm%n2o_lay(nlev-1))
      allocate(atm%o2_lay (nlev-1))

      atm%play    = meanvec(atm%plev   )
      atm%zm      = meanvec(atm%zt     )
      atm%dz      = -gradient(atm%zt   )
      atm%tlay    = meanvec(atm%tlev   )
      atm%h2o_lay = meanvec(atm%h2o_lev)
      atm%o3_lay  = meanvec(atm%o3_lev )
      atm%co2_lay = meanvec(atm%co2_lev)
      atm%ch4_lay = meanvec(atm%ch4_lev)
      atm%n2o_lay = meanvec(atm%n2o_lev)
      atm%o2_lay  = meanvec(atm%o2_lev )
    end subroutine

    pure elemental function hydrostat_dz_real32(dp, p, T)
      real(REAL32), intent(in) :: dp, p, T
      real(REAL32) :: hydrostat_dz_real32, rho
      rho = p / Ra32 / T
      hydrostat_dz_real32 = dp / rho / grav32
    end function
    pure elemental function hydrostat_dz_real64(dp, p, T)
      real(REAL64), intent(in) :: dp, p, T
      real(REAL64) :: hydrostat_dz_real64, rho
      rho = p / Ra64 / T
      hydrostat_dz_real64 = dp / rho / grav64
    end function

    subroutine hydrostat_lev(plev,tlay, hsrfc, hhl, dz)
      ! Integrate vertical height profile hydrostatically. arrays start at bottom(surface)
      real(ireals),intent(in) :: plev(:),tlay(:)
      real(ireals),intent(in) :: hsrfc
      real(ireals),intent(out) :: hhl(:), dz(:)
      integer(im) :: k
      hhl(1) = hsrfc
      do k=1,size(tlay)
        dz(k) = hydrostat_dz(abs(plev(k+1)-plev(k)), (plev(k+1)+plev(k))/2, tlay(k))
        hhl(k+1) = hhl(k) + dz(k)
      enddo
      if(any(dz.le.zero)) then
        print *,'plev',plev
        print *,'tlay',tlay
        print *,'dz',dz
        print *,'hhl',hhl
        call CHKERR(1_mpiint, 'error in dz, dz.le.zero minval:'//ftoa(minval(dz)))
      endif
    end subroutine

    pure elemental function hydrostat_dp_real32(dz, p, T)
      real(REAL32), intent(in) :: dz, p, T
      real(REAL32) :: hydrostat_dp_real32, rho
      rho = p / Ra32 / T
      hydrostat_dp_real32 = dz * rho * grav32
    end function
    pure elemental function hydrostat_dp_real64(dz, p, T)
      real(REAL64), intent(in) :: dz, p, T
      real(REAL64) :: hydrostat_dp_real64, rho
      rho = p / Ra64 / T
      hydrostat_dp_real64 = dz * rho * grav64
    end function

    subroutine hydrostat_plev(psrfc, tlay, hhl, plev, dp)
      ! Integrate vertical height profile hydrostatically. arrays start at bottom(surface)
      real(ireals),intent(in) :: psrfc, tlay(:), hhl(:)
      real(ireals),intent(out) :: plev(:), dp(:)
      integer(im) :: k
      plev(1:2) = psrfc
      dp(1) = zero
      do k=1,size(tlay)
        dp(k) = hydrostat_dp(abs(hhl(k+1)-hhl(k)), plev(k)-dp(max(1,k-1))/2, tlay(k))
        plev(k+1) = plev(k) - dp(k)
      enddo
      if(any(dp.le.zero)) then
        do k = 1, size(tlay)
          print *,'k', k, 'hhl',hhl(k), 'plev',plev(k), 'tlay',tlay(k), 'dp',dp(k)
        enddo
        k = size(hhl)
        print *,'k', k, 'hhl', hhl(k), 'plev', plev(k)
        stop 'error in dp'
      endif
    end subroutine

    ! Compute effective radius from liquid water content and droplet density
    ! after Martin et al (1993)
    ! and Bugliaro et al (2011)
    pure elemental function reff_from_lwc_and_N(lwc, N, k) result(reff)
      real(ireals), intent(in) :: lwc ! liquid water content [g m-3],  typically .1
      real(ireals), intent(in) :: N   ! droplet density      [1 cm-3], typically 100
      real(ireals), intent(in), optional :: k   ! airmass factor (maritime ~ .8 +- 0.07) (continental ~ .67 +- 0.07)
      real(ireals) :: reff            ! effective radius [1e-6 m]

      real(ireals), parameter :: &
        rho = 1e3_ireals, &    ! water density @ 4degC
        default_k = .75_ireals ! airmass factor in between

      if(present(k)) then
        reff = 1e3_ireals * (3._ireals * lwc / (4._ireals*pi*rho*k*max(tiny(N),N))) ** (one/3._ireals)
      else
        reff = 1e3_ireals * (3._ireals * lwc / (4._ireals*pi*rho*default_k*max(tiny(N),N))) ** (one/3._ireals)
      endif
    end function

    ! Compute effective radius from liquid water content and droplet density
    ! from the implementation in ICON/newcld_optics
    ! see ECHAM5 documentation (Roeckner et al, MPI report 349)
    pure elemental function reff_from_lwc_and_N_after_ICON(lwc, N, l_liquid, zkap) result(reff)
      real(ireals), intent(in) :: lwc ! liquid water content [g m-3],  typically .1
      real(ireals), intent(in) :: N   ! droplet density      [1 cm-3], typically 100
      logical, intent(in) :: l_liquid
      ! zkap => breadth parameter e.g. continental=1.143 (Martin et al.), maritime=1.077(Martin et al.)
      real(ireals), intent(in), optional :: zkap

      real(ireals) :: reff            ! effective radius [1e-6 m]

      real(ireals), parameter :: &
        rho = 1e3_ireals, &       ! water density @ 4degC
        zfact = 1e6_ireals * (3e-9_ireals / (4*pi*rho))**(1._ireals/3._ireals) ! conversion factor

      real(ireals) :: zkap_default

      if(l_liquid) then
        zkap_default = 1.1_ireals
        if(present(zkap)) zkap_default = zkap
        reff = zfact * zkap_default * (lwc / max(tiny(N),N) )**(1._ireals/3._ireals)
      else
        reff = 83.8_ireals*lwc**0.216_ireals
      endif
    end function

    ! convert e.g. from lwc [g/kg] to lwp[g/m2] with: lwp = lwc * vert_integral_coeff(p0, p1, T)
    elemental function vert_integral_coeff(p0, p1) result(c)
      real(ireals), intent(in) :: p0, p1 ! pressure at bottom/top of layer [hPa]
      real(ireals) :: c                  ! coeff to convert from [g/kg] to [g m**-2]

      real(ireals) :: dp

      dp    = abs(p1-p0) * 1e2_ireals
      c = dp / grav
    end function

  end module
