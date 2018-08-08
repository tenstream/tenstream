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

  use iso_fortran_env, only: REAL32, REAL64
  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb

  use m_data_parameters, only : iintegers, mpiint, ireals, default_str_len, &
    zero, i2, i9

  use m_helper_functions, only: CHKERR, search_sorted_bisection, reverse, &
    imp_allreduce_min, imp_allreduce_max, meanvec, imp_bcast, read_ascii_file_2d, &
    gradient

  use m_tenstream_interpolation, only : interp_1d

  implicit none

  logical,parameter :: ldebug=.True.
  !logical,parameter :: ldebug=.False.

  interface
    real function PLKINT(WVLLO, WVLHI, T)
      real :: WVLLO, WVLHI, T
    end function
  end interface

  interface hydrostat_dz
    module procedure hydrostat_dz_real32, hydrostat_dz_real64
  end interface

  type t_bg_atm
    real(ireals),allocatable :: plev   (:) ! dim(nlay+1 of background profile)
    real(ireals),allocatable :: tlev   (:) !
    real(ireals),allocatable :: zt     (:) !
    real(ireals),allocatable :: h2o_lev(:) !
    real(ireals),allocatable :: o3_lev (:) !
    real(ireals),allocatable :: co2_lev(:) !
    real(ireals),allocatable :: ch4_lev(:) !
    real(ireals),allocatable :: n2o_lev(:) !
    real(ireals),allocatable :: o2_lev (:) !

    real(ireals),allocatable :: play   (:) ! dim(nlay of background profile)
    real(ireals),allocatable :: zm     (:) !
    real(ireals),allocatable :: dz     (:) !
    real(ireals),allocatable :: tlay   (:) !
    real(ireals),allocatable :: h2o_lay(:) !
    real(ireals),allocatable :: o3_lay (:) !
    real(ireals),allocatable :: co2_lay(:) !
    real(ireals),allocatable :: ch4_lay(:) !
    real(ireals),allocatable :: n2o_lay(:) !
    real(ireals),allocatable :: o2_lay (:) !
  end type

  type t_tenstr_atm
    type(t_bg_atm),allocatable :: bg_atm
    real(ireals),allocatable :: plev   (:,:) ! dim(nlay+1 of merged grid, ncol)
    real(ireals),allocatable :: tlev   (:,:) !
    real(ireals),allocatable :: zt     (:,:) !
    real(ireals),allocatable :: h2o_lev(:,:) !
    real(ireals),allocatable :: o3_lev (:,:) !
    real(ireals),allocatable :: co2_lev(:,:) !
    real(ireals),allocatable :: ch4_lev(:,:) !
    real(ireals),allocatable :: n2o_lev(:,:) !
    real(ireals),allocatable :: o2_lev (:,:) !

    real(ireals),allocatable :: play   (:,:) ! dim(nlay of merged grid, ncol)
    real(ireals),allocatable :: zm     (:,:) !
    real(ireals),allocatable :: dz     (:,:) !
    real(ireals),allocatable :: tlay   (:,:) !
    real(ireals),allocatable :: h2o_lay(:,:) !
    real(ireals),allocatable :: o3_lay (:,:) !
    real(ireals),allocatable :: co2_lay(:,:) !
    real(ireals),allocatable :: ch4_lay(:,:) !
    real(ireals),allocatable :: n2o_lay(:,:) !
    real(ireals),allocatable :: o2_lay (:,:) !
    real(ireals),allocatable :: lwc    (:,:)
    real(ireals),allocatable :: reliq  (:,:)
    real(ireals),allocatable :: iwc    (:,:)
    real(ireals),allocatable :: reice  (:,:)

    logical :: lTOA_to_srfc
    ! index of lowermost layer in atm: search for level where height is bigger and pressure is lower:
    integer(iintegers), allocatable :: atm_ke ! number of vertical levels of background profile, i.e. 1:atm_ke will be put on top of the dynamics grid
  end type

  contains

  subroutine setup_tenstr_atm(comm, lTOA_to_srfc, atm_filename, d_plev, d_tlev, atm, &
      d_tlay, d_h2ovmr, d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr,  d_o2vmr, &
      d_lwc, d_reliq, d_iwc, d_reice)
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
    real(ireals),intent(in),optional :: d_tlay   (:,:) ! layer mean temperature [K]
    real(ireals),intent(in),optional :: d_h2ovmr (:,:) ! watervapor volume mixing ratio [e.g. 1e-3]
    real(ireals),intent(in),optional :: d_o3vmr  (:,:) ! ozone volume mixing ratio      [e.g. .1e-6]
    real(ireals),intent(in),optional :: d_co2vmr (:,:) ! CO2 volume mixing ratio        [e.g. 407e-6]
    real(ireals),intent(in),optional :: d_ch4vmr (:,:) ! methane volume mixing ratio    [e.g. 2e-6]
    real(ireals),intent(in),optional :: d_n2ovmr (:,:) ! n2o volume mixing ratio        [e.g. .32]
    real(ireals),intent(in),optional :: d_o2vmr  (:,:) ! oxygen volume mixing ratio     [e.g. .2]
    real(ireals),intent(in),optional :: d_lwc    (:,:) ! liq water content              [g/kg]
    real(ireals),intent(in),optional :: d_reliq  (:,:) ! effective radius               [micron]
    real(ireals),intent(in),optional :: d_iwc    (:,:) ! ice water content              [g/kg]
    real(ireals),intent(in),optional :: d_reice  (:,:) ! ice effective radius           [micron]

    integer(iintegers) :: i

    if(lTOA_to_srfc) call CHKERR(1_mpiint, 'currently not possible to supply dynamics input starting at the TOP, input should be starting at the surface')

    if(.not.allocated(atm%bg_atm)) then
      call load_atmfile(comm, atm_filename, atm%bg_atm)
      call sanitize_input(atm%bg_atm%plev, atm%bg_atm%tlev, atm%bg_atm%tlay)
    endif

    do i=lbound(d_plev,2),ubound(d_plev,2)
      if(present(d_tlay)) then
        call sanitize_input(d_plev(:,i), d_tlev(:,i), d_tlay(:,i))
      else
        call sanitize_input(d_plev(:,i), d_tlev(:,i))
      endif
    enddo

    call merge_dyn_rad_grid(comm, atm, &
      d_plev, d_tlev, d_tlay, d_h2ovmr, &
      d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr, &
      d_o2vmr, d_lwc, d_reliq, d_iwc, d_reice )

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
    if(allocated(atm%h2o_lev)) atm%h2o_lev = reverse(atm%h2o_lev)
    if(allocated(atm%o3_lev )) atm%o3_lev  = reverse(atm%o3_lev )
    if(allocated(atm%co2_lev)) atm%co2_lev = reverse(atm%co2_lev)
    if(allocated(atm%ch4_lev)) atm%ch4_lev = reverse(atm%ch4_lev)
    if(allocated(atm%n2o_lev)) atm%n2o_lev = reverse(atm%n2o_lev)
    if(allocated(atm%o2_lev )) atm%o2_lev  = reverse(atm%o2_lev )

    if(allocated(atm%play   )) atm%play    = reverse(atm%play   )
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

    atm%lTOA_to_srfc = lTOA_to_srfc
  end subroutine

  !> Concatenate the dynamical grid with the background profile
  subroutine merge_dyn_rad_grid(comm, atm,   &
      in_d_plev, d_tlev, d_tlay, d_h2ovmr,   &
      d_o3vmr, d_co2vmr, d_ch4vmr, d_n2ovmr, &
      d_o2vmr, d_lwc, d_reliq, d_iwc, d_reice )

    integer(mpiint), intent(in) :: comm
    type(t_tenstr_atm),intent(inout) :: atm

    real(ireals),intent(in) :: in_d_plev (:,:), d_tlev(:,:) ! dim(nlay_dynamics+1, ncol)

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

    real(ireals) :: d_plev (ubound(in_d_plev,1), ubound(in_d_plev,2))

    integer(iintegers) :: d_ke, d_ke1 ! number of vertical levels of dynamics grid

    integer(iintegers) :: ke, ke1 ! number of vertical levels of merged grid
    integer(iintegers) :: is, ie, icol, l,m

    real(ireals),allocatable :: d_hhl(:,:), d_dz(:)
    real(ireals) :: global_maxheight, global_minplev

    if(.not.allocated(atm%bg_atm)) call CHKERR(1_mpiint, 'bg_atm has to be allocated before merging with dynamics grid variables')
    associate( bg_atm => atm%bg_atm )

    is = lbound(d_plev,2); ie = ubound(d_plev,2)

    d_ke1 = ubound(d_plev,1); d_ke = d_ke1-1

    ! find out how many layers we have to put on top of the dynamics grid
    if(.not.allocated(atm%atm_ke)) then
      ! first put a tiny increment on the top of dynamics pressure value,
      ! to handle a cornercase where dynamics and background profile are the same
      d_plev = in_d_plev
      !d_plev(d_ke1, :) = d_plev(d_ke1, :) - 1e-3_ireals

      ! First get top height of dynamics grid
      allocate(d_hhl(d_ke1, ie))
      allocate(d_dz(d_ke))

      do icol = is, ie
        if(present(d_tlay)) then
          call hydrostat_lev(d_plev(:,icol),d_tlay(:,icol), zero, d_hhl(:, icol), d_dz)
        else
          call hydrostat_lev(d_plev(:,icol),(d_tlev(1:d_ke,icol)+d_tlev(2:d_ke1,icol))/2, zero, d_hhl(:, icol), d_dz)
        endif
      enddo

      ! index of lowermost layer in atm: search for level where height is bigger and
      ! pressure is lower
      call imp_allreduce_max(comm, maxval(d_hhl), global_maxheight)
      call imp_allreduce_min(comm, minval(d_plev), global_minplev)

      l = floor(search_sorted_bisection(bg_atm%zt, global_maxheight))
      m = floor(search_sorted_bisection(bg_atm%plev, global_minplev))
      allocate(atm%atm_ke)
      atm%atm_ke = min(l,m)
    endif

    associate(atm_ke => atm%atm_ke)
    ke  = atm_ke + d_ke
    ke1 = atm_ke + d_ke1

    ! then from there on couple background atm data on top of that
    if(.not.allocated(atm%plev   )) allocate(atm%plev   (ke1, ie))
    if(.not.allocated(atm%tlev   )) allocate(atm%tlev   (ke1, ie))
    if(.not.allocated(atm%zt     )) allocate(atm%zt     (ke1, ie))
    if(.not.allocated(atm%tlay   )) allocate(atm%tlay   (ke,  ie))
    if(.not.allocated(atm%dz     )) allocate(atm%dz     (ke1, ie))
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

    do icol=is,ie

      ! First merge pressure levels .. pressure is always given..
      atm%plev(ke1-atm_ke+1:ke1, icol) = reverse(bg_atm%plev(1:atm_ke))
      atm%plev(1:d_ke1, icol) = d_plev(:, icol)
      if(atm%plev(ke1-atm_ke+1, icol) .gt. atm%plev(d_ke1, icol)) then
        print *,'background profile pressure is .ge. than uppermost pressure &
          & level of dynamics grid -- this suggests the dynamics grid is way &
          & off hydrostatic balance... please check', atm%plev(:, icol)
        call CHKERR(1_mpiint, 'error in rrtm_lw merging grids')
      endif

      ! And also Tlev has to be present always
      atm%tlev(ke1-atm_ke+1:ke1, icol) = reverse(bg_atm%tlev(1:atm_ke))
      atm%tlev(1:d_ke1, icol) = d_tlev(:,icol)

      if(present(d_tlay)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%tlay, bg_atm%tlev, atm%tlay(:, icol), d_tlay(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%tlay, bg_atm%tlev, atm%tlay(:, icol), &
          (d_tlev(1:d_ke,icol)+d_tlev(2:d_ke1,icol))/2)
      endif

      ! compute dz
      call hydrostat_lev(atm%plev(:,icol),atm%tlay(:,icol), zero, atm%zt(:, icol), atm%dz(:, icol))

      if(present(d_lwc)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, zero*bg_atm%tlay, zero*bg_atm%tlev, atm%lwc(:,icol), d_lwc(:,icol))
      else
        atm%lwc(:,icol) = zero
      endif
      if(present(d_reliq)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, zero*bg_atm%tlay, zero*bg_atm%tlev, atm%reliq(:,icol), d_reliq(:,icol))
      else
        atm%reliq = zero
      endif

      if(present(d_iwc)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, zero*bg_atm%tlay, zero*bg_atm%tlev, atm%iwc(:,icol), d_iwc(:,icol))
      else
        atm%iwc(:,icol) = zero
      endif
      if(present(d_reice)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, zero*bg_atm%tlay, zero*bg_atm%tlev, atm%reice(:,icol), d_reice(:,icol))
      else
        atm%reice = zero
      endif

      if(present(d_h2ovmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%h2o_lay, bg_atm%h2o_lev, atm%h2o_lay(:,icol), d_h2ovmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%h2o_lay, bg_atm%h2o_lev, atm%h2o_lay(:,icol))
      endif
      if(present(d_o3vmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%o3_lay, bg_atm%o3_lev, atm%o3_lay(:,icol), d_o3vmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%o3_lay, bg_atm%o3_lev, atm%o3_lay(:,icol))
      endif

      if(present(d_co2vmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%co2_lay, bg_atm%co2_lev, atm%co2_lay(:,icol), d_co2vmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%co2_lay, bg_atm%co2_lev, atm%co2_lay(:,icol))
      endif
      if(present(d_ch4vmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%ch4_lay, bg_atm%ch4_lev, atm%ch4_lay(:,icol), d_ch4vmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%ch4_lay, bg_atm%ch4_lev, atm%ch4_lay(:,icol))
      endif
      if(present(d_n2ovmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%n2o_lay, bg_atm%n2o_lev, atm%n2o_lay(:,icol), d_n2ovmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%n2o_lay, bg_atm%n2o_lev, atm%n2o_lay(:,icol))
      endif
      if(present(d_o2vmr)) then
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%o2_lay, bg_atm%o2_lev, atm%o2_lay(:,icol), d_o2vmr(:,icol))
      else
        call merge_grid_var(bg_atm%zt, d_hhl(:,icol), atm_ke, bg_atm%o2_lay, bg_atm%o2_lev, atm%o2_lay(:,icol))
      endif
    enddo
  end associate
  end associate
  atm%lTOA_to_srfc = .False.
  end subroutine

  ! merge the dynamics grid and the background profile together at lvl atm_ke
  ! NOTE! Only use with variables on layer
  subroutine merge_grid_var(a_hhl, d_hhl, atm_ke, a_lay, a_lev, col_var, d_var)
    integer(iintegers),intent(in) :: atm_ke
    real(ireals),intent(in) :: a_hhl(:), d_hhl(:), a_lay(:), a_lev(:) ! a_arr is from atm%, d_arr corresponds to dynamics grids
    real(rb),intent(out) :: col_var(:)
    real(ireals),intent(in),optional :: d_var(:)
    integer(iintegers) :: k, kt ! kt is reverse index
    real(ireals) :: h

    ! Top of atmosphere layers are always given by background profile
    do k=1,atm_ke
      kt = size(col_var)-k+1
      col_var(kt) = a_lay(k)
    enddo

    if(present(d_var)) then ! dynamics grid variable is provided, use that
      do k=1,size(d_var)
        col_var(k) = d_var(k)
      enddo
    else ! we may still use atmospheric grid file instead...
      do k=1,size(col_var)-atm_ke
        h = (d_hhl(k+1) + d_hhl(k)) / 2
        col_var(k) = interp_1d(search_sorted_bisection(a_hhl, h), a_lev)
      enddo
    endif
  end subroutine

  pure elemental function hydrostat_dz_real32(dp, p, T)
    real(REAL32), intent(in) :: dp, p, T
    real(REAL32) :: hydrostat_dz_real32, rho
    rho = p / 287.058_REAL32 / T
    hydrostat_dz_real32 = dp / rho / 9.8065_REAL32
  end function
  pure elemental function hydrostat_dz_real64(dp, p, T)
    real(REAL64), intent(in) :: dp, p, T
    real(REAL64) :: hydrostat_dz_real64, rho
    rho = p / 287.058_REAL64 / T
    hydrostat_dz_real64 = dp / rho / 9.8065_REAL64
  end function

  subroutine hydrostat_lev(plev,tlay, hsrfc, hhl, dz)
    ! Integrate vertical height profile hydrostatically. arrays start at bottom(surface)
    real(ireals),intent(in) :: plev(:),tlay(:)
    real(ireals),intent(in) :: hsrfc
    real(ireals),intent(out) :: hhl(size(plev))
    real(ireals),intent(out) :: dz(size(tlay))
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
      stop 'error in dz'
    endif
  end subroutine

  subroutine sanitize_input(plev, tlev, tlay)
    real(ireals),intent(in),dimension(:) :: plev, tlev
    real(ireals),intent(in),dimension(:),optional :: tlay

    integer(mpiint) :: errcnt
    integer(iintegers) :: k
    logical :: lerr

    errcnt = 0
    lerr = maxval(plev) .gt. 1050
    if(lerr) then
      print *,'Pressure above 1050 hPa -- are you sure this is earth?', maxval(plev)
      errcnt = errcnt+1
    endif

    lerr = minval(plev) .lt. zero
    if(lerr) then
      print *,'Pressure negative -- are you sure this is physically correct?', minval(plev)
      errcnt = errcnt+1
    endif

    lerr = minval(tlev) .lt. 180
    if(lerr) then
      print *,'Temperature is very low -- are you sure RRTMG can handle that?', minval(tlev)
      errcnt = errcnt+1
    endif

    lerr = maxval(tlev) .gt. 400
    if(lerr) then
      print *,'Temperature is very high -- are you sure RRTMG can handle that?', maxval(tlev)
      errcnt = errcnt+1
    endif

    if(present(tlay) .and. ldebug) then
      do k=lbound(tlay,1), ubound(tlay,1)
        lerr = (tlay(k)-tlev(k).ge.zero) .eqv. (tlay(k)-tlev(k+1).gt.zero) ! different sign says its in between
        if(lerr) then
          print *,'Layer Temperature not between level temps?', k, tlev(k), '|', tlay(k), '|', tlev(k+1)
          errcnt = errcnt+1
        endif
      enddo
    endif

    if(errcnt.gt.0) then
      print *,'Found wonky input to pprts_rrtm_lw -- please check! -- will abort now.'
      call CHKERR(errcnt)
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
      call read_ascii_file_2d(atm_filename, prof, i9, i2, ierr)
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
          print *,k,'zt', atm%zt(k), 'plev', atm%plev(k), 'T', atm%tlev(k), 'CO2', atm%co2_lev(k), 'H2O', atm%h2o_lev(k), 'O3', atm%o3_lev(k),'N2O' , atm%n2o_lev(k), 'O2', atm%o2_lev(k)
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

end module