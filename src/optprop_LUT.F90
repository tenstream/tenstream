!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
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

module m_optprop_LUT
  use iso_fortran_env, only: INT64

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only: lLUT_mockup

  use m_helper_functions, only : approx,  &
    rel_approx, imp_bcast,                &
    get_arg, ftoa, itoa,                  &
    cstr, char_arr_to_str,                &
    mpi_logical_and, mpi_logical_or,      &
    CHKERR,                               &
    triangle_area_by_vertices,            &
    rad2deg, deg2rad,                     &
    ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets, &
    linspace

  use m_search, only: find_real_location

  use m_data_parameters, only : iintegers, mpiint,           &
    ireals, ireal_dp, irealLUT, ireal_params,                &
    one, zero, i0, i1, i2, i3, i10, nil, inil,               &
    imp_iinteger, imp_ireals, imp_irealLUT, imp_logical,     &
    default_str_len

  use m_optprop_parameters, only:         &
    luse_memory_map,                      &
    ldebug_optprop, lut_basename,         &
    LUT_dump_interval, LUT_max_create_jobtime, &
    LUT_MAX_DIM,                          &
    ldelta_scale,delta_scale_truncate,    &
    stddev_atol, stddev_rtol,             &
    wedge_sphere_radius,                  &
    preset_g6,                            &
    preset_param_phi6,                    &
    preset_param_phi11,                   &
    preset_param_phi19,                   &
    preset_param_theta4,                  &
    preset_param_theta13,                 &
    preset_aspect5,                       &
    preset_aspect7,                       &
    preset_aspect11,                      &
    preset_aspect13,                      &
    preset_aspect17,                      &
    preset_aspect18,                      &
    preset_aspect23,                      &
    preset_w010,                          &
    preset_w020,                          &
    preset_w021,                          &
    preset_tau15,                         &
    preset_tau20,                         &
    preset_tau31

  use m_boxmc, only: t_boxmc, &
    t_boxmc_1_2, &
    t_boxmc_3_6, t_boxmc_3_10, t_boxmc_3_16, &
    t_boxmc_8_10, t_boxmc_8_12, t_boxmc_8_16, &
    t_boxmc_8_18, &
    t_boxmc_wedge_5_8, t_boxmc_wedge_18_8
  use m_tenstream_interpolation, only: interp_4d, interp_vec_simplex_nd
  use m_netcdfio

  use m_mmap, only : arr_to_mmap, munmap_mmap_ptr

  use m_boxmc_geometry, only : setup_default_wedge_geometry, setup_default_unit_cube_geometry

  use m_LUT_param_phi, only: param_phi_from_azimuth, azimuth_from_param_phi, &
    iterative_phi_theta_from_param_phi_and_param_theta, &
    param_phi_param_theta_from_phi_and_theta_withcoords, &
    LUT_wedge_dz

  implicit none

  private
  public :: t_optprop_LUT, &
    t_optprop_LUT_1_2, &
    t_optprop_LUT_3_6, &
    t_optprop_LUT_3_10, &
    t_optprop_LUT_3_16, &
    t_optprop_LUT_8_10, &
    t_optprop_LUT_8_12, &
    t_optprop_LUT_8_16, &
    t_optprop_LUT_8_18, &
    t_optprop_LUT_wedge_5_8, &
    t_optprop_LUT_rectilinear_wedge_5_8, &
    t_optprop_LUT_wedge_18_8, &
    find_lut_dim_by_name, t_LUT_config, &
    azimuth_from_param_phi, param_phi_from_azimuth
  ! This module loads and generates the LUT-tables for Tenstream Radiation
  ! computations.
  ! It also holds functions for interpolation on the regular LUT grid.

  integer(mpiint) :: iierr
  integer(mpiint) :: mpierr

  type t_LUT_dim
    integer(iintegers) :: N      ! size of dimension
    character(len=default_str_len) :: dimname
    real(irealLUT) :: vrange(2) ! min / max of dimension
    real(irealLUT), allocatable :: v(:) ! sampling points of dimension, size N
  end type

  type t_LUT_config
    type(t_LUT_dim), allocatable :: dims(:)
    integer(iintegers),allocatable :: offsets(:) ! offsets of respective dimensions (starts at 0 and next one is dim(1)%N ... etc... )
  end type

  type t_table
    real(irealLUT), contiguous, pointer :: c(:,:) => NULL() ! depending on config has Ndim_1*Ndim_2*etc. many entries
    real(irealLUT), allocatable :: stddev_tol(:) ! maxval of tolerance for a given entry
    character(default_str_len), allocatable :: table_name_c(:)
    character(default_str_len), allocatable :: table_name_tol(:)
  end type

  type,abstract :: t_optprop_LUT
    class(t_boxmc), allocatable :: bmc
    type(t_table), allocatable :: Sdiff, Sdir, Tdir
    type(t_LUT_config), allocatable :: dirconfig, diffconfig
    integer(iintegers) :: dir_streams = inil, diff_streams = inil
    logical :: LUT_initialized=.False., optprop_LUT_debug=ldebug_optprop
    character(default_str_len) :: lutbasename

    contains
      procedure :: init
      procedure :: destroy
      procedure :: LUT_get_dir2dir
      procedure :: LUT_get_dir2diff
      procedure :: LUT_get_diff2diff
      procedure :: LUT_bmc_wrapper_determine_sample_pts
      procedure :: LUT_bmc_wrapper
      procedure :: bmc_wrapper
      procedure :: scatter_LUTtables
      procedure :: createLUT
      procedure :: loadLUT_dir
      procedure :: loadLUT_diff
      procedure :: set_parameter_space
      procedure :: print_configs
  end type

  type,extends(t_optprop_LUT) :: t_optprop_LUT_1_2
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_3_6
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_3_10
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_3_16
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_8_10
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_8_12
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_8_16
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_8_18
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_wedge_5_8
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_rectilinear_wedge_5_8
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_wedge_18_8
  end type

  ! You should not need to change this... but feel free to play around...
  ! interp_mode 1 == nearest neighbour interpolation
  ! interp_mode 2 == linear interpolation
  integer(iintegers), parameter :: interp_mode=2


  logical, parameter :: ldebug=.False.

contains

  subroutine init(OPP, comm, skip_load_LUT)
      class(t_optprop_LUT) :: OPP
      integer(mpiint) ,intent(in) :: comm
      logical, intent(in), optional :: skip_load_LUT

      integer(mpiint) :: comm_size, myid, ierr
      logical :: lskip_load_LUT, load_diffuse_LUT_first, lflg

      if(OPP%LUT_initialized) return

      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Initializing LUT`s...'

      if(.not.allocated(OPP%bmc)) then
        select type (OPP)
          class is (t_optprop_LUT_1_2)
            OPP%dir_streams  =  1
            OPP%diff_streams =  2
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_1_2::OPP%bmc)

          class is (t_optprop_LUT_3_6)
            OPP%dir_streams  = 3
            OPP%diff_streams = 6
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_3_6::OPP%bmc)

          class is (t_optprop_LUT_3_10)
            OPP%dir_streams  = 3
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_3_10::OPP%bmc)

          class is (t_optprop_LUT_3_16)
            OPP%dir_streams  = 3
            OPP%diff_streams = 16
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_3_16::OPP%bmc)

          class is (t_optprop_LUT_8_10)
            OPP%dir_streams  = 8
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_8_10::OPP%bmc)

          class is (t_optprop_LUT_8_12)
            OPP%dir_streams  = 8
            OPP%diff_streams = 12
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_8_12::OPP%bmc)

          class is (t_optprop_LUT_8_16)
            OPP%dir_streams  = 8
            OPP%diff_streams = 16
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_8_16::OPP%bmc)

          class is (t_optprop_LUT_8_18)
            OPP%dir_streams  = 8
            OPP%diff_streams = 18
            OPP%lutbasename=trim(lut_basename)
            allocate(t_boxmc_8_18::OPP%bmc)

          class is (t_optprop_LUT_wedge_5_8)
            OPP%dir_streams  = 5
            OPP%diff_streams = 8
            OPP%lutbasename=trim(lut_basename)//'_wedge.Rsphere'//itoa(int(wedge_sphere_radius))
            allocate(t_boxmc_wedge_5_8::OPP%bmc)

          class is (t_optprop_LUT_rectilinear_wedge_5_8)
            OPP%dir_streams  = 5
            OPP%diff_streams = 8
            OPP%lutbasename=trim(lut_basename)//'_rectilinear_wedge'
            allocate(t_boxmc_wedge_5_8::OPP%bmc)

          class is (t_optprop_LUT_wedge_18_8)
            OPP%dir_streams  = 18
            OPP%diff_streams = 8
            OPP%lutbasename=trim(lut_basename)//'_wedge.Rsphere'//itoa(int(wedge_sphere_radius))
            allocate(t_boxmc_wedge_18_8::OPP%bmc)

          class default
            stop 'initialize LUT: unexpected type for optprop_LUT object!'
        end select

        call OPP%bmc%init(comm)
      endif

      call OPP%set_parameter_space()

      load_diffuse_LUT_first = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        '-load_diffuse_LUT_first', load_diffuse_LUT_first, lflg, ierr); call CHKERR(ierr)

      lskip_load_LUT = get_arg(.True., skip_load_LUT)
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        '-skip_load_LUT', lskip_load_LUT, lflg, ierr); call CHKERR(ierr)
      if(.not.lskip_load_LUT .and. myid.eq.0) &
        print *,'loading and checking LUT`s from netCDF', .not.lskip_load_LUT

      if(load_diffuse_LUT_first) then
        call OPP%loadLUT_diff(comm, lskip_load_LUT)
        call OPP%loadLUT_dir(comm, lskip_load_LUT)
      else
        call OPP%loadLUT_dir(comm, lskip_load_LUT)
        call OPP%loadLUT_diff(comm, lskip_load_LUT)
      endif

      call OPP%scatter_LUTtables(comm)

      OPP%LUT_initialized=.True.
      if(ldebug.and.myid.eq.0) call print_configs(OPP)
      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Initializing LUT`s... finished'
  end subroutine

  subroutine destroy(OPP)
      class(t_optprop_LUT) :: OPP
      integer(mpiint) :: ierr
      if(allocated(OPP%Tdir)) then
        if(luse_memory_map) then
          call munmap_mmap_ptr(OPP%Tdir%c, ierr); call CHKERR(ierr)
        endif
        deallocate(OPP%Tdir)
      endif
      if(allocated(OPP%Sdir )) then
        if(luse_memory_map) then
          call munmap_mmap_ptr(OPP%Sdir%c, ierr); call CHKERR(ierr)
        endif
        deallocate(OPP%Sdir)
      endif
      if(allocated(OPP%Sdiff)) then
        if(luse_memory_map) then
          call munmap_mmap_ptr(OPP%Sdiff%c, ierr); call CHKERR(ierr)
        endif
        deallocate(OPP%Sdiff)
      endif
      if(allocated(OPP%bmc  )) deallocate(OPP%bmc)
      if(allocated(OPP%dirconfig)) deallocate(OPP%dirconfig)
      if(allocated(OPP%diffconfig)) deallocate(OPP%diffconfig)
      OPP%LUT_initialized=.False.
      if(OPP%optprop_LUT_debug) print *,'Destroyed LUTs', OPP%LUT_initialized
  end subroutine

  function gen_lut_basename(prefix, config) result(lutname)
    character(len=default_str_len) :: lutname
    character(len=*), intent(in) :: prefix
    type(t_lut_config), intent(in) :: config
    integer(iintegers) :: k
    lutname = trim(prefix)
    do k=1,size(config%dims)
      lutname = trim(lutname)//'.'//trim(config%dims(k)%dimname)//itoa(config%dims(k)%N)
    enddo
    lutname = trim(lutname)//'.ds'//itoa(int(delta_scale_truncate*1000))
  end function

subroutine load_table_from_netcdf(table, istat)
  type(t_table), intent(inout) :: table
  integer(mpiint), intent(out) :: istat

  integer(iintegers) :: k
  integer(mpiint) :: ierr
  istat = 0

  call ncload(table%table_name_tol, table%stddev_tol, ierr); istat = istat + ierr
  if(ierr.eq.0) then
    if(ldebug) print *,'we were able to load stddev from file '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_tol,'/'), 'green')//new_line('')// &
      ' but still have to check if they all have a good enough precision / noise level...'
    do k = 1, size(table%stddev_tol, kind=iintegers)
      if(table%stddev_tol(k).gt.stddev_atol+10*epsilon(stddev_atol) &
        .or. table%stddev_tol(k).lt.0._irealLUT) then
        istat = istat + 1
        if(istat.ne.0) then
          print *,'coefficients stddev not good enough', k, table%stddev_tol(k),&
            'should be positive and less than', stddev_atol+10*epsilon(stddev_atol)
          exit
        endif
      endif
    enddo
  else
    print *,'loading stddev_tolerances from '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_tol,'/'), 'purple')//new_line('')// &
      ' failed', istat, ierr
  endif
  if(istat.ne.0) &
    print *,'Test if coeffs in '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_tol,'/'), 'purple')//new_line('')// &
      ' are good results in:', istat

  call ncload(table%table_name_c, table%c, ierr); istat = istat + ierr
  if(ierr.eq.0) then ! we were able to load coeffs but still have to check if they are all ok...
    do k = 1, size(table%c, dim=2, kind=iintegers)
      if(any(table%c(:,k).gt.one) .or. any(table%c(:,k).lt.zero)) then
        istat = istat + 2
        if(istat.ne.0) then
          print *,'coefficients not good enough', k, table%c(:,k), 'should be between [0,1]'
          exit
        endif
      endif
    enddo
  else
    print *,'loading coeffs from '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_c,'/'), 'purple')//new_line('')// &
      ' failed', istat, ierr
  endif

  if(istat.ne.0) then
    print *,'Test if coeffs in '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_c,'/'), 'purple')//new_line('')// &
      ' are good results in:', istat
  endif
end subroutine

subroutine loadLUT_diff(OPP, comm, skip_load_LUT)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    logical, intent(in) :: skip_load_LUT

    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(allocated(OPP%Sdiff)) return ! already loaded
    allocate(OPP%Sdiff)
    OPP%Sdiff%c => NULL()
    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('_diffuse_'//itoa(OPP%diff_streams), OPP%diffconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'diffuse'

    str(3) = 'Stol'; allocate(OPP%Sdiff%table_name_tol(size(str)), source=str)
    str(3) = 'S'   ; allocate(OPP%Sdiff%table_name_c(size(str)), source=str)

    if(skip_load_LUT) return

    if(myid.eq.0) then
      call load_table_from_netcdf(OPP%Sdiff, iierr); errcnt = errcnt + iierr
      call write_pspace(str(1), OPP%diffconfig)
    endif

    call mpi_bcast(errcnt, 1_mpiint, imp_iinteger, 0_mpiint ,comm ,mpierr); call CHKERR(mpierr)

    if(errcnt.ne.0) then ! something went wrong loading the LUT
      call OPP%createLUT(comm, OPP%diffconfig, OPP%Sdiff)
    endif

    if(allocated(OPP%Sdiff%stddev_tol)) deallocate(OPP%Sdiff%stddev_tol)
    if(ldebug) print *,'Loaded Diff2Diff LUT with shape:', shape(OPP%Sdiff%c)
end subroutine

subroutine loadLUT_dir(OPP, comm, skip_load_LUT)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    logical, intent(in) :: skip_load_LUT

    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(allocated(OPP%Sdir).and.allocated(OPP%Tdir)) return ! already loaded
    allocate(OPP%Sdir)
    allocate(OPP%Tdir)
    OPP%Sdir%c => NULL()
    OPP%Tdir%c => NULL()

    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('_direct_'//itoa(OPP%dir_streams)//'_'//itoa(OPP%diff_streams), OPP%dirconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'direct'

    str(3) = 'Stol'; allocate(OPP%Sdir%table_name_tol(size(str)), source=str)
    str(3) = 'S'   ; allocate(OPP%Sdir%table_name_c(size(str)), source=str)

    str(3) = 'Ttol'; allocate(OPP%Tdir%table_name_tol(size(str)), source=str)
    str(3) = 'T'   ; allocate(OPP%Tdir%table_name_c(size(str)), source=str)

    if(skip_load_LUT) return

    if(myid.eq.0) then
      call load_table_from_netcdf(OPP%Sdir, iierr); errcnt = errcnt + iierr
      call load_table_from_netcdf(OPP%Tdir, iierr); errcnt = errcnt + iierr
      call write_pspace(str(1), OPP%dirconfig)
    endif

    call mpi_bcast(errcnt, 1_mpiint, imp_iinteger, 0_mpiint ,comm ,mpierr); call CHKERR(mpierr)

    if(errcnt.ne.0) then ! something went wrong loading the LUT
      call OPP%createLUT(comm, OPP%dirconfig, OPP%Sdir, OPP%Tdir)
    endif

    if(allocated(OPP%Tdir%stddev_tol)) deallocate(OPP%Tdir%stddev_tol)
    if(allocated(OPP%Sdir%stddev_tol)) deallocate(OPP%Sdir%stddev_tol)
    if(ldebug) print *,'Loaded Dir2Dir LUT with shape:', shape(OPP%Tdir%c)
    if(ldebug) print *,'Loaded Dir2Diff LUT with shape:', shape(OPP%Sdir%c)
end subroutine

subroutine write_pspace(fname, config)
  character(len=*), intent(in) :: fname
  type(t_lut_config), intent(in) :: config
  character(len=default_str_len) :: groups(4)

  integer(iintegers) :: kdim
  integer(mpiint) :: ierr
  real(irealLUT), allocatable :: existing_values(:)

  groups(1) = trim(fname)
  groups(2) = 'pspace'

  do kdim = 1, size(config%dims)
    groups(3) = trim(config%dims(kdim)%dimname)

    ! First try to load em from an existing file and compare
    groups(4) = 'values'
    if(allocated(existing_values)) deallocate(existing_values)
    call ncload(groups, existing_values, ierr)
    if(ierr.eq.0) then
      if(.not.all(approx(existing_values, config%dims(kdim)%v, sqrt(epsilon(1._irealLUT))*10))) then
        print *, kdim, trim(groups(1)), trim(groups(3)), new_line(''), &
          ':existing', existing_values, new_line(''), &
          ':new', config%dims(kdim)%v
        call CHKERR(1_mpiint, 'Dimensions of LUT and in optprop_parameters definition do not match!')
      endif
    else ! Otherwise, just save the current ones
      if(ldebug) print *,kdim,'Writing pspace for ',trim(groups(1)),':', trim(groups(2)), &
        ':',trim(groups(3)), trim(groups(4)), '::', config%dims(kdim)%v
      call ncwrite(groups, config%dims(kdim)%v, ierr); call CHKERR(ierr)
    endif

    groups(4) = 'range'
    if(allocated(existing_values)) deallocate(existing_values)
    call ncload(groups, existing_values, ierr)
    if(ierr.eq.0) then
      if(.not.all(approx(existing_values, config%dims(kdim)%vrange, sqrt(epsilon(1._irealLUT))*10))) then
        call CHKERR(1_mpiint, 'Range of dimensions of LUT and in optprop_parameters definition do not match!')
      endif
    else ! Otherwise, just save the current ones
      if(ldebug) print *,kdim,'Writing pspace for ',trim(groups(1)),':', trim(groups(2)), &
        ':',trim(groups(3)), trim(groups(4)), '::', config%dims(kdim)%vrange
      call ncwrite(groups, config%dims(kdim)%vrange, ierr); call CHKERR(ierr)
    endif
  enddo
end subroutine

subroutine createLUT(OPP, comm, config, S, T)
    class(t_optprop_LUT), intent(in) :: OPP
    integer(mpiint),intent(in) :: comm
    type(t_lut_config),intent(in) :: config
    type(t_table) :: S
    type(t_table), optional :: T

    logical :: gotmsg
    integer(mpiint) :: status(MPI_STATUS_SIZE)

    integer(mpiint), parameter :: READYMSG=1,HAVERESULTSMSG=2, WORKMSG=3, FINALIZEMSG=4, RESULTMSG=5

    integer(iintegers) :: idummy, lutindex
    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(myid.eq.0) then
      call prepare_table_space(OPP, config, S, T)
    endif

    if(myid.le.0 .and. comm_size.le.1) &
      call CHKERR(1_mpiint, 'At the moment creation of direct Lookuptable needs at least two mpi-ranks to work...'// &
        'please run with more ranks.')

    if(myid.eq.0) then
      call master(S, T)
      print *,'done calculating coefficients. present(T)?', present(T)
    else
      call worker(config)
    endif

    contains
      subroutine master(S, T)
        type(t_table),intent(inout) :: S
        type(t_table),intent(inout), optional :: T

        integer(iintegers) :: total_size, cnt, finalizedworkers

        integer(iintegers) :: Nsrc
        integer(iintegers) :: isrc, idst, ind

        logical :: ldoneS, ldoneT

        real(irealLUT),allocatable, dimension(:,:) :: S_diff, T_dir, S_tol, T_tol

        real :: starttime, lastsavetime, now
        integer :: clock_count, clock_count_rate
        call system_clock(clock_count, clock_count_rate)
        starttime = clock_count / clock_count_rate
        lastsavetime = starttime

        finalizedworkers=0
        if(present(T)) then
          Nsrc = OPP%dir_streams
        else
          Nsrc = OPP%diff_streams
        endif
        total_size = size(S%c, dim=2) ! we have Ndims work packages with Nsrc bmc calls each

        allocate(S_diff(OPP%diff_streams, Nsrc), S_tol(OPP%diff_streams, Nsrc), &
          T_dir(OPP%dir_streams, Nsrc), T_tol(OPP%dir_streams, Nsrc))

        cnt = 1
        do ! endless loop until everything has been done

          ! Check if we already calculated the coefficients
          if(cnt.le.total_size) then

            ldoneS = all( [ &
              all(S%c(:,cnt).ge.zero), &
              all(S%c(:,cnt).le.one), &
              S%stddev_tol(cnt).le.stddev_atol, &
              S%stddev_tol(cnt).ge.0._irealLUT ] )
            !print *,'ldoneS', cnt, ldoneS, ':', S%stddev_tol(cnt)

            if(present(T)) then
              ldoneT = all( [ &
                all(T%c(:,cnt).ge.zero), &
                all(T%c(:,cnt).le.one), &
                T%stddev_tol(cnt).le.stddev_atol, &
                T%stddev_tol(cnt).ge.0._irealLUT ] )
            else
              ldoneT = .True.
            endif

            if( ldoneS .and. ldoneT ) then
              if( mod(cnt-1, total_size*1/100).eq.0 ) & !every 1 percent report status
                  print *,'Resuming from LUT... ',cnt/(total_size/100),'%'
              cnt=cnt+1
              cycle
            endif
          endif

          ! Now that we know we got something to do, lets find a suitable worker
          gotmsg=.False.
          call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

          if (gotmsg) then

            select case (status(MPI_TAG))

            case(READYMSG)

              ! capture the READY MSG -- we should not leave messages hanging around.
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), READYMSG, &
                      comm, status, mpierr) ; call CHKERR(mpierr)

              if(cnt.le.total_size) then ! we got something to do for a worker -- send him...
                lutindex = cnt
                call mpi_send(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), WORKMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                call mpi_send(present(T), 1_mpiint, imp_logical, status(MPI_SOURCE), WORKMSG, &
                        comm, mpierr); call CHKERR(mpierr)

              else ! no more work to do... tell the worker to quit
                call mpi_send(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                        comm, mpierr); call CHKERR(mpierr)
              endif
              cnt = cnt+1

            case(HAVERESULTSMSG)
              call mpi_recv(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)

              call mpi_recv(S_diff, size(S_diff, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_tol , size(S_tol , kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_dir , size(T_dir , kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_tol , size(T_tol , kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)

                ! Sort coefficients into destination ordering and put em in LUT
              do isrc = 1, Nsrc
                do idst = 1, OPP%diff_streams
                  ind = (idst-1) * Nsrc + isrc
                  S%c(ind, lutindex) = S_diff(idst, isrc)
                enddo

                if(present(T)) then
                  do idst = 1, OPP%dir_streams
                    ind = (idst-1) * Nsrc + isrc
                    T%c(ind, lutindex) = T_dir(idst, isrc)
                  enddo
                endif
              enddo ! isrc

              S%stddev_tol(lutindex) = maxval(S_tol)
              if(present(T)) then
                T%stddev_tol(lutindex) = maxval(T_tol)
              endif

              if (.False. .and. ldebug) call random_print_coeffs(lutindex, S_diff, T_dir, S_tol, T_tol)

              if( mod(lutindex-1, max(i1, total_size/1000_iintegers)).eq.0 ) & !every .1 percent report status
                print *,'Calculated LUT...', lutindex, &
                        real(lutindex-1, irealLUT)*100._irealLUT/real(total_size, irealLUT),'%'

              call system_clock(clock_count, clock_count_rate)
              now = clock_count / clock_count_rate
              if( (now-lastsavetime).gt.LUT_dump_interval .or. (now-starttime).gt.LUT_max_create_jobtime ) then !every 30 minutes wall clock time, dump the LUT.
                print *,'Dumping LUT after ',(now-lastsavetime)/60,'minutes'
                if(present(T)) then
                  print *,'Writing table to file... ', char_arr_to_str(T%table_name_c)
                  call ncwrite(T%table_name_c, T%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
                  print *,'Writing table to file... ', char_arr_to_str(T%table_name_tol)
                  call ncwrite(T%table_name_tol, T%stddev_tol, iierr); call CHKERR(iierr, 'Could not write Table to file')
                endif
                print *,'Writing table to file... ', char_arr_to_str(S%table_name_c)
                call ncwrite(S%table_name_c, S%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
                print *,'Writing table to file... ', char_arr_to_str(S%table_name_tol)
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr); call CHKERR(iierr, 'Could not write Table to file')
                print *,'done writing!',iierr
                lastsavetime = now ! reset the countdown
                if((now-starttime).gt.LUT_max_create_jobtime) then
                  call CHKERR(int(now-starttime, mpiint),'Maximum duration of create_LUT reached ... exiting')
                endif
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                      comm, status, mpierr); call CHKERR(mpierr)
              finalizedworkers = finalizedworkers+1
              if(finalizedworkers.eq.comm_size-1) exit ! all work is done

            end select
          endif
        enddo

        print *,'Writing table to file... '
        call ncwrite(S%table_name_c  , S%c, iierr)
        call ncwrite(S%table_name_tol, S%stddev_tol, iierr)
        if(present(T)) then
          call ncwrite(T%table_name_c  , T%c, iierr)
          call ncwrite(T%table_name_tol, T%stddev_tol, iierr)
          print *,'done writing!',iierr,':: max_atol S',maxval(S%stddev_tol),'max_atol T',maxval(T%stddev_tol)
        else
          print *,'done writing!',iierr,':: max_atol S',maxval(S%stddev_tol)
        endif
      end subroutine
      subroutine worker(config)
          type(t_lut_config), intent(in) :: config
          integer(iintegers) :: isrc, Nsrc
          logical :: ldir

          real(irealLUT),allocatable, dimension(:,:) :: S_diff, T_dir, S_tol, T_tol

          ! workers send READY message to master
          call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, READYMSG, comm, mpierr); call CHKERR(mpierr)

          do
            ! ask what to do
            gotmsg=.False.
            call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

            if (gotmsg) then

              select case (status(MPI_TAG))

              case(WORKMSG)
                ! wait for work to arrive
                call mpi_recv( lutindex, 1_mpiint, imp_iinteger, 0_mpiint, WORKMSG, &
                        comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( ldir, 1_mpiint, imp_logical, 0_mpiint, WORKMSG, &
                        comm, status, mpierr); call CHKERR(mpierr)

                if(ldir) then
                  Nsrc = OPP%dir_streams
                else
                  Nsrc = OPP%diff_streams
                endif
                allocate(S_diff(OPP%diff_streams, Nsrc), S_tol(OPP%diff_streams, Nsrc), &
                         T_dir(OPP%dir_streams, Nsrc), T_tol(OPP%dir_streams, Nsrc))

                do isrc = 1, Nsrc
                  call OPP%LUT_bmc_wrapper(config, lutindex, isrc, ldir, &
                    mpi_comm_self, S_diff(:,isrc), T_dir(:,isrc), S_tol(:,isrc), T_tol(:,isrc))

                  if ( maxval(S_tol(:,isrc)).gt.stddev_atol+sqrt(epsilon(stddev_atol)) &
                    .or. maxval(T_tol(:,isrc)).gt.stddev_atol+sqrt(epsilon(stddev_atol))) then
                    print *,'Ttol', T_tol(:,isrc)
                    print *,'Stol', S_tol(:,isrc)
                    call CHKERR(1_mpiint, 'BOXMC violated stddev constraints!')
                  endif
                enddo

                !print *,'Computed isrc',isrc,'aspect',aspect_zx,'tau',tau_z, w0, g,':', phi, theta
                !print *,myid,'Computed values for ',lutindex, isrc, ldir

                call mpi_send(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                call mpi_send(S_diff, size(S_diff,kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                call mpi_send(S_tol , size(S_tol ,kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                call mpi_send(T_dir , size(T_dir ,kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                call mpi_send(T_tol , size(T_tol ,kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                        comm, mpierr); call CHKERR(mpierr)

                call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, READYMSG, &
                        comm, mpierr); call CHKERR(mpierr)

                deallocate(S_diff, S_tol, T_dir, T_tol)

              case(FINALIZEMSG)
                call mpi_recv(idummy, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, &
                        comm, status, mpierr); call CHKERR(mpierr)
                call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, &
                        comm, mpierr); call CHKERR(mpierr)
                exit

              end select

            endif !gotmsg
          enddo
      end subroutine
      subroutine random_print_coeffs(lutindex, S_diff, T_dir, S_tol, T_tol)
        integer(iintegers) :: lutindex
        real(irealLUT), dimension(:,:) :: S_diff, T_dir, S_tol, T_tol
        real(ireals) :: R
        real(ireals), parameter :: chance=1e-4
        integer(iintegers) :: isrc
        call random_number(R)
        if(R.lt.chance) then
          do isrc=1,size(T_dir,dim=1)
            print *,'lutindex '//itoa(lutindex)//' src '//itoa(isrc)//' :T', T_dir(:, isrc)
          enddo
          print *,'Min/Max T_tol ', minval(T_tol), maxval(T_tol)
          do isrc=1,size(T_dir,dim=1)
            print *,'lutindex '//itoa(lutindex)//' src '//itoa(isrc)//' :S', S_diff(:, isrc)
          enddo
          print *,'Min/Max S_tol ', minval(S_tol), maxval(S_tol)
        endif
      end subroutine
end subroutine createLUT

subroutine prepare_table_space(OPP, config, S, T)
  class(t_optprop_LUT) :: OPP
  type(t_lut_config), intent(in) :: config
  type(t_table),intent(inout) :: S
  type(t_table),intent(inout), optional :: T

  integer(INT64) :: entries, bytesize
  integer(iintegers) :: errcnt

  if(present(T)) then
    entries = &
      int(OPP%dir_streams**2, int64)               * int(product(config%dims(:)%N), int64) + &
      int(OPP%diff_streams*OPP%dir_streams, int64) * int(product(config%dims(:)%N), int64) + &
      2_int64 * int(product(config%dims(:)%N), int64) ! tolerances
  else
    entries = int(OPP%diff_streams**2+1, int64) * int(product(config%dims(:)%N), int64)
  endif
  bytesize = int(C_SIZEOF(1._irealLUT), INT64) * entries
  print *,'Allocating Space for LUTs '//itoa(entries)// &
    ' entries ( '//ftoa(real(bytesize, irealLUT)/1024._irealLUT**3)//' Gb) ...'

  errcnt = 0
  if(present(T)) then
    if(.not.associated(S%c)) &
      allocate(S%c(OPP%diff_streams*OPP%dir_streams, product(config%dims(:)%N)), source=-1._irealLUT)
    if(.not.associated(T%c)) &
      allocate(T%c(OPP%dir_streams*OPP%dir_streams , product(config%dims(:)%N)), source=-1._irealLUT)

    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(product(config%dims(:)%N)), source=huge(-1._irealLUT))
    if(.not.allocated (T%stddev_tol)) allocate(T%stddev_tol(product(config%dims(:)%N)), source=huge(-1._irealLUT))
  else
    if(.not.associated(S%c)) allocate(S%c(OPP%diff_streams**2, product(config%dims(:)%N)), source=-1._irealLUT)
    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(product(config%dims(:)%N)), source=huge(-1._irealLUT))
  endif
  print *,'Allocating Space for LUTs '//itoa(entries)// &
    ' entries ( '//ftoa(real(bytesize, irealLUT)/1024._irealLUT**3)//' Gb) ... done'
end subroutine

! return the integer in config%dims that corresponds to the given dimension
function find_lut_dim_by_name(config, dimname) result(kdim)
  type(t_lut_config), intent(in) :: config
  character(len=*), intent(in) :: dimname
  integer(iintegers) :: kdim

  integer(iintegers) :: k
  do k=1,size(config%dims)
    if(trim(dimname).eq.trim(config%dims(k)%dimname)) then
      kdim = k
      return
    endif
  enddo
  kdim=-1
end function

subroutine get_sample_pnt_by_name_and_index(config, dimname, index_1d, sample_pnt, ierr)
  type(t_lut_config), intent(in) :: config
  character(len=*), intent(in) :: dimname
  integer(iintegers), intent(in) :: index_1d
  real(irealLUT), intent(inout) :: sample_pnt
  integer(mpiint), intent(out) :: ierr

  integer(iintegers) :: kdim, nd_indices(size(config%dims))

  call ind_1d_to_nd(config%offsets, index_1d, nd_indices)

  kdim = find_lut_dim_by_name(config, trim(dimname))
  if(kdim.lt.i1) then ! could not find the corresponding dimension
    ierr = 1
    return
  endif
  if(nd_indices(kdim).gt.size(config%dims(kdim)%v)) then
    print *,index_1d,'nd_indices', nd_indices
    call CHKERR(1_mpiint, 'wrong indices in kdim')
  endif
  sample_pnt = config%dims(kdim)%v(nd_indices(kdim))
  ierr = 0
end subroutine

subroutine LUT_bmc_wrapper_determine_sample_pts(OPP, config, index_1d, dir, &
    vertices, tauz, w0, g, phi, theta, lvalid_entry)
  class(t_optprop_LUT) :: OPP
  type(t_lut_config), intent(in) :: config
  integer(iintegers), intent(in) :: index_1d
  logical, intent(in) :: dir
  real(ireal_dp), intent(out), allocatable :: vertices(:)
  real(irealLUT), intent(out) :: tauz, w0, g, phi, theta
  logical, intent(out) :: lvalid_entry

  real(irealLUT) :: aspect_zx, aspect_zy
  integer(mpiint) :: ierr

  lvalid_entry = .True.

  call get_sample_pnt_by_name_and_index(config, 'tau', index_1d, tauz, ierr)
  call CHKERR(ierr, 'tauz has to be present')

  call get_sample_pnt_by_name_and_index(config, 'w0', index_1d, w0, ierr)
  call CHKERR(ierr, 'w0 has to be present')

  call get_sample_pnt_by_name_and_index(config, 'g', index_1d, g, ierr)
  if(ierr.ne.0) then
    g = zero ! defaults to isotropic scattering
  endif

  call get_sample_pnt_by_name_and_index(config, 'aspect_zx', index_1d, aspect_zx, ierr)
  call CHKERR(ierr, 'aspect_zx has to be present')

  call get_sample_pnt_by_name_and_index(config, 'aspect_zy', index_1d, aspect_zy, ierr)
  if(ierr.ne.0) then
    aspect_zy = aspect_zx ! set dy = dx
  endif

  ! default values for angles in case of diff2diff and such
  phi = 0
  theta = -1

  ! First define Default vertices for Cube Geometry
  select type(OPP)
  class is (t_optprop_LUT_1_2)
    call prep_pprts()

  class is (t_optprop_LUT_3_6)
    call prep_pprts()

  class is (t_optprop_LUT_3_10)
    call prep_pprts()

  class is (t_optprop_LUT_3_16)
    call prep_pprts()

  class is (t_optprop_LUT_8_10)
    call prep_pprts()

  class is (t_optprop_LUT_8_12)
    call prep_pprts()

  class is (t_optprop_LUT_8_16)
    call prep_pprts()

  class is (t_optprop_LUT_8_18)
    call prep_pprts()

  class is (t_optprop_LUT_wedge_5_8)
    call prep_plexrt(sphere_radius=real(wedge_sphere_radius, ireal_dp))

  class is (t_optprop_LUT_rectilinear_wedge_5_8)
    call prep_plexrt()

  class is (t_optprop_LUT_wedge_18_8)
    call prep_plexrt(sphere_radius=real(wedge_sphere_radius, ireal_dp))

  class default
    call CHKERR(1_mpiint, 'unexpected type for optprop_LUT object!')
end select
contains
  subroutine prep_pprts()
    call setup_default_unit_cube_geometry(1._ireal_dp, &
      real(aspect_zx/aspect_zy, ireal_dp), real(aspect_zx, ireal_dp), &
      vertices)
    if(dir) then
      call get_sample_pnt_by_name_and_index(config, 'phi', index_1d, phi, ierr)
      call CHKERR(ierr, 'phi has to be present for direct calculations')
      call get_sample_pnt_by_name_and_index(config, 'theta', index_1d, theta, ierr)
      call CHKERR(ierr, 'theta has to be present for direct calculations')
    endif
  end subroutine
  subroutine prep_plexrt(sphere_radius)
    real(ireal_dp), intent(in), optional :: sphere_radius
    real(irealLUT) :: param_phi, param_theta
    real(irealLUT), dimension(2) :: wedge_C
    real(ireal_params) :: rphi, rtheta
    real(irealLUT) :: Atop, dz

    wedge_C(1) = .5_irealLUT
    call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cx', index_1d, wedge_C(1), ierr)
    call CHKERR(ierr, 'wedge_coord_Cx has to be present for wedge calculations')

    wedge_C(2) = 0.8660254_irealLUT
    call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cy', index_1d, wedge_C(2), ierr)
    call CHKERR(ierr, 'wedge_coord_Cy has to be present for wedge calculations')


    Atop = triangle_area_by_vertices([0._irealLUT, 0._irealLUT], [1._irealLUT, 0._irealLUT], wedge_C)
    dz = LUT_wedge_dz(Atop, aspect_zx)

    call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], &
                                      [1._ireal_dp, 0._ireal_dp], &
                                      real(wedge_C, ireal_dp), &
                                      real(dz, ireal_dp), &
                                      vertices, sphere_radius=sphere_radius)

    if(dir) then
      call get_sample_pnt_by_name_and_index(config, 'param_phi', index_1d, param_phi, ierr)
      call CHKERR(ierr, 'param_phi has to be present for wedge calculations')
      call get_sample_pnt_by_name_and_index(config, 'param_theta', index_1d, param_theta, ierr)
      call CHKERR(ierr, 'param_theta has to be present for wedge calculations')
      call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
        real(param_phi, ireal_params), real(param_theta, ireal_params), rphi, rtheta, ierr); call CHKERR(ierr)
      phi   = real(rad2deg(rphi), irealLUT)
      theta = real(rad2deg(rtheta), irealLUT)

      if(param_theta.le.0._ireals) then
        if(param_phi.lt.-1._ireals .or. param_phi.gt.1._ireals) then
          lvalid_entry = .False. ! this is a dead zone in the param space. Should never happen in a mesh
        endif
      endif
    endif
  end subroutine
end subroutine

subroutine LUT_bmc_wrapper_validate(OPP, config, tauz, w0, index_1d, src, dir, T, S, ierr)
  class(t_optprop_LUT) :: OPP
  type(t_lut_config), intent(in) :: config
  real(irealLUT), intent(in) :: tauz, w0
  integer(iintegers), intent(in) :: index_1d
  logical, intent(in) :: dir
  integer(iintegers), intent(in) :: src
  real(irealLUT), intent(in), dimension(:) :: T, S

  integer(mpiint), intent(out) :: ierr
  ierr=0

  if(any(T.lt.0._ireals)) ierr = ierr+3
  if(any(S.lt.0._ireals)) ierr = ierr+5
  if(any(T.gt.1._ireals)) ierr = ierr+7
  if(any(S.gt.1._ireals)) ierr = ierr+11

  select type(OPP)

  class is (t_optprop_LUT_wedge_5_8)
    call validate_plexrt()

  class is (t_optprop_LUT_rectilinear_wedge_5_8)
    call validate_plexrt()

  class is (t_optprop_LUT_wedge_18_8)
    call validate_plexrt()

  class default
    continue
  end select

  if(.not.dir .and. tauz*(1._irealLUT - w0) .lt. 1._irealLUT) then
    if(all(S.lt.tiny(S))) then
      call CHKERR(1_mpiint, 'Found all diff2diff coefficients to be zero but absorption is not that high.'// &
        'Something weird is going on. Please investigate!')
    endif
  endif

contains
  subroutine src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)
    real(irealLUT), intent(in) :: param_phi, param_theta
    logical, intent(out) :: lsrc(:)
    if(param_theta.le.0._irealLUT) then
      lsrc(1) = .True.
      lsrc(2:5) = .False.
      return
    endif
    if(param_phi.le.-1._irealLUT) then
      lsrc(1:3) = .True.
      lsrc(4:5) = .False.
      return
    endif
    if(param_phi.ge.1._irealLUT) then
      lsrc([1,2,4]) = .True.
      lsrc([3,5]) = .False.
      return
    endif
    if(param_phi.gt.-1._irealLUT .and. param_phi.lt.1._irealLUT) then
      lsrc(1:2) = .True.
      lsrc(3:5) = .False.
      return
    endif
    call CHKERR(1_mpiint, 'should not be here? '// &
      'param_phi '//ftoa(param_phi)// &
      'param_theta '//ftoa(param_theta))
  end subroutine
  subroutine validate_plexrt()
    real(irealLUT) :: param_phi, param_theta
    logical :: lsrc5(5)

    if(dir) then
      call get_sample_pnt_by_name_and_index(config, 'param_phi', index_1d, param_phi, ierr)
      call CHKERR(ierr, 'param_phi has to be present for wedge calculations')
      call get_sample_pnt_by_name_and_index(config, 'param_theta', index_1d, param_theta, ierr)
      call CHKERR(ierr, 'param_theta has to be present for wedge calculations')

      select type(OPP)
      class is (t_optprop_LUT_wedge_5_8)
        associate(lsrc => lsrc5)
          call src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)

          if(lsrc(src)) then ! I am a src face
            if(any(lsrc .and. T.gt.0._irealLUT)) ierr = ierr + 13 ! all srcs cannot have any incoming energy
          else
            if(any(T.gt.0._irealLUT)) ierr = ierr+17 ! if I am a dst, nobody should get any radiation anyway
          endif
        end associate
      class is (t_optprop_LUT_rectilinear_wedge_5_8)
        associate(lsrc => lsrc5)
          call src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)

          if(lsrc(src)) then ! I am a src face
            if(any(lsrc .and. T.gt.0._irealLUT)) ierr = ierr+13 ! all srcs cannot have any incoming energy
          else
            if(any(T.gt.0._irealLUT)) then
              ierr = ierr+17 ! if I am a dst, nobody should get any radiation anyway
            endif
          endif
        end associate
      end select
    endif
  end subroutine
end subroutine

subroutine LUT_bmc_wrapper(OPP, config, index_1d, src, dir, comm, S_diff, T_dir, S_tol, T_tol)
    class(t_optprop_LUT) :: OPP
    type(t_lut_config), intent(in) :: config
    integer(iintegers), intent(in) :: index_1d
    integer(iintegers), intent(in) :: src
    logical, intent(in) :: dir
    integer(mpiint), intent(in) :: comm

    real(irealLUT),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(irealLUT),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    real(ireal_dp), allocatable :: vertices(:)
    real(irealLUT) :: tauz, w0, g, phi, theta
    real(ireal_params) :: param_phi, param_theta

    logical :: lvalid
    integer(mpiint) :: ierr, ierr2

    call OPP%LUT_bmc_wrapper_determine_sample_pts(config, index_1d, dir, &
      vertices, tauz, w0, g, phi, theta, lvalid)

    if(.not.lvalid) then
      if(.not.dir) call CHKERR(1_mpiint, 'all diff2diff coeffs should be valid?')
      !if(ldebug) print *,'have an invalid coeff here, skipping computations', index_1d
      S_diff = 0
      T_dir = 0
      S_tol = 0
      T_tol = 0
      return
    endif

    call bmc_wrapper(OPP, src, vertices, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)

    call LUT_bmc_wrapper_validate(OPP, config, tauz, w0, index_1d, src, dir, T_dir, S_diff, ierr)
    if(ierr.ne.0) then
      ierr2 = ierr
      call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, ireal_params), &
        real(deg2rad(phi), ireal_params), real(deg2rad(theta), ireal_params), &
        param_phi, param_theta, ierr)
      print *,'LUTBMC :: calling bmc_get_coeff src',src,'tauz',tauz,w0,g,'angles',phi,theta,&
        ':ind1d',index_1d, ': T', T_dir, ':: S ', S_diff, ': verts', vertices, &
        'param phi/theta_afterwards', param_phi, param_theta
      call CHKERR(ierr2, 'found bad results for a given geometry')
    endif

    !if(ldebug) print *,'LUTBMC :: calling bmc_get_coeff src', src, &
    !  'tauz',tauz,w0,g,'angles',phi,theta, ':ind1d',index_1d, ':', &
    !  T_dir, ':(', T_tol, ') //', S_diff, ':(', S_tol, ')'
end subroutine

subroutine bmc_wrapper(OPP, src, vertices, tauz, w0, g, dir, phi, theta, comm, &
    S_diff, T_dir, S_tol, T_tol, inp_atol, inp_rtol)
    class(t_optprop_LUT) :: OPP
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: dir
    integer(mpiint),intent(in) :: comm
    real(irealLUT), intent(in) :: tauz, w0, g, phi, theta
    real(ireal_dp) :: vertices(:)

    real(irealLUT),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(irealLUT),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)
    real(irealLUT),intent(in),optional :: inp_atol, inp_rtol
    real(ireals) :: rS_diff(OPP%diff_streams),rT_dir(OPP%dir_streams)
    real(ireals) :: rS_tol (OPP%diff_streams),rT_tol(OPP%dir_streams)

    real(ireal_dp) :: bg(3), dz, atol, rtol

    atol = get_arg(stddev_atol, inp_atol)
    rtol = get_arg(stddev_rtol, inp_rtol)

    atol = atol-epsilon(atol)*10
    rtol = rtol-epsilon(rtol)*10

    dz = real(vertices(size(vertices)), kind(dz))

    bg(1) = tauz / dz * (1._irealLUT-w0)
    bg(2) = tauz / dz * w0
    bg(3) = g

    !print *,comm,'BMC :: calling bmc_get_coeff tauz',tauz,'w0,g',w0,g,phi,theta
    !print *,comm,'BMC :: calling bmc_get_coeff dz bg',vertices(size(vertices)),bg, '=>', sum(bg(1:2))*vertices(size(vertices)),'/',tauz
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(1:3) ,':', vertices(10:12)
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(4:6) ,':', vertices(13:15)
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(7:9) ,':', vertices(16:18)
    !print *,'area bot', triangle_area_by_vertices(vertices(1:3), vertices(4:6), vertices(7:9))
    !print *,'area top', triangle_area_by_vertices(vertices(10:12), vertices(13:15), vertices(16:18))
    !print *,'area_ratio:', triangle_area_by_vertices(vertices(1:3), vertices(4:6), vertices(7:9)) &
    !  / triangle_area_by_vertices(vertices(10:12), vertices(13:15), vertices(16:18))
    !call CHKERR(1_mpiint, 'DEBUG')

    call OPP%bmc%get_coeff(comm, bg, &
      src, dir, &
      real(phi, ireal_dp), real(theta, ireal_dp), &
      vertices,   &
      rS_diff, rT_dir, &
      rS_tol, rT_tol, &
      inp_atol=atol, inp_rtol=rtol)

    S_diff = real(rS_diff, irealLUT)
    T_dir  = real(rT_dir, irealLUT)
    S_tol  = real(rS_tol, irealLUT)
    T_tol  = real(rT_tol, irealLUT)
end subroutine

!function lin_index_to_param(idx,rng,N)
!    real(irealLUT) :: lin_index_to_param
!    real(irealLUT),intent(in) :: idx,rng(2)
!    integer(iintegers),intent(in) :: N
!    if(N.gt.i1) then
!      lin_index_to_param = rng(1) + (idx-1) * ( rng(2)-rng(1) ) / real(N-1, irealLUT)
!    else
!      lin_index_to_param = rng(1)
!    endif
!end function

subroutine populate_LUT_dim(dimname, N, lut_dim, vrange, preset)
  character(len=*),intent(in) :: dimname
  integer(iintegers), intent(in) :: N
  type(t_LUT_dim),intent(out) :: lut_dim
  real(irealLUT), optional :: vrange(:), preset(:)
  integer(iintegers) :: k
  if(allocated(lut_dim%v)) return ! already done
  allocate(lut_dim%v(N))

  if(present(vrange)) then
    do k=1,N
      lut_dim%v(k) = linspace(k, vrange, N)
    enddo
  elseif(present(preset)) then
    if(size(preset).ne.N) &
      call CHKERR(1_mpiint, 'Given preset size does not conform to proposed size N ' &
      //itoa(N)//' vs '//itoa(size(preset, kind=iintegers)))
    lut_dim%v = preset
  else
    call CHKERR(1_mpiint, 'Have to provide either a number of a preset for LUT dimension')
  endif

  lut_dim%vrange = [lut_dim%v(1), lut_dim%v(N)]
  lut_dim%dimname = trim(dimname)
  lut_dim%N = size(lut_dim%v)
end subroutine

subroutine set_parameter_space(OPP)
    class(t_optprop_LUT) :: OPP
    integer(mpiint) :: k, myid, ierr

    allocate(OPP%dirconfig)
    allocate(OPP%diffconfig)

    if(.not.lLUT_mockup) then
      select type(OPP)
        class is (t_optprop_LUT_1_2)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       i1, OPP%dirconfig%dims(5), vrange=real([0], irealLUT)) ! azimithally average in 1D
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_8_10)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_8_12)
            allocate(OPP%dirconfig%dims(6))
            !call populate_LUT_dim('tau',       i2, OPP%dirconfig%dims(1), vrange=real([1e-5,10.], irealLUT))
            !call populate_LUT_dim('w0',        i2, OPP%dirconfig%dims(2), vrange=real([.0,.99], irealLUT))
            !call populate_LUT_dim('aspect_zx', i2, OPP%dirconfig%dims(3), vrange=real([.1,2.], irealLUT))
            !call populate_LUT_dim('g',         i2, OPP%dirconfig%dims(4), vrange=real([0.,.5], irealLUT))
            !call populate_LUT_dim('phi',       i2, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            !call populate_LUT_dim('theta',     i2, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       i2, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     i3, OPP%dirconfig%dims(6), vrange=real([55,60,65], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            !call populate_LUT_dim('tau',       i2, OPP%diffconfig%dims(1), vrange=real([1e-5,10.], irealLUT))
            !call populate_LUT_dim('w0',        i2, OPP%diffconfig%dims(2), vrange=real([.0,.99], irealLUT))
            !call populate_LUT_dim('aspect_zx', i2, OPP%diffconfig%dims(3), vrange=real([.1,2.], irealLUT))
            !call populate_LUT_dim('g',         i2, OPP%diffconfig%dims(4), vrange=real([0.,.5], irealLUT))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_8_16)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_8_18)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_3_10)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_3_16)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_3_6)
            allocate(OPP%dirconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
            call populate_LUT_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
            allocate(OPP%diffconfig%dims(4))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)

        class is (t_optprop_LUT_wedge_5_8)
            allocate(OPP%dirconfig%dims(8))
            call populate_LUT_dim('tau',       size(preset_tau15,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau15)
            call populate_LUT_dim('w0',        size(preset_w010,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w010)
            call populate_LUT_dim('aspect_zx', size(preset_aspect18,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect18)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('wedge_coord_Cx', 7_iintegers, OPP%dirconfig%dims(5), &
                    vrange=real([.35,.65], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 7_iintegers, OPP%dirconfig%dims(6), &
                    vrange=real([0.7760254, 0.9560254], irealLUT))

            call populate_LUT_dim('param_phi', size(preset_param_phi19, kind=iintegers), &
                    OPP%dirconfig%dims(7), preset=preset_param_phi19)
            call populate_LUT_dim('param_theta', size(preset_param_theta13, kind=iintegers), &
                    OPP%dirconfig%dims(8), preset=preset_param_theta13)

            allocate(OPP%diffconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w021,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w021)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('wedge_coord_Cx', 7_iintegers, OPP%diffconfig%dims(5), &
                    vrange=real([.35,.65], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 7_iintegers, OPP%diffconfig%dims(6), &
                    vrange=real([0.7760254, 0.9560254], irealLUT))

        class is (t_optprop_LUT_rectilinear_wedge_5_8)
            allocate(OPP%dirconfig%dims(8))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w021,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w021)
            call populate_LUT_dim('aspect_zx', size(preset_aspect18,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect18)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('wedge_coord_Cx', 3_iintegers, OPP%dirconfig%dims(5), &
                    vrange=real([-0.000001,1.000001], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 2_iintegers, OPP%dirconfig%dims(6), &
                    vrange=real([.499999,1.000001], irealLUT))

            call populate_LUT_dim('param_phi', size(preset_param_phi19, kind=iintegers), &
                    OPP%dirconfig%dims(7), preset=preset_param_phi19)
            call populate_LUT_dim('param_theta', size(preset_param_theta13, kind=iintegers), &
                    OPP%dirconfig%dims(8), preset=preset_param_theta13)

            allocate(OPP%diffconfig%dims(6))
            call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
            call populate_LUT_dim('w0',        size(preset_w021,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w021)
            call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('wedge_coord_Cx', 3_iintegers, OPP%diffconfig%dims(5), &
                    vrange=real([-0.000001,1.000001], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 2_iintegers, OPP%diffconfig%dims(6), &
                    vrange=real([.499999,1.000001], irealLUT))

        class is (t_optprop_LUT_wedge_18_8)
            allocate(OPP%dirconfig%dims(7))
            call populate_LUT_dim('tau',       size(preset_tau15,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau15)
            call populate_LUT_dim('w0',        size(preset_w010,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w010)
            call populate_LUT_dim('aspect_zx', size(preset_aspect18,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect18)
            call populate_LUT_dim('wedge_coord_Cx', 7_iintegers, OPP%dirconfig%dims(4), &
                    vrange=real([.35,.65], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 7_iintegers, OPP%dirconfig%dims(5), &
                    vrange=real([0.7760254, 0.9560254], irealLUT))

            call populate_LUT_dim('param_phi', size(preset_param_phi19, kind=iintegers), &
                    OPP%dirconfig%dims(6), preset=preset_param_phi19)
            call populate_LUT_dim('param_theta', size(preset_param_theta13, kind=iintegers), &
                    OPP%dirconfig%dims(7), preset=preset_param_theta13)

            allocate(OPP%diffconfig%dims(5))
            call populate_LUT_dim('tau',       size(preset_tau20,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau20)
            call populate_LUT_dim('w0',        i2, OPP%diffconfig%dims(2), vrange=real([.0,.99999], irealLUT))
            call populate_LUT_dim('aspect_zx', size(preset_aspect7,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect7)
            call populate_LUT_dim('g',         size(preset_g6,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g6)
            call populate_LUT_dim('wedge_coord_Cx', 7_iintegers, OPP%diffconfig%dims(5), &
                    vrange=real([.35,.65], irealLUT))
            call populate_LUT_dim('wedge_coord_Cy', 7_iintegers, OPP%diffconfig%dims(6), &
                    vrange=real([0.7760254, 0.9560254], irealLUT))

        class default
          call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop_LUT object!')
      end select

    else ! do a mockup of a LUT
      select type(OPP)
        class is (t_optprop_LUT_1_2)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_8_10)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_8_12)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_8_16)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_8_18)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_3_10)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_3_16)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_3_6)
          call setup_pprts_mockup()
        class is (t_optprop_LUT_wedge_5_8)
          call setup_wedge_mockup()
        class is (t_optprop_LUT_rectilinear_wedge_5_8)
          call setup_wedge_mockup()
        class is (t_optprop_LUT_wedge_18_8)
          call setup_wedge_mockup()
        class default
          call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop_LUT object!')
      end select
    endif

    if(size(OPP%dirconfig%dims).gt.LUT_MAX_DIM) call CHKERR(1_mpiint, 'Parameter LUT_MAX_DIM is too small '//&
      'for the LUT you are using... please increase it to '// &
      itoa(size(OPP%dirconfig%dims))//' in src/optprop_parameters.f90')
    if(size(OPP%diffconfig%dims).gt.LUT_MAX_DIM) call CHKERR(1_mpiint, 'Parameter LUT_MAX_DIM is too small '//&
      'for the LUT you are using... please increase it to '// &
      itoa(size(OPP%diffconfig%dims))//' in src/optprop_parameters.f90')

    ! Determine offsets
    allocate(OPP%dirconfig%offsets(size(OPP%dirconfig%dims)))
    call ndarray_offsets(OPP%dirconfig%dims(:)%N, OPP%dirconfig%offsets)

    allocate(OPP%diffconfig%offsets(size(OPP%diffconfig%dims)))
    call ndarray_offsets(OPP%diffconfig%dims(:)%N, OPP%diffconfig%offsets)

    if(ldebug.and..False.) then
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      print *,myid,'set_parameter space dims:', size(OPP%diffconfig%dims), size(OPP%dirconfig%dims)
      do k=1,size(OPP%dirconfig%dims)
        print *,myid,'dim ',trim(OPP%dirconfig%dims(k)%dimname), OPP%dirconfig%offsets(k), &
                OPP%dirconfig%dims(k)%vrange, ':', OPP%dirconfig%dims(k)%v
      enddo
    endif
    contains
    subroutine setup_pprts_mockup()
      allocate(OPP%dirconfig%dims(6))
      call populate_LUT_dim('tau',       2_iintegers, OPP%dirconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_LUT_dim('w0',        2_iintegers, OPP%dirconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_LUT_dim('aspect_zx', 2_iintegers, OPP%dirconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_LUT_dim('g',         2_iintegers, OPP%dirconfig%dims(4), vrange=real([0,1], irealLUT))
      call populate_LUT_dim('phi',       2_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
      call populate_LUT_dim('theta',     2_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
      allocate(OPP%diffconfig%dims(4))
      call populate_LUT_dim('tau',       2_iintegers, OPP%diffconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_LUT_dim('w0',        2_iintegers, OPP%diffconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_LUT_dim('aspect_zx', 2_iintegers, OPP%diffconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_LUT_dim('g',         2_iintegers, OPP%diffconfig%dims(4), vrange=real([0,1], irealLUT))
    end subroutine
    subroutine setup_wedge_mockup()
      allocate(OPP%dirconfig%dims(7))
      call populate_LUT_dim('tau',       2_iintegers, OPP%dirconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_LUT_dim('w0',        2_iintegers, OPP%dirconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_LUT_dim('aspect_zx', 2_iintegers, OPP%dirconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_LUT_dim('wedge_coord_Cx', 2_iintegers, OPP%dirconfig%dims(4), &
        vrange=real([.35,.65], irealLUT))
      call populate_LUT_dim('wedge_coord_Cy', 2_iintegers, OPP%dirconfig%dims(5), &
        vrange=real([0.7760254, 0.9560254], irealLUT))
      call populate_LUT_dim('param_phi', size(preset_param_phi6, kind=iintegers), &
        OPP%dirconfig%dims(6), preset=preset_param_phi6)
      call populate_LUT_dim('param_theta', size(preset_param_theta4, kind=iintegers), &
        OPP%dirconfig%dims(7), preset=preset_param_theta4)
      allocate(OPP%diffconfig%dims(5))
      call populate_LUT_dim('tau',       2_iintegers, OPP%diffconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_LUT_dim('w0',        2_iintegers, OPP%diffconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_LUT_dim('aspect_zx', 2_iintegers, OPP%diffconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_LUT_dim('wedge_coord_Cx', 2_iintegers, OPP%diffconfig%dims(4), &
        vrange=real([.35,.65], irealLUT))
      call populate_LUT_dim('wedge_coord_Cy', 2_iintegers, OPP%diffconfig%dims(5), &
        vrange=real([0.7760254, 0.9560254], irealLUT))
    end subroutine
end subroutine

  subroutine scatter_LUTtables(OPP, comm)
      integer(mpiint) ,intent(in) :: comm
      class(t_optprop_LUT) :: OPP

      integer(mpiint) :: myid, ierr
      real(irealLUT), contiguous, pointer :: mmap_ptr(:,:)

      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)

      if (luse_memory_map) then
        mmap_ptr => NULL()
        if(associated(OPP%Sdiff%c)) then
          call arr_to_mmap(comm, trim(OPP%Sdiff%table_name_c(1))//'.Sdiff.mmap', mmap_ptr, ierr, OPP%Sdiff%c)
          deallocate(OPP%Sdiff%c)
        else
          call arr_to_mmap(comm, trim(OPP%Sdiff%table_name_c(1))//'.Sdiff.mmap', mmap_ptr, ierr)
        endif
        OPP%Sdiff%c => mmap_ptr

        call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
          '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
          '  You probably need to run the createLUT program first '//new_line('')// &
          '  or run with the additional option:'//new_line('')// &
          '    -skip_load_LUT no')

        mmap_ptr => NULL()
        if(associated(OPP%Sdir%c)) then
          call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr, OPP%Sdir%c)
          deallocate(OPP%Sdir%c)
        else
          call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr)
        endif
        OPP%Sdir%c => mmap_ptr

        call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
          '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
          '  You probably need to run the createLUT program first '//new_line('')// &
          '  or run with the additional option:'//new_line('')// &
          '    -skip_load_LUT no')

        mmap_ptr => NULL()
        if(associated(OPP%Tdir%c)) then
          call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr, OPP%Tdir%c)
          deallocate(OPP%Tdir%c)
        else
          call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr)
        endif
        OPP%Tdir%c => mmap_ptr

        call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
          '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
          '  You probably need to run the createLUT program first '//new_line('')// &
          '  or run with the additional option:'//new_line('')// &
          '    -skip_load_LUT no')

      else
        if(myid.eq.0) then
          if(.not.associated(OPP%Sdir%c)) then
            call CHKERR(1_mpiint, 'LUT data is not loaded. '//new_line('')// &
            '  Somehow we ended up in a situation where we would like to distribute'//new_line('')// &
            '  the LUT`s via MPI but rank 0 does not have the info. '//new_line('')// &
            '  Did you use the option: -skip_load_LUT ? '//new_line('')// &
            '  Maybe try to set it to -skip_load_LUT no')
          endif
        endif

        if( mpi_logical_or(comm, .not.associated(OPP%Sdir%c) )) &
          call imp_bcast(comm, OPP%Sdir%c, 0_mpiint)  ! DIRECT 2 DIRECT

        if( mpi_logical_or(comm, .not.associated(OPP%Tdir%c) )) &
          call imp_bcast(comm, OPP%Tdir%c, 0_mpiint)  ! DIRECT 2 DIFFUSE

        if( mpi_logical_or(comm, .not.associated(OPP%Sdiff%c) )) &
          call imp_bcast(comm, OPP%Sdiff%c, 0_mpiint)

      endif
  end subroutine

  subroutine print_configs(OPP)
    class(t_optprop_LUT) :: OPP
    integer(mpiint) :: myid, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)
    if(myid.eq.0) then
      print *,'Diffuse LUT config:'
      call print_LUT_config(OPP%diffconfig)
      print *,'Direct LUT config:'
      call print_LUT_config(OPP%dirconfig)
      print *,'----------------------'
    endif
  end subroutine

  subroutine print_LUT_config(config)
    type(t_LUT_config), allocatable, intent(in) :: config
    integer(iintegers) :: i
    integer(mpiint) :: myid, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)

    print *,myid,'LUT Config initialized', allocated(config)
    if(.not.allocated(config)) return

    print *,myid,'LUT Config Ndim', size(config%dims)
    do i = 1, size(config%dims)
      print *,myid,'Dimension '//itoa(i)//' '//trim(config%dims(i)%dimname)// &
              ' size '//itoa(config%dims(i)%N), '(', config%dims(i)%vrange, ')'
    enddo
  end subroutine

  subroutine check_if_samplepts_in_LUT_bounds(sample_pts, config)
    real(irealLUT),intent(in) :: sample_pts(:)
    type(t_LUT_config), allocatable, intent(in) :: config
    integer(mpiint) :: ierr, kdim

    ierr = 0
    do kdim = 1,size(sample_pts)
      if(sample_pts(kdim).lt.config%dims(kdim)%vrange(1).or.sample_pts(kdim).gt.config%dims(kdim)%vrange(2)) then
        print *,'ERROR value in dimension '//trim(config%dims(kdim)%dimname)// &
                ' ('//itoa(kdim)//') is outside of LUT range', &
                sample_pts(kdim), 'not in:', config%dims(kdim)%vrange
        ierr = ierr +1
      endif
    enddo

    if(ierr.ne.0) then
      call print_LUT_config(config)
      call CHKERR(ierr, 'Out of Bounds ERROR in LUT retrieval')
    endif
  end subroutine

  subroutine LUT_get_dir2dir(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(irealLUT),intent(in) :: sample_pts(:)
    real(irealLUT),intent(out):: C(:) ! dimension(OPP%dir_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%dirconfig%dims)) then
        call print_LUT_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%dirconfig%dims)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%dirconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Tdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Tdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%dir_streams
        norm = sum(C( src:size(C):OPP%dir_streams))
        if(real(norm).gt.one+1e-5_irealLUT) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        print *,'Error in dir2dir coeffs :: ierr',iierr,size(C),OPP%dir_streams,'::',C
        do src=1,OPP%dir_streams
          print *,'SUM dir2dir coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),&
                  ' :: coeff',C( src:size(C):OPP%dir_streams )
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
    !call CHKERR(1_mpiint, 'DEBUG')
  end subroutine

  subroutine LUT_get_dir2diff(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(irealLUT),intent(in) :: sample_pts(:)
    real(irealLUT),intent(out):: C(:) ! dimension(OPP%dir_streams*OPP%diff_streams)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%dirconfig%dims)) then
        call print_LUT_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%dirconfig%dims)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%dirconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Sdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Sdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%dir_streams ) )
        if(real(norm).gt.one+1e-5_irealLUT) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        do src=1,OPP%diff_streams
          print *,'SUM dir2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)), &
                  ' :: coeff',C( src:size(C):OPP%dir_streams )
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
  end subroutine

  subroutine LUT_get_diff2diff(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(irealLUT),intent(in) :: sample_pts(:)
    real(irealLUT),intent(out):: C(:) ! dimension(OPP%diff_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%diffconfig%dims)) then
        call print_LUT_config(OPP%diffconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%diffconfig%dims, kind=iintegers)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%diffconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%diffconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%diffconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Sdiff%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Sdiff%c, OPP%diffconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%diff_streams ) )
        if(norm.gt.one+1e-5_irealLUT) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        do src=1,OPP%diff_streams
          print *,'SUM diff2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%diff_streams)), &
                  ' :: coeff',C(src:size(C):OPP%diff_streams)
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
  end subroutine

end module
