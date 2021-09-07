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
  use iso_fortran_env, only: int64

#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: approx, &
                                rel_approx, imp_bcast, &
                                get_arg, toStr, &
                                cstr, char_arr_to_str, &
                                mpi_logical_and, mpi_logical_or, &
                                CHKERR, &
                                triangle_area_by_vertices, &
                                rad2deg, deg2rad, &
                                ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets, &
                                linspace

  use m_search, only: find_real_location

  use m_data_parameters, only: iintegers, mpiint, &
                               ireals, ireal_dp, irealLUT, ireal_params, &
                               one, zero, i0, i1, i2, i3, i10, nil, inil, &
                               imp_iinteger, imp_ireals, imp_irealLUT, imp_logical, &
                               default_str_len

  use m_optprop_parameters, only: &
    luse_memory_map, &
    ldebug_optprop, lut_basename, &
    LUT_dump_interval, &
    LUT_max_create_jobtime, &
    LUT_MAX_DIM, &
    delta_scale_truncate, &
    stddev_atol, stddev_rtol, &
    wedge_sphere_radius

  use m_boxmc, only: t_boxmc, &
                     t_boxmc_1_2, &
                     t_boxmc_3_6, t_boxmc_3_10, t_boxmc_3_16, &
                     t_boxmc_8_10, t_boxmc_8_12, t_boxmc_8_16, &
                     t_boxmc_8_18, &
                     t_boxmc_wedge_5_8, t_boxmc_wedge_18_8
  use m_tenstream_interpolation, only: interp_4d, interp_vec_simplex_nd
  use m_netcdfio

  use m_mmap, only: arr_to_mmap, munmap_mmap_ptr

  use m_boxmc_geometry, only: setup_default_wedge_geometry, setup_default_unit_cube_geometry

  use m_LUT_param_phi, only: param_phi_from_azimuth, azimuth_from_param_phi, &
                             iterative_phi_theta_from_param_phi_and_param_theta, &
                             param_phi_param_theta_from_phi_and_theta_withcoords, &
                             LUT_wedge_dz

  use m_optprop_base, only: &
    & t_optprop_base,                   &
    & t_op_config,                      &
    & find_op_dim_by_name,              &
    & set_op_param_space,               &
    & get_sample_pnt_by_name_and_index, &
    & check_if_samplepts_in_bounds,     &
    & print_op_config

  implicit none

  private
  public :: t_optprop_LUT, &
            t_optprop_LUT_1_2, &
            t_optprop_LUT_3_6, &
            t_optprop_LUT_3_10, &
            t_optprop_LUT_3_10_for_ANN, &
            t_optprop_LUT_3_16, &
            t_optprop_LUT_8_10, &
            t_optprop_LUT_8_12, &
            t_optprop_LUT_8_16, &
            t_optprop_LUT_8_18, &
            t_optprop_LUT_wedge_5_8, &
            t_optprop_LUT_rectilinear_wedge_5_8, &
            t_optprop_LUT_wedge_18_8, &
            azimuth_from_param_phi, param_phi_from_azimuth
  ! This module loads and generates the LUT-tables for Tenstream Radiation
  ! computations.
  ! It also holds functions for interpolation on the regular LUT grid.

  integer(mpiint) :: iierr
  integer(mpiint) :: mpierr

  type t_table
    real(irealLUT), pointer :: c(:, :) => null() ! depending on config has Ndim_1*Ndim_2*etc. many entries
    real(irealLUT), allocatable :: stddev_tol(:) ! maxval of tolerance for a given entry
    character(default_str_len), allocatable :: table_name_c(:)
    character(default_str_len), allocatable :: table_name_tol(:)
  end type

  type, abstract, extends(t_optprop_base) :: t_optprop_LUT
    type(t_table), allocatable :: Sdiff, Sdir, Tdir
    logical :: initialized = .false., optprop_LUT_debug = ldebug_optprop
    character(default_str_len) :: lutbasename

  contains
    procedure :: init
    procedure :: destroy
    procedure :: get_dir2dir => LUT_get_dir2dir
    procedure :: get_dir2diff => LUT_get_dir2diff
    procedure :: get_diff2diff => LUT_get_diff2diff
    procedure :: LUT_bmc_wrapper_determine_sample_pts
    procedure :: LUT_bmc_wrapper
    procedure :: scatter_LUTtables
    procedure :: createLUT
    procedure :: loadLUT_dir
    procedure :: loadLUT_diff
    procedure :: set_parameter_space
    procedure :: print_configs
  end type

  type, extends(t_optprop_LUT) :: t_optprop_LUT_1_2
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_3_6
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_3_10
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_3_10_for_ANN
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_3_16
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_8_10
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_8_12
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_8_16
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_8_18
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_wedge_5_8
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_rectilinear_wedge_5_8
  end type
  type, extends(t_optprop_LUT) :: t_optprop_LUT_wedge_18_8
  end type

  ! You should not need to change this... but feel free to play around...
  ! interp_mode 1 == nearest neighbour interpolation
  ! interp_mode 2 == linear interpolation
  integer(iintegers), parameter :: interp_mode = 2

  logical, parameter :: ldebug = .false.

contains

  subroutine init(OPP, comm, skip_load_LUT)
    class(t_optprop_LUT) :: OPP
    integer(mpiint), intent(in) :: comm
    logical, intent(in), optional :: skip_load_LUT

    integer(mpiint) :: comm_size, myid, ierr
    logical :: lskip_load_LUT, load_diffuse_LUT_first
    logical :: lskip_load_LUT_dir, lskip_load_LUT_diff
    logical :: lshow_LUT, lflg

    if (OPP%initialized) return

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if (OPP%optprop_LUT_debug .and. myid .eq. 0) print *, 'Initializing LUT`s...'

    if (.not. allocated(OPP%bmc)) then
      select type (OPP)
      class is (t_optprop_LUT_1_2)
        OPP%dir_streams = 1
        OPP%diff_streams = 2
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_1_2 :: OPP%bmc)

      class is (t_optprop_LUT_3_6)
        OPP%dir_streams = 3
        OPP%diff_streams = 6
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_3_6 :: OPP%bmc)

      class is (t_optprop_LUT_3_10_for_ANN)
        OPP%dir_streams = 3
        OPP%diff_streams = 10
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_3_10 :: OPP%bmc)

      class is (t_optprop_LUT_3_10)
        OPP%dir_streams = 3
        OPP%diff_streams = 10
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_3_10 :: OPP%bmc)

      class is (t_optprop_LUT_3_16)
        OPP%dir_streams = 3
        OPP%diff_streams = 16
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_3_16 :: OPP%bmc)

      class is (t_optprop_LUT_8_10)
        OPP%dir_streams = 8
        OPP%diff_streams = 10
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_8_10 :: OPP%bmc)

      class is (t_optprop_LUT_8_12)
        OPP%dir_streams = 8
        OPP%diff_streams = 12
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_8_12 :: OPP%bmc)

      class is (t_optprop_LUT_8_16)
        OPP%dir_streams = 8
        OPP%diff_streams = 16
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_8_16 :: OPP%bmc)

      class is (t_optprop_LUT_8_18)
        OPP%dir_streams = 8
        OPP%diff_streams = 18
        OPP%lutbasename = trim(lut_basename)
        allocate (t_boxmc_8_18 :: OPP%bmc)

      class is (t_optprop_LUT_wedge_5_8)
        OPP%dir_streams = 5
        OPP%diff_streams = 8
        OPP%lutbasename = trim(lut_basename)//'_wedge.Rsphere'//toStr(int(wedge_sphere_radius))
        allocate (t_boxmc_wedge_5_8 :: OPP%bmc)

      class is (t_optprop_LUT_rectilinear_wedge_5_8)
        OPP%dir_streams = 5
        OPP%diff_streams = 8
        OPP%lutbasename = trim(lut_basename)//'_rectilinear_wedge'
        allocate (t_boxmc_wedge_5_8 :: OPP%bmc)

      class is (t_optprop_LUT_wedge_18_8)
        OPP%dir_streams = 18
        OPP%diff_streams = 8
        OPP%lutbasename = trim(lut_basename)//'_wedge.Rsphere'//toStr(int(wedge_sphere_radius))
        allocate (t_boxmc_wedge_18_8 :: OPP%bmc)

      class default
        call CHKERR(1_mpiint, 'initialize LUT: unexpected type for optprop_LUT object!')
      end select

      call OPP%bmc%init(comm)
    end if

    call OPP%set_parameter_space()

    load_diffuse_LUT_first = .false.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-load_diffuse_LUT_first', load_diffuse_LUT_first, lflg, ierr); call CHKERR(ierr)

    lskip_load_LUT = get_arg(.true., skip_load_LUT)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-skip_load_LUT', lskip_load_LUT, lflg, ierr); call CHKERR(ierr)

    lskip_load_LUT_dir = lskip_load_LUT
    lskip_load_LUT_diff = lskip_load_LUT
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-skip_load_LUT_dir', lskip_load_LUT_dir, lflg, ierr); call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-skip_load_LUT_diff', lskip_load_LUT_diff, lflg, ierr); call CHKERR(ierr)

    if ((.not. lskip_load_LUT_dir .or. .not. lskip_load_LUT_diff) .and. myid .eq. 0) &
      & print *, 'loading and checking LUT`s from netCDF', &
      & .not. lskip_load_LUT_dir,.not. lskip_load_LUT_diff

    if (load_diffuse_LUT_first) then
      call OPP%loadLUT_diff(comm, lskip_load_LUT_diff)
      call OPP%loadLUT_dir(comm, lskip_load_LUT_dir)
    else
      call OPP%loadLUT_dir(comm, lskip_load_LUT_dir)
      call OPP%loadLUT_diff(comm, lskip_load_LUT_diff)
    end if

    call OPP%scatter_LUTtables(comm)

    OPP%initialized = .true.

    if (OPP%optprop_LUT_debug .and. myid .eq. 0) print *, 'Initializing LUT`s... finished'

    lshow_LUT = ldebug
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-LUT_view', lshow_LUT, lflg, ierr); call CHKERR(ierr)
    if (lshow_LUT .and. myid .eq. 0) call print_configs(OPP)
  end subroutine

  subroutine destroy(OPP, ierr)
    class(t_optprop_LUT) :: OPP
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if (allocated(OPP%bmc)) then
      call OPP%bmc%destroy(ierr); call CHKERR(ierr)
      deallocate (OPP%bmc)
    end if
    if (allocated(OPP%Tdir)) then
      if (luse_memory_map) then
        call munmap_mmap_ptr(OPP%Tdir%c, ierr); call CHKERR(ierr)
      end if
      deallocate (OPP%Tdir)
    end if
    if (allocated(OPP%Sdir)) then
      if (luse_memory_map) then
        call munmap_mmap_ptr(OPP%Sdir%c, ierr); call CHKERR(ierr)
      end if
      deallocate (OPP%Sdir)
    end if
    if (allocated(OPP%Sdiff)) then
      if (luse_memory_map) then
        call munmap_mmap_ptr(OPP%Sdiff%c, ierr); call CHKERR(ierr)
      end if
      deallocate (OPP%Sdiff)
    end if
    if (allocated(OPP%dirconfig)) deallocate (OPP%dirconfig)
    if (allocated(OPP%diffconfig)) deallocate (OPP%diffconfig)
    OPP%initialized = .false.
    if (OPP%optprop_LUT_debug) print *, 'Destroyed LUTs'
  end subroutine

  function gen_lut_basename(prefix, config) result(lutname)
    character(len=default_str_len) :: lutname
    character(len=*), intent(in) :: prefix
    type(t_op_config), intent(in) :: config
    integer(iintegers) :: k
    lutname = trim(prefix)
    do k = 1, size(config%dims)
      lutname = trim(lutname)//'.'//trim(config%dims(k)%dimname)//toStr(config%dims(k)%N)
    end do
    lutname = trim(lutname)//'.ds'//toStr(int(delta_scale_truncate * 1000))
  end function

  subroutine load_table_from_netcdf(table, istat)
    type(t_table), intent(inout) :: table
    integer(mpiint), intent(out) :: istat

    integer(iintegers) :: k
    integer(mpiint) :: ierr
    istat = 0

    call ncload(table%table_name_tol, table%stddev_tol, ierr); istat = istat + ierr
    if (ierr .eq. 0) then
      if (ldebug) print *, 'we were able to load stddev from file '//new_line('')// &
        cstr(char_arr_to_str(table%table_name_tol, '/'), 'green')//new_line('')// &
        ' but still have to check if they all have a good enough precision / noise level...'
      do k = 1, size(table%stddev_tol, kind=iintegers)
        if (table%stddev_tol(k) .gt. stddev_atol + 10 * epsilon(stddev_atol) &
            .or. table%stddev_tol(k) .lt. 0._ireallut) then
          istat = istat + 1
          if (istat .ne. 0) then
            print *, 'coefficients stddev not good enough', k, table%stddev_tol(k), &
              'should be positive and less than', stddev_atol + 10 * epsilon(stddev_atol)
            exit
          end if
        end if
      end do
    else
      print *, 'loading stddev_tolerances from '//new_line('')// &
        cstr(char_arr_to_str(table%table_name_tol, '/'), 'purple')//new_line('')// &
        ' failed', istat, ierr
    end if
    if (istat .ne. 0) &
      print *, 'Test if coeffs in '//new_line('')// &
      cstr(char_arr_to_str(table%table_name_tol, '/'), 'purple')//new_line('')// &
      ' are good results in:', istat

    call ncload(table%table_name_c, table%c, ierr); istat = istat + ierr
    if (ierr .eq. 0) then ! we were able to load coeffs but still have to check if they are all ok...
      do k = 1, size(table%c, dim=2, kind=iintegers)
        if (any(table%c(:, k) .gt. one) .or. any(table%c(:, k) .lt. zero)) then
          istat = istat + 2
          if (istat .ne. 0) then
            print *, 'coefficients not good enough', k, table%c(:, k), 'should be between [0,1]'
            exit
          end if
        end if
      end do
    else
      print *, 'loading coeffs from '//new_line('')// &
        cstr(char_arr_to_str(table%table_name_c, '/'), 'purple')//new_line('')// &
        ' failed', istat, ierr
    end if

    if (istat .ne. 0) then
      print *, 'Test if coeffs in '//new_line('')// &
        cstr(char_arr_to_str(table%table_name_c, '/'), 'purple')//new_line('')// &
        ' are good results in:', istat
    end if
  end subroutine

  subroutine loadLUT_diff(OPP, comm, skip_load_LUT)
    class(t_optprop_LUT) :: OPP
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: skip_load_LUT

    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if (allocated(OPP%Sdiff)) return ! already loaded
    allocate (OPP%Sdiff)
    OPP%Sdiff%c => null()
    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('_diffuse_'//toStr(OPP%diff_streams), OPP%diffconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'diffuse'

    str(3) = 'Stol'; allocate (OPP%Sdiff%table_name_tol(size(str)), source=str)
    str(3) = 'S'; allocate (OPP%Sdiff%table_name_c(size(str)), source=str)

    if (skip_load_LUT) return

    if (myid .eq. 0) then
      call load_table_from_netcdf(OPP%Sdiff, iierr); errcnt = errcnt + iierr
      call write_pspace(str(1), OPP%diffconfig)
    end if

    call mpi_bcast(errcnt, 1_mpiint, imp_iinteger, 0_mpiint, comm, mpierr); call CHKERR(mpierr)

    if (errcnt .ne. 0) then ! something went wrong loading the LUT
      call OPP%createLUT(comm, OPP%diffconfig, OPP%Sdiff)
    end if

    if (allocated(OPP%Sdiff%stddev_tol)) deallocate (OPP%Sdiff%stddev_tol)
    if (ldebug) print *, 'Loaded Diff2Diff LUT with shape:', shape(OPP%Sdiff%c)
  end subroutine

  subroutine loadLUT_dir(OPP, comm, skip_load_LUT)
    class(t_optprop_LUT) :: OPP
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: skip_load_LUT

    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if (allocated(OPP%Sdir) .and. allocated(OPP%Tdir)) return ! already loaded
    allocate (OPP%Sdir)
    allocate (OPP%Tdir)
    OPP%Sdir%c => null()
    OPP%Tdir%c => null()

    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('_direct_'//toStr(OPP%dir_streams)//'_'//toStr(OPP%diff_streams), OPP%dirconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'direct'

    str(3) = 'Stol'; allocate (OPP%Sdir%table_name_tol(size(str)), source=str)
    str(3) = 'S'; allocate (OPP%Sdir%table_name_c(size(str)), source=str)

    str(3) = 'Ttol'; allocate (OPP%Tdir%table_name_tol(size(str)), source=str)
    str(3) = 'T'; allocate (OPP%Tdir%table_name_c(size(str)), source=str)

    if (skip_load_LUT) return

    if (myid .eq. 0) then
      call load_table_from_netcdf(OPP%Sdir, iierr); errcnt = errcnt + iierr
      call load_table_from_netcdf(OPP%Tdir, iierr); errcnt = errcnt + iierr
      call write_pspace(str(1), OPP%dirconfig)
    end if

    call mpi_bcast(errcnt, 1_mpiint, imp_iinteger, 0_mpiint, comm, mpierr); call CHKERR(mpierr)

    if (errcnt .ne. 0) then ! something went wrong loading the LUT
      call OPP%createLUT(comm, OPP%dirconfig, OPP%Sdir, OPP%Tdir)
    end if

    if (allocated(OPP%Tdir%stddev_tol)) deallocate (OPP%Tdir%stddev_tol)
    if (allocated(OPP%Sdir%stddev_tol)) deallocate (OPP%Sdir%stddev_tol)
    if (ldebug) print *, 'Loaded Dir2Dir LUT with shape:', shape(OPP%Tdir%c)
    if (ldebug) print *, 'Loaded Dir2Diff LUT with shape:', shape(OPP%Sdir%c)
  end subroutine

  subroutine write_pspace(fname, config)
    character(len=*), intent(in) :: fname
    type(t_op_config), intent(in) :: config
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
      if (allocated(existing_values)) deallocate (existing_values)
      call ncload(groups, existing_values, ierr)
      if (ierr .eq. 0) then
        if (.not. all(approx(existing_values, config%dims(kdim)%v, sqrt(epsilon(1._ireallut)) * 10))) then
          print *, kdim, trim(groups(1)), trim(groups(3)), new_line(''), &
            ':existing', existing_values, new_line(''), &
            ':new', config%dims(kdim)%v
          call CHKERR(1_mpiint, 'Dimensions of LUT and in optprop_parameters definition do not match!')
        end if
      else ! Otherwise, just save the current ones
        if (ldebug) print *, kdim, 'Writing pspace for ', trim(groups(1)), ':', trim(groups(2)), &
          ':', trim(groups(3)), trim(groups(4)), '::', config%dims(kdim)%v
        call ncwrite(groups, config%dims(kdim)%v, ierr); call CHKERR(ierr)
      end if

      groups(4) = 'range'
      if (allocated(existing_values)) deallocate (existing_values)
      call ncload(groups, existing_values, ierr)
      if (ierr .eq. 0) then
        if (.not. all(approx(existing_values, config%dims(kdim)%vrange, sqrt(epsilon(1._ireallut)) * 10))) then
          call CHKERR(1_mpiint, 'Range of dimensions of LUT and in optprop_parameters definition do not match!')
        end if
      else ! Otherwise, just save the current ones
        if (ldebug) print *, kdim, 'Writing pspace for ', trim(groups(1)), ':', trim(groups(2)), &
          ':', trim(groups(3)), trim(groups(4)), '::', config%dims(kdim)%vrange
        call ncwrite(groups, config%dims(kdim)%vrange, ierr); call CHKERR(ierr)
      end if
    end do
  end subroutine

  subroutine createLUT(OPP, comm, config, S, T)
    class(t_optprop_LUT), intent(in) :: OPP
    integer(mpiint), intent(in) :: comm
    type(t_op_config), intent(in) :: config
    type(t_table) :: S
    type(t_table), optional :: T

    logical :: gotmsg
    integer(mpiint) :: status(MPI_STATUS_SIZE)

    integer(mpiint), parameter :: READYMSG = 1, HAVERESULTSMSG = 2, WORKMSG = 3, FINALIZEMSG = 4, RESULTMSG = 5

    integer(iintegers) :: idummy, lutindex
    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if (myid .eq. 0) then
      call prepare_table_space(OPP, config, S, T)
    end if

    if (myid .le. 0 .and. comm_size .le. 1) &
      call CHKERR(1_mpiint, 'At the moment creation of direct Lookuptable needs at least two mpi-ranks to work...'// &
                  'please run with more ranks.')

    if (myid .eq. 0) then
      call master(S, T)
      print *, 'done calculating coefficients. present(T)?', present(T)
    else
      call worker(config)
    end if

  contains
    subroutine master(S, T)
      type(t_table), intent(inout) :: S
      type(t_table), intent(inout), optional :: T

      integer(iintegers) :: total_size, cnt, finalizedworkers

      integer(iintegers) :: Nsrc
      integer(iintegers) :: isrc, idst, ind

      logical :: ldoneS, ldoneT

      real(irealLUT), allocatable, dimension(:, :) :: S_diff, T_dir, S_tol, T_tol

      real :: starttime, lastsavetime, now
      integer :: clock_count, clock_count_rate
      call system_clock(clock_count, clock_count_rate)
      starttime = clock_count / clock_count_rate
      lastsavetime = starttime

      finalizedworkers = 0
      if (present(T)) then
        Nsrc = OPP%dir_streams
      else
        Nsrc = OPP%diff_streams
      end if
      total_size = size(S%c, dim=2) ! we have Ndims work packages with Nsrc bmc calls each

      allocate (S_diff(OPP%diff_streams, Nsrc), S_tol(OPP%diff_streams, Nsrc), &
                T_dir(OPP%dir_streams, Nsrc), T_tol(OPP%dir_streams, Nsrc))

      cnt = 1
      do ! endless loop until everything has been done

        ! Check if we already calculated the coefficients
        if (cnt .le. total_size) then

          ldoneS = all([ &
                       all(S%c(:, cnt) .ge. zero), &
                       all(S%c(:, cnt) .le. one), &
                       S%stddev_tol(cnt) .le. stddev_atol, &
                       S%stddev_tol(cnt) .ge. 0._ireallut])
          !print *,'ldoneS', cnt, ldoneS, ':', S%stddev_tol(cnt)

          if (present(T)) then
            ldoneT = all([ &
                         all(T%c(:, cnt) .ge. zero), &
                         all(T%c(:, cnt) .le. one), &
                         T%stddev_tol(cnt) .le. stddev_atol, &
                         T%stddev_tol(cnt) .ge. 0._ireallut])
          else
            ldoneT = .true.
          end if

          if (ldoneS .and. ldoneT) then
            if (mod(cnt - 1, max(i1, total_size / 100)) .eq. 0) & !every 1 percent report status
              print *, 'Resuming from LUT... ', cnt / max(i1, total_size / 100), '%'
            cnt = cnt + 1
            cycle
          end if
        end if

        ! Now that we know we got something to do, lets find a suitable worker
        gotmsg = .false.
        call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

        if (gotmsg) then

          select case (status(MPI_TAG))

          case (READYMSG)

            ! capture the READY MSG -- we should not leave messages hanging around.
            call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), READYMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)

            if (cnt .le. total_size) then ! we got something to do for a worker -- send him...
              lutindex = cnt
              call mpi_send(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), WORKMSG, &
                            comm, mpierr); call CHKERR(mpierr)
              call mpi_send(present(T), 1_mpiint, imp_logical, status(MPI_SOURCE), WORKMSG, &
                            comm, mpierr); call CHKERR(mpierr)

            else ! no more work to do... tell the worker to quit
              call mpi_send(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                            comm, mpierr); call CHKERR(mpierr)
            end if
            cnt = cnt + 1

          case (HAVERESULTSMSG)
            call mpi_recv(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)

            call mpi_recv(S_diff, size(S_diff, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            call mpi_recv(S_tol, size(S_tol, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            call mpi_recv(T_dir, size(T_dir, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            call mpi_recv(T_tol, size(T_tol, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)

            ! Sort coefficients into destination ordering and put em in LUT
            do isrc = 1, Nsrc
              do idst = 1, OPP%diff_streams
                ind = (idst - 1) * Nsrc + isrc
                S%c(ind, lutindex) = S_diff(idst, isrc)
              end do

              if (present(T)) then
                do idst = 1, OPP%dir_streams
                  ind = (idst - 1) * Nsrc + isrc
                  T%c(ind, lutindex) = T_dir(idst, isrc)
                end do
              end if
            end do ! isrc

            S%stddev_tol(lutindex) = maxval(S_tol)
            if (present(T)) then
              T%stddev_tol(lutindex) = maxval(T_tol)
            end if

            if (.false. .and. ldebug) call random_print_coeffs(lutindex, S_diff, T_dir, S_tol, T_tol)

            if (mod(lutindex - 1, max(i1, total_size / 1000_iintegers)) .eq. 0) & !every .1 percent report status
              print *, 'Calculated LUT...', lutindex, &
              real(lutindex - 1, irealLUT) * 100._ireallut / real(total_size, irealLUT), '%'

            call system_clock(clock_count, clock_count_rate)
            now = clock_count / clock_count_rate
            if ((now - lastsavetime) .gt. LUT_dump_interval .or. (now - starttime) .gt. LUT_max_create_jobtime) then !every 30 minutes wall clock time, dump the LUT.
              print *, 'Dumping LUT after ', (now - lastsavetime) / 60, 'minutes'
              if (present(T)) then
                print *, 'Writing table to file... ', char_arr_to_str(T%table_name_c)
                call ncwrite(T%table_name_c, T%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
                print *, 'Writing table to file... ', char_arr_to_str(T%table_name_tol)
                call ncwrite(T%table_name_tol, T%stddev_tol, iierr); call CHKERR(iierr, 'Could not write Table to file')
              end if
              print *, 'Writing table to file... ', char_arr_to_str(S%table_name_c)
              call ncwrite(S%table_name_c, S%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
              print *, 'Writing table to file... ', char_arr_to_str(S%table_name_tol)
              call ncwrite(S%table_name_tol, S%stddev_tol, iierr); call CHKERR(iierr, 'Could not write Table to file')
              print *, 'done writing!', iierr
              lastsavetime = now ! reset the countdown
              if ((now - starttime) .gt. LUT_max_create_jobtime) then
                call CHKERR(int(now - starttime, mpiint), 'Maximum duration of create_LUT reached ... exiting')
              end if
            end if

          case (FINALIZEMSG)
            call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            finalizedworkers = finalizedworkers + 1
            if (finalizedworkers .eq. comm_size - 1) exit ! all work is done

          end select
        end if
      end do

      print *, 'Writing table to file... '
      call ncwrite(S%table_name_c, S%c, iierr)
      call ncwrite(S%table_name_tol, S%stddev_tol, iierr)
      if (present(T)) then
        call ncwrite(T%table_name_c, T%c, iierr)
        call ncwrite(T%table_name_tol, T%stddev_tol, iierr)
        print *, 'done writing!', iierr, ':: max_atol S', maxval(S%stddev_tol), 'max_atol T', maxval(T%stddev_tol)
      else
        print *, 'done writing!', iierr, ':: max_atol S', maxval(S%stddev_tol)
      end if
    end subroutine
    subroutine worker(config)
      type(t_op_config), intent(in) :: config
      integer(iintegers) :: isrc, Nsrc
      logical :: ldir

      real(irealLUT), allocatable, dimension(:, :) :: S_diff, T_dir, S_tol, T_tol

      ! workers send READY message to master
      call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, READYMSG, comm, mpierr); call CHKERR(mpierr)

      do
        ! ask what to do
        gotmsg = .false.
        call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

        if (gotmsg) then

          select case (status(MPI_TAG))

          case (WORKMSG)
            ! wait for work to arrive
            call mpi_recv(lutindex, 1_mpiint, imp_iinteger, 0_mpiint, WORKMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            call mpi_recv(ldir, 1_mpiint, imp_logical, 0_mpiint, WORKMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)

            if (ldir) then
              Nsrc = OPP%dir_streams
            else
              Nsrc = OPP%diff_streams
            end if
            allocate (S_diff(OPP%diff_streams, Nsrc), S_tol(OPP%diff_streams, Nsrc), &
                      T_dir(OPP%dir_streams, Nsrc), T_tol(OPP%dir_streams, Nsrc))

            do isrc = 1, Nsrc
              call OPP%LUT_bmc_wrapper(config, lutindex, isrc, ldir, &
                                       mpi_comm_self, S_diff(:, isrc), T_dir(:, isrc), S_tol(:, isrc), T_tol(:, isrc))

              if (maxval(S_tol(:, isrc)) .gt. stddev_atol + sqrt(epsilon(stddev_atol)) &
                  .or. maxval(T_tol(:, isrc)) .gt. stddev_atol + sqrt(epsilon(stddev_atol))) then
                print *, 'Ttol', T_tol(:, isrc)
                print *, 'Stol', S_tol(:, isrc)
                call CHKERR(1_mpiint, 'BOXMC violated stddev constraints!')
              end if
            end do

            !print *,'Computed isrc',isrc,'aspect',aspect_zx,'tau',tau_z, w0, g,':', phi, theta
            !print *,myid,'Computed values for ',lutindex, isrc, ldir

            call mpi_send(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, &
                          comm, mpierr); call CHKERR(mpierr)
            call mpi_send(S_diff, size(S_diff, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, mpierr); call CHKERR(mpierr)
            call mpi_send(S_tol, size(S_tol, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, mpierr); call CHKERR(mpierr)
            call mpi_send(T_dir, size(T_dir, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, mpierr); call CHKERR(mpierr)
            call mpi_send(T_tol, size(T_tol, kind=mpiint), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, &
                          comm, mpierr); call CHKERR(mpierr)

            call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, READYMSG, &
                          comm, mpierr); call CHKERR(mpierr)

            deallocate (S_diff, S_tol, T_dir, T_tol)

          case (FINALIZEMSG)
            call mpi_recv(idummy, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, &
                          comm, status, mpierr); call CHKERR(mpierr)
            call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, &
                          comm, mpierr); call CHKERR(mpierr)
            exit

          end select

        end if !gotmsg
      end do
    end subroutine
    subroutine random_print_coeffs(lutindex, S_diff, T_dir, S_tol, T_tol)
      integer(iintegers) :: lutindex
      real(irealLUT), dimension(:, :) :: S_diff, T_dir, S_tol, T_tol
      real(ireals) :: R
      real(ireals), parameter :: chance = 1e-4
      integer(iintegers) :: isrc
      call random_number(R)
      if (R .lt. chance) then
        do isrc = 1, size(T_dir, dim=1)
          print *, 'lutindex '//toStr(lutindex)//' src '//toStr(isrc)//' :T', T_dir(:, isrc)
        end do
        print *, 'Min/Max T_tol ', minval(T_tol), maxval(T_tol)
        do isrc = 1, size(T_dir, dim=1)
          print *, 'lutindex '//toStr(lutindex)//' src '//toStr(isrc)//' :S', S_diff(:, isrc)
        end do
        print *, 'Min/Max S_tol ', minval(S_tol), maxval(S_tol)
      end if
    end subroutine
  end subroutine createLUT

  subroutine prepare_table_space(OPP, config, S, T)
    class(t_optprop_LUT) :: OPP
    type(t_op_config), intent(in) :: config
    type(t_table), intent(inout) :: S
    type(t_table), intent(inout), optional :: T

    integer(int64) :: entries, bytesize
    integer(iintegers) :: errcnt

    if (present(T)) then
      entries = &
        int(OPP%dir_streams**2, int64) * int(product(config%dims(:)%N), int64) + &
        int(OPP%diff_streams * OPP%dir_streams, int64) * int(product(config%dims(:)%N), int64) + &
        2_int64 * int(product(config%dims(:)%N), int64) ! tolerances
    else
      entries = int(OPP%diff_streams**2 + 1, int64) * int(product(config%dims(:)%N), int64)
    end if
    bytesize = int(c_sizeof(1._ireallut), int64) * entries
    print *, 'Allocating Space for LUTs '//toStr(entries)// &
      ' entries ( '//toStr(real(bytesize, irealLUT) / 1024._ireallut**3)//' Gb) ...'

    errcnt = 0
    if (present(T)) then
      if (.not. associated(S%c)) &
        allocate (S%c(OPP%diff_streams * OPP%dir_streams, product(config%dims(:)%N)), source=-1._ireallut)
      if (.not. associated(T%c)) &
        allocate (T%c(OPP%dir_streams * OPP%dir_streams, product(config%dims(:)%N)), source=-1._ireallut)

      if (.not. allocated(S%stddev_tol)) allocate (S%stddev_tol(product(config%dims(:)%N)), source=huge(-1._ireallut))
      if (.not. allocated(T%stddev_tol)) allocate (T%stddev_tol(product(config%dims(:)%N)), source=huge(-1._ireallut))
    else
      if (.not. associated(S%c)) allocate (S%c(OPP%diff_streams**2, product(config%dims(:)%N)), source=-1._ireallut)
      if (.not. allocated(S%stddev_tol)) allocate (S%stddev_tol(product(config%dims(:)%N)), source=huge(-1._ireallut))
    end if
    print *, 'Allocating Space for LUTs '//toStr(entries)// &
      ' entries ( '//toStr(real(bytesize, irealLUT) / 1024._ireallut**3)//' Gb) ... done'
  end subroutine

  subroutine LUT_bmc_wrapper_determine_sample_pts(OPP, config, index_1d, dir, &
                                                  vertices, tauz, w0, g, phi, theta, lvalid_entry)
    class(t_optprop_LUT) :: OPP
    type(t_op_config), intent(in) :: config
    integer(iintegers), intent(in) :: index_1d
    logical, intent(in) :: dir
    real(ireal_dp), intent(out), allocatable :: vertices(:)
    real(irealLUT), intent(out) :: tauz, w0, g, phi, theta
    logical, intent(out) :: lvalid_entry

    real(irealLUT) :: aspect_zx, aspect_zy
    integer(mpiint) :: ierr

    lvalid_entry = .true.

    call get_sample_pnt_by_name_and_index(config, 'tau', index_1d, tauz, ierr)
    call CHKERR(ierr, 'tauz has to be present')

    call get_sample_pnt_by_name_and_index(config, 'w0', index_1d, w0, ierr)
    call CHKERR(ierr, 'w0 has to be present')

    call get_sample_pnt_by_name_and_index(config, 'g', index_1d, g, ierr)
    if (ierr .ne. 0) then
      g = zero ! defaults to isotropic scattering
    end if

    call get_sample_pnt_by_name_and_index(config, 'aspect_zx', index_1d, aspect_zx, ierr)
    call CHKERR(ierr, 'aspect_zx has to be present')

    call get_sample_pnt_by_name_and_index(config, 'aspect_zy', index_1d, aspect_zy, ierr)
    if (ierr .ne. 0) then
      aspect_zy = aspect_zx ! set dy = dx
    end if

    ! default values for angles in case of diff2diff and such
    phi = 0
    theta = -1

    ! First define Default vertices for Cube Geometry
    select type (OPP)
    class is (t_optprop_LUT_1_2)
      call prep_pprts()

    class is (t_optprop_LUT_3_6)
      call prep_pprts()

    class is (t_optprop_LUT_3_10_for_ANN)
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
                                            real(aspect_zx / aspect_zy, ireal_dp), real(aspect_zx, ireal_dp), &
                                            vertices)
      if (dir) then
        call get_sample_pnt_by_name_and_index(config, 'phi', index_1d, phi, ierr)
        call CHKERR(ierr, 'phi has to be present for direct calculations')
        call get_sample_pnt_by_name_and_index(config, 'theta', index_1d, theta, ierr)
        call CHKERR(ierr, 'theta has to be present for direct calculations')
      end if
    end subroutine
    subroutine prep_plexrt(sphere_radius)
      real(ireal_dp), intent(in), optional :: sphere_radius
      real(irealLUT) :: param_phi, param_theta
      real(irealLUT), dimension(2) :: wedge_C
      real(ireal_params) :: rphi, rtheta
      real(irealLUT) :: Atop, dz

      wedge_C(1) = .5_ireallut
      call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cx', index_1d, wedge_C(1), ierr)
      call CHKERR(ierr, 'wedge_coord_Cx has to be present for wedge calculations')

      wedge_C(2) = 0.8660254_ireallut
      call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cy', index_1d, wedge_C(2), ierr)
      call CHKERR(ierr, 'wedge_coord_Cy has to be present for wedge calculations')

      Atop = triangle_area_by_vertices([0._ireallut, 0._ireallut], [1._ireallut, 0._ireallut], wedge_C)
      dz = LUT_wedge_dz(Atop, aspect_zx)

      call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], &
                                        [1._ireal_dp, 0._ireal_dp], &
                                        real(wedge_C, ireal_dp), &
                                        real(dz, ireal_dp), &
                                        vertices, sphere_radius=sphere_radius)

      if (dir) then
        call get_sample_pnt_by_name_and_index(config, 'param_phi', index_1d, param_phi, ierr)
        call CHKERR(ierr, 'param_phi has to be present for wedge calculations')
        call get_sample_pnt_by_name_and_index(config, 'param_theta', index_1d, param_theta, ierr)
        call CHKERR(ierr, 'param_theta has to be present for wedge calculations')
        call iterative_phi_theta_from_param_phi_and_param_theta( &
          real(vertices, ireal_params), &
          real(param_phi, ireal_params), &
          real(param_theta, ireal_params), &
          rphi, rtheta, ierr); call CHKERR(ierr)
        phi = real(rad2deg(rphi), irealLUT)
        theta = real(rad2deg(rtheta), irealLUT)

        if (param_theta .le. 0._ireals) then
          if (param_phi .lt. -1._ireals .or. param_phi .gt. 1._ireals) then
            lvalid_entry = .false. ! this is a dead zone in the param space. Should never happen in a mesh
          end if
        end if
      end if
    end subroutine
  end subroutine

  subroutine LUT_bmc_wrapper_validate(OPP, config, tauz, w0, index_1d, src, dir, T, S, ierr)
    class(t_optprop_LUT) :: OPP
    type(t_op_config), intent(in) :: config
    real(irealLUT), intent(in) :: tauz, w0
    integer(iintegers), intent(in) :: index_1d
    logical, intent(in) :: dir
    integer(iintegers), intent(in) :: src
    real(irealLUT), intent(in), dimension(:) :: T, S

    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if (any(T .lt. 0._ireals)) ierr = ierr + 3
    if (any(S .lt. 0._ireals)) ierr = ierr + 5
    if (any(T .gt. 1._ireals)) ierr = ierr + 7
    if (any(S .gt. 1._ireals)) ierr = ierr + 11

    select type (OPP)

    class is (t_optprop_LUT_wedge_5_8)
      call validate_plexrt()

    class is (t_optprop_LUT_rectilinear_wedge_5_8)
      call validate_plexrt()

    class is (t_optprop_LUT_wedge_18_8)
      call validate_plexrt()

    class default
      continue
    end select

    if (.not. dir .and. tauz * (1._ireallut - w0) .lt. 1._ireallut) then
      if (all(S .lt. tiny(S))) then
        call CHKERR(1_mpiint, 'Found all diff2diff coefficients to be zero but absorption is not that high.'// &
                    'Something weird is going on. Please investigate!')
      end if
    end if

  contains
    subroutine src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)
      real(irealLUT), intent(in) :: param_phi, param_theta
      logical, intent(out) :: lsrc(:)
      if (param_theta .le. 0._ireallut) then
        lsrc(1) = .true.
        lsrc(2:5) = .false.
        return
      end if
      if (param_phi .le. -1._ireallut) then
        lsrc(1:3) = .true.
        lsrc(4:5) = .false.
        return
      end if
      if (param_phi .ge. 1._ireallut) then
        lsrc([1, 2, 4]) = .true.
        lsrc([3, 5]) = .false.
        return
      end if
      if (param_phi .gt. -1._ireallut .and. param_phi .lt. 1._ireallut) then
        lsrc(1:2) = .true.
        lsrc(3:5) = .false.
        return
      end if
      call CHKERR(1_mpiint, 'should not be here? '// &
                  'param_phi '//toStr(param_phi)// &
                  'param_theta '//toStr(param_theta))
    end subroutine
    subroutine validate_plexrt()
      real(irealLUT) :: param_phi, param_theta
      logical :: lsrc5(5)

      if (dir) then
        call get_sample_pnt_by_name_and_index(config, 'param_phi', index_1d, param_phi, ierr)
        call CHKERR(ierr, 'param_phi has to be present for wedge calculations')
        call get_sample_pnt_by_name_and_index(config, 'param_theta', index_1d, param_theta, ierr)
        call CHKERR(ierr, 'param_theta has to be present for wedge calculations')

        select type (OPP)
        class is (t_optprop_LUT_wedge_5_8)
          associate (lsrc => lsrc5)
            call src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)

            if (lsrc(src)) then ! I am a src face
              if (any(lsrc .and. T .gt. 0._ireallut)) ierr = ierr + 13 ! all srcs cannot have any incoming energy
            else
              if (any(T .gt. 0._ireallut)) ierr = ierr + 17 ! if I am a dst, nobody should get any radiation anyway
            end if
          end associate
        class is (t_optprop_LUT_rectilinear_wedge_5_8)
          associate (lsrc => lsrc5)
            call src_or_dst_by_param_phi_param_theta5(param_phi, param_theta, lsrc)

            if (lsrc(src)) then ! I am a src face
              if (any(lsrc .and. T .gt. 0._ireallut)) ierr = ierr + 13 ! all srcs cannot have any incoming energy
            else
              if (any(T .gt. 0._ireallut)) then
                ierr = ierr + 17 ! if I am a dst, nobody should get any radiation anyway
              end if
            end if
          end associate
        end select
      end if
    end subroutine
  end subroutine

  subroutine LUT_bmc_wrapper(OPP, config, index_1d, src, dir, comm, S_diff, T_dir, S_tol, T_tol)
    class(t_optprop_LUT) :: OPP
    type(t_op_config), intent(in) :: config
    integer(iintegers), intent(in) :: index_1d
    integer(iintegers), intent(in) :: src
    logical, intent(in) :: dir
    integer(mpiint), intent(in) :: comm

    real(irealLUT), intent(out) :: S_diff(OPP%diff_streams), T_dir(OPP%dir_streams)
    real(irealLUT), intent(out) :: S_tol(OPP%diff_streams), T_tol(OPP%dir_streams)

    real(ireal_dp), allocatable :: vertices(:)
    real(irealLUT) :: tauz, w0, g, phi, theta
    real(ireal_params) :: param_phi, param_theta

    logical :: lvalid
    integer(mpiint) :: ierr, ierr2

    call OPP%LUT_bmc_wrapper_determine_sample_pts(config, index_1d, dir, &
                                                  vertices, tauz, w0, g, phi, theta, lvalid)

    if (.not. lvalid) then
      if (.not. dir) call CHKERR(1_mpiint, 'all diff2diff coeffs should be valid?')
      !if(ldebug) print *,'have an invalid coeff here, skipping computations', index_1d
      S_diff = 0
      T_dir = 0
      S_tol = 0
      T_tol = 0
      return
    end if

    call OPP%bmc_wrapper(src, vertices, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)

    call LUT_bmc_wrapper_validate(OPP, config, tauz, w0, index_1d, src, dir, T_dir, S_diff, ierr)
    if (ierr .ne. 0) then
      ierr2 = ierr
      call param_phi_param_theta_from_phi_and_theta_withcoords( &
        real(vertices, ireal_params), &
        real(deg2rad(phi), ireal_params), &
        real(deg2rad(theta), ireal_params), &
        param_phi, param_theta, ierr)
      print *, 'LUTBMC :: calling bmc_get_coeff src', src, 'tauz', tauz, w0, g, 'angles', phi, theta, &
        ':ind1d', index_1d, ': T', T_dir, ':: S ', S_diff, ': verts', vertices, &
        'param phi/theta_afterwards', param_phi, param_theta
      call CHKERR(ierr2, 'found bad results for a given geometry')
    end if

    !if(ldebug) print *,'LUTBMC :: calling bmc_get_coeff src', src, &
    !  'tauz',tauz,w0,g,'angles',phi,theta, ':ind1d',index_1d, ':', &
    !  T_dir, ':(', T_tol, ') //', S_diff, ':(', S_tol, ')'
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

  subroutine set_parameter_space(OPP)
    class(t_optprop_lut) :: OPP
    integer(mpiint) :: ierr

    select type (OPP)
    class is (t_optprop_LUT_1_2)
      call set_op_param_space(OPP, 'LUT_1_2', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_8_10)
      call set_op_param_space(OPP, 'LUT_8_10', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_8_12)
      call set_op_param_space(OPP, 'LUT_8_12', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_8_16)
      call set_op_param_space(OPP, 'LUT_8_16', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_8_18)
      call set_op_param_space(OPP, 'LUT_8_18', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_3_10)
      call set_op_param_space(OPP, 'LUT_3_10', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_3_10_for_ANN)
      call set_op_param_space(OPP, 'LUT_3_10_for_ANN', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_3_16)
      call set_op_param_space(OPP, 'LUT_3_16', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_3_6)
      call set_op_param_space(OPP, 'LUT_3_6', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_wedge_5_8)
      call set_op_param_space(OPP, 'LUT_wedge_5_8', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_rectilinear_wedge_5_8)
      call set_op_param_space(OPP, 'LUT_rectilinear_wedge_5_8', ierr); call CHKERR(ierr)

    class is (t_optprop_LUT_wedge_18_8)
      call set_op_param_space(OPP, 'LUT_wedge_18_8', ierr); call CHKERR(ierr)

    class default
      call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop_LUT object!')
    end select
  end subroutine

  subroutine scatter_LUTtables(OPP, comm)
    integer(mpiint), intent(in) :: comm
    class(t_optprop_LUT) :: OPP

    integer(mpiint) :: myid, ierr
    real(irealLUT), pointer :: mmap_ptr(:, :)

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)

    if (luse_memory_map) then
      mmap_ptr => null()
      if (associated(OPP%Sdiff%c)) then
        call arr_to_mmap(comm, trim(OPP%Sdiff%table_name_c(1))//'.Sdiff.mmap', mmap_ptr, ierr, OPP%Sdiff%c)
        deallocate (OPP%Sdiff%c)
      else
        call arr_to_mmap(comm, trim(OPP%Sdiff%table_name_c(1))//'.Sdiff.mmap', mmap_ptr, ierr)
      end if
      OPP%Sdiff%c => mmap_ptr

      call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
                  '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
                  '  You probably need to run the createLUT program first '//new_line('')// &
                  '  or run with the additional option:'//new_line('')// &
                  '    -skip_load_LUT no')

      mmap_ptr => null()
      if (associated(OPP%Sdir%c)) then
        call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr, OPP%Sdir%c)
        deallocate (OPP%Sdir%c)
      else
        call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr)
      end if
      OPP%Sdir%c => mmap_ptr

      call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
                  '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
                  '  You probably need to run the createLUT program first '//new_line('')// &
                  '  or run with the additional option:'//new_line('')// &
                  '    -skip_load_LUT no')

      mmap_ptr => null()
      if (associated(OPP%Tdir%c)) then
        call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr, OPP%Tdir%c)
        deallocate (OPP%Tdir%c)
      else
        call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr)
      end if
      OPP%Tdir%c => mmap_ptr

      call CHKERR(ierr, 'Could not generate mmap for LUT. '//new_line('')// &
                  '  If the binary dump (*.mmap) file does not yet exist, '//new_line('')// &
                  '  You probably need to run the createLUT program first '//new_line('')// &
                  '  or run with the additional option:'//new_line('')// &
                  '    -skip_load_LUT no')

    else
      if (myid .eq. 0) then
        if (.not. associated(OPP%Sdir%c)) then
          call CHKERR(1_mpiint, 'LUT data is not loaded. '//new_line('')// &
                      '  Somehow we ended up in a situation where we would like to distribute'//new_line('')// &
                      '  the LUT`s via MPI but rank 0 does not have the info. '//new_line('')// &
                      '  Did you use the option: -skip_load_LUT ? '//new_line('')// &
                      '  Maybe try to set it to -skip_load_LUT no')
        end if
      end if

      if (mpi_logical_or(comm,.not. associated(OPP%Sdir%c))) &
        call imp_bcast(comm, OPP%Sdir%c, 0_mpiint)  ! DIRECT 2 DIRECT

      if (mpi_logical_or(comm,.not. associated(OPP%Tdir%c))) &
        call imp_bcast(comm, OPP%Tdir%c, 0_mpiint)  ! DIRECT 2 DIFFUSE

      if (mpi_logical_or(comm,.not. associated(OPP%Sdiff%c))) &
        call imp_bcast(comm, OPP%Sdiff%c, 0_mpiint)

    end if
  end subroutine

  subroutine print_configs(OPP)
    class(t_optprop_LUT) :: OPP
    integer(mpiint) :: myid, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)
    if (myid .eq. 0) then
      print *, 'Diffuse LUT config:'
      call print_op_config(OPP%diffconfig)
      print *, 'Direct LUT config:'
      call print_op_config(OPP%dirconfig)
      print *, '----------------------'
    end if
  end subroutine

  subroutine LUT_get_dir2dir(OPP, sample_pts, C)
    class(t_optprop_LUT), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(OPP%dir_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if (ldebug_optprop) then
      if (size(sample_pts) .ne. size(OPP%dirconfig%dims)) then
        call print_op_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
                    //toStr(size(sample_pts, kind=iintegers))//'/'//toStr(size(OPP%dirconfig%dims)))
      end if
      call check_if_samplepts_in_bounds(sample_pts, OPP%dirconfig)
    end if

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    end do

    select case (interp_mode)
    case (1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Tdir%c(:, ind1d)
    case (2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Tdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//toStr(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if (ldebug_optprop) then
      !Check for energy conservation:
      iierr = 0
      do src = 1, OPP%dir_streams
        norm = sum(C(src:size(C):OPP%dir_streams))
        if (real(norm) .gt. one + 1e-5_ireallut) iierr = iierr + 1
      end do
      if (iierr .ne. 0) then
        print *, 'Error in dir2dir coeffs :: ierr', iierr, size(C), OPP%dir_streams, '::', C
        do src = 1, OPP%dir_streams
          print *, 'SUM dir2dir coeff for src ', src, ' :: sum ', sum(C(src:size(C):OPP%dir_streams)), &
            ' :: coeff', C(src:size(C):OPP%dir_streams)
        end do
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      end if
    end if
    !call CHKERR(1_mpiint, 'DEBUG')
  end subroutine

  subroutine LUT_get_dir2diff(OPP, sample_pts, C)
    class(t_optprop_LUT), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(OPP%dir_streams*OPP%diff_streams)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if (ldebug_optprop) then
      if (size(sample_pts) .ne. size(OPP%dirconfig%dims)) then
        call print_op_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
                    //toStr(size(sample_pts, kind=iintegers))//'/'//toStr(size(OPP%dirconfig%dims)))
      end if
      call check_if_samplepts_in_bounds(sample_pts, OPP%dirconfig)
    end if

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    end do

    select case (interp_mode)
    case (1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Sdir%c(:, ind1d)
    case (2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Sdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//toStr(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if (ldebug_optprop) then
      !Check for energy conservation:
      iierr = 0
      do src = 1, OPP%diff_streams
        norm = sum(C(src:size(C):OPP%dir_streams))
        if (real(norm) .gt. one + 1e-5_ireallut) iierr = iierr + 1
      end do
      if (iierr .ne. 0) then
        do src = 1, OPP%diff_streams
          print *, 'SUM dir2diff coeff for src ', src, ' :: sum ', sum(C(src:size(C):OPP%dir_streams)), &
            ' :: coeff', C(src:size(C):OPP%dir_streams)
        end do
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      end if
    end if
  end subroutine

  subroutine LUT_get_diff2diff(OPP, sample_pts, C)
    class(t_optprop_LUT), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(OPP%diff_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(irealLUT) :: pti_buffer(LUT_MAX_DIM)
    real(irealLUT) :: norm

    if (ldebug_optprop) then
      if (size(sample_pts) .ne. size(OPP%diffconfig%dims)) then
        call print_op_config(OPP%diffconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
                    //toStr(size(sample_pts, kind=iintegers))//'/'//toStr(size(OPP%diffconfig%dims, kind=iintegers)))
      end if
      call check_if_samplepts_in_bounds(sample_pts, OPP%diffconfig)
    end if

    do kdim = 1, size(sample_pts)
      pti_buffer(kdim) = find_real_location(OPP%diffconfig%dims(kdim)%v, sample_pts(kdim))
    end do
    !print *,'LUT_get_diff2diff sample', sample_pts, 'idx', pti_buffer(1:size(sample_pts)), &
    !  & 'i1d', ind_nd_to_1d(OPP%diffconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))

    select case (interp_mode)
    case (1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%diffconfig%offsets, nint(pti_buffer(1:size(sample_pts)), kind=iintegers))
      C = OPP%Sdiff%c(:, ind1d)
    case (2)
      call interp_vec_simplex_nd(pti_buffer(1:size(sample_pts)), OPP%Sdiff%c, OPP%diffconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//toStr(interp_mode)//' not implemented yet! please choose something else!')
    end select

    if (ldebug_optprop) then
      !Check for energy conservation:
      iierr = 0
      do src = 1, OPP%diff_streams
        norm = sum(C(src:size(C):OPP%diff_streams))
        if (norm .gt. one + 1e-5_ireallut) iierr = iierr + 1
      end do
      if (iierr .ne. 0) then
        do src = 1, OPP%diff_streams
          print *, 'SUM diff2diff coeff for src ', src, ' :: sum ', sum(C(src:size(C):OPP%diff_streams)), &
            ' :: coeff', C(src:size(C):OPP%diff_streams)
        end do
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      end if
    end if
  end subroutine

end module
