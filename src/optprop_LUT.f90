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

#include "petsc/finclude/petsc.h"
  use petsc

  use mpi!, only: MPI_BCAST,MPI_LAND,MPI_LOR

  use m_helper_functions, only : approx,  &
    rel_approx, imp_bcast, get_arg, ftoa, &
    mpi_logical_and, mpi_logical_or,      &
    search_sorted_bisection, CHKERR, itoa,&
    triangle_area_by_vertices,            &
    ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets

  use m_data_parameters, only : ireals, iintegers, irealLUT, &
    one, zero, i0, i1, i2, i3, i10, mpiint, nil, inil,       &
    imp_iinteger, imp_ireals, imp_irealLUT, imp_logical,     &
    default_str_len

  use m_optprop_parameters, only:         &
    ldebug_optprop, lut_basename,         &
    LUT_dump_interval, LUT_max_create_jobtime, &
    interp_mode_pprts,                    &
    interp_mode_wedge_5_8,                &
    ldelta_scale,delta_scale_truncate,    &
    stddev_atol, stddev_rtol,             &
    wedge_sphere_radius,                  &
    preset_g2, preset_g3, preset_g4,      &
    preset_aspect5,                       &
    preset_aspect7,                       &
    preset_aspect13,                      &
    preset_aspect23,                      &
    preset_w010,                          &
    preset_w020,                          &
    preset_tau15,                         &
    preset_tau20,                         &
    preset_tau31

  use m_boxmc, only: t_boxmc, &
    t_boxmc_1_2, &
    t_boxmc_3_6, t_boxmc_3_10, t_boxmc_3_16, &
    t_boxmc_8_10, t_boxmc_8_12, t_boxmc_8_16, &
    t_boxmc_8_18, &
    t_boxmc_wedge_5_8
  use m_tenstream_interpolation, only: interp_4d, interp_vec_simplex_nd
  use m_netcdfio

  use m_mmap, only : arr_to_mmap, munmap_mmap_ptr

  use m_boxmc_geometry, only : setup_default_wedge_geometry, setup_default_unit_cube_geometry

  implicit none

  private
  public :: t_optprop_LUT, t_optprop_LUT_1_2, &
    t_optprop_LUT_3_6, t_optprop_LUT_3_10, t_optprop_LUT_3_16, &
    t_optprop_LUT_8_10, t_optprop_LUT_8_12, t_optprop_LUT_8_16, t_optprop_LUT_8_18, &
    t_optprop_LUT_wedge_5_8, &
    find_lut_dim_by_name
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
    real(irealLUT), pointer :: c(:,:) => NULL() ! depending on config has Ndim_1*Ndim_2*etc. many entries
    real(irealLUT), allocatable :: stddev_tol(:,:)
    character(default_str_len), allocatable :: table_name_c(:)
    character(default_str_len), allocatable :: table_name_tol(:)
  end type

  type,abstract :: t_optprop_LUT
    class(t_boxmc), allocatable :: bmc
    type(t_table), allocatable :: Sdiff, Sdir, Tdir
    type(t_LUT_config), allocatable :: dirconfig, diffconfig
    integer(iintegers) :: dir_streams = inil, diff_streams = inil
    integer(iintegers) :: interp_mode
    logical :: LUT_initialized=.False., optprop_LUT_debug=ldebug_optprop
    character(default_str_len) :: lutbasename

    contains
      procedure :: init
      procedure :: destroy
      procedure :: LUT_get_dir2dir
      procedure :: LUT_get_dir2diff
      procedure :: LUT_get_diff2diff
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

  logical, parameter :: ldebug=.True.

contains

  subroutine init(OPP, comm)
      class(t_optprop_LUT) :: OPP
      integer(mpiint) ,intent(in) :: comm

      integer(mpiint) :: comm_size, myid

      if(OPP%LUT_initialized) return

      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Initializing LUT`s...'

      if(.not.allocated(OPP%bmc)) then
        select type (OPP)
          class is (t_optprop_LUT_1_2)
            OPP%dir_streams  =  1
            OPP%diff_streams =  2
            OPP%lutbasename=trim(lut_basename)//'_1_2.'
            allocate(t_boxmc_1_2::OPP%bmc)

          class is (t_optprop_LUT_3_6)
            OPP%dir_streams  = 3
            OPP%diff_streams = 6
            OPP%lutbasename=trim(lut_basename)//'_3_6.'
            allocate(t_boxmc_3_6::OPP%bmc)

          class is (t_optprop_LUT_3_10)
            OPP%dir_streams  = 3
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)//'_3_10.'
            allocate(t_boxmc_3_10::OPP%bmc)

          class is (t_optprop_LUT_3_16)
            OPP%dir_streams  = 3
            OPP%diff_streams = 16
            OPP%lutbasename=trim(lut_basename)//'_3_16.'
            allocate(t_boxmc_3_16::OPP%bmc)

          class is (t_optprop_LUT_8_10)
            OPP%dir_streams  = 8
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)//'_8_10.'
            allocate(t_boxmc_8_10::OPP%bmc)

          class is (t_optprop_LUT_8_12)
            OPP%dir_streams  = 8
            OPP%diff_streams = 12
            OPP%lutbasename=trim(lut_basename)//'_8_12.'
            allocate(t_boxmc_8_12::OPP%bmc)

          class is (t_optprop_LUT_8_16)
            OPP%dir_streams  = 8
            OPP%diff_streams = 16
            OPP%lutbasename=trim(lut_basename)//'_8_16.'
            allocate(t_boxmc_8_16::OPP%bmc)

          class is (t_optprop_LUT_8_18)
            OPP%dir_streams  = 8
            OPP%diff_streams = 18
            OPP%lutbasename=trim(lut_basename)//'_8_18.'
            allocate(t_boxmc_8_18::OPP%bmc)

          class is (t_optprop_LUT_wedge_5_8)
            OPP%dir_streams  = 5
            OPP%diff_streams = 8
            OPP%lutbasename=trim(lut_basename)//'_wedge_5_8.Rsphere'//itoa(int(wedge_sphere_radius))//'.'
            allocate(t_boxmc_wedge_5_8::OPP%bmc)

          class default
            stop 'initialize LUT: unexpected type for optprop_LUT object!'
        end select

        call OPP%bmc%init(comm)
      endif

      call OPP%set_parameter_space()

      call OPP%loadLUT_dir(comm)
      call OPP%loadLUT_diff(comm)

      call OPP%scatter_LUTtables(comm)

      OPP%LUT_initialized=.True.
      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Initializing LUT`s... finished'
  end subroutine

  subroutine destroy(OPP)
      use m_optprop_parameters, only: luse_memory_map
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

  integer(mpiint) :: ierr
  istat = 0

  call ncload(table%table_name_tol, table%stddev_tol, ierr); istat = istat + ierr
  if(ierr.eq.0) then ! we were able to load stddev but still have to check if they all have a good enough noise...
    if(any(table%stddev_tol.gt.stddev_atol+10*epsilon(stddev_atol))) then
      istat = istat + 1
      print *,'coefficients not good enough', maxval(table%stddev_tol),'should be less than', stddev_atol+10*epsilon(stddev_atol)
    endif
  endif

  call ncload(table%table_name_c, table%c, ierr); istat = istat + ierr
  if(ierr.eq.0) then ! we were able to load coeffs but still have to check if they are all ok...
    if(any( table%c.gt.one ).or.any(table%c.lt.zero) ) then
      istat = istat + 2
      print *,'coefficients not good enough', minval(table%c), maxval(table%c),'should be between [0,1]'
    endif
  endif
  if(istat.ne.0) print *,'Test if coeffs are good results in:', istat
end subroutine

subroutine loadLUT_diff(OPP, comm)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid, ierr
    logical :: lskip_load_LUT, lflg

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(allocated(OPP%Sdiff)) return ! already loaded
    allocate(OPP%Sdiff)
    OPP%Sdiff%c => NULL()
    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('diffuse', OPP%diffconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'diffuse'

    str(3) = 'Stol'; allocate(OPP%Sdiff%table_name_tol(size(str)), source=str)
    str(3) = 'S'   ; allocate(OPP%Sdiff%table_name_c(size(str)), source=str)

    lskip_load_LUT = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-skip_load_LUT', lskip_load_LUT, lflg, ierr); call CHKERR(ierr)
    if(lskip_load_LUT) return

    if(myid.eq.0) then
      call load_table_from_netcdf(OPP%Sdiff, iierr); errcnt = errcnt + iierr
      call write_pspace(str(1), OPP%diffconfig)
    endif

    call mpi_bcast(errcnt, 1_mpiint, imp_iinteger, 0_mpiint ,comm ,mpierr); call CHKERR(mpierr)

    if(errcnt.ne.0) then ! something went wrong loading the LUT
      call OPP%createLUT(comm, OPP%diffconfig, OPP%Sdiff)
    endif

    if(allocated(OPP%Sdiff%stddev_tol)) deallocate(OPP%Sdiff%stddev_tol)
end subroutine

subroutine loadLUT_dir(OPP, comm)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt
    character(default_str_len) :: descr, str(3)

    integer(mpiint) :: comm_size, myid, ierr
    logical :: lskip_load_LUT, lflg

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(allocated(OPP%Sdir).and.allocated(OPP%Tdir)) return ! already loaded
    allocate(OPP%Sdir)
    allocate(OPP%Tdir)
    OPP%Sdir%c => NULL()
    OPP%Tdir%c => NULL()

    errcnt = 0

    ! Set filename of LUT
    descr = gen_lut_basename('direct', OPP%dirconfig)

    str(1) = trim(OPP%lutbasename)//trim(descr)//'.nc'
    str(2) = 'direct'

    str(3) = 'Stol'; allocate(OPP%Sdir%table_name_tol(size(str)), source=str)
    str(3) = 'S'   ; allocate(OPP%Sdir%table_name_c(size(str)), source=str)

    str(3) = 'Ttol'; allocate(OPP%Tdir%table_name_tol(size(str)), source=str)
    str(3) = 'T'   ; allocate(OPP%Tdir%table_name_c(size(str)), source=str)

    lskip_load_LUT = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-skip_load_LUT', lskip_load_LUT, lflg, ierr); call CHKERR(ierr)
    if(lskip_load_LUT) return

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
        print *, kdim, trim(groups(1)), trim(groups(3)), ':existing', existing_values, ':new', config%dims(kdim)%v
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
    real(irealLUT) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(irealLUT) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    integer(mpiint), parameter :: READYMSG=1,HAVERESULTSMSG=2, WORKMSG=3, FINALIZEMSG=4, RESULTMSG=5

    integer(iintegers) :: idummy, lutindex
    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(myid.eq.0) then
      call prepare_table_space(OPP, config, S, T)
    endif

    if(myid.le.0 .and. comm_size.le.1) &
      stop 'At the moment creation of direct Lookuptable needs at least two mpi-ranks to work... please run with more ranks.'

    if(myid.eq.0) then
      call master(S, T)
      print *,'done calculating direct coefficients'
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

        logical :: ldoneS(OPP%diff_streams), ldoneT(OPP%dir_streams)

        real :: starttime, lastsavetime, now
        call cpu_time(starttime)
        lastsavetime = starttime

        finalizedworkers=0
        if(present(T)) then
          Nsrc = OPP%dir_streams
        else
          Nsrc = OPP%diff_streams
        endif
        total_size = size(S%c, dim=2) * Nsrc ! we have Ndims * Nsrc work packages

        cnt=1
        do
          ! Check if we already calculated the coefficients
          if(cnt.le.total_size) then
            isrc = modulo(cnt-1, Nsrc)+1
            lutindex = (cnt-1) / Nsrc +1

            do idst = 1,OPP%diff_streams
                ind = (idst-1) * Nsrc + isrc
                ldoneS(idst) = ( ( S%c         ( ind, lutindex ).ge.zero)            &
                           .and. ( S%c         ( ind, lutindex ).le.one )            &
                           .and. ( S%stddev_tol( ind, lutindex ).le.stddev_atol ) )
            enddo
            if(present(T)) then
              do idst = 1,OPP%dir_streams
                ind = (idst-1)*OPP%dir_streams + isrc
                ldoneT(idst) = ( ( T%c         ( ind, lutindex ).ge.zero)            &
                  .and. ( T%c         ( ind, lutindex ).le.one )            &
                  .and. ( T%stddev_tol( ind, lutindex ).le.stddev_atol ) )
              enddo
            else
              ldoneT = .True.
            endif

            if( all(ldoneS) .and. all(ldoneT) ) then
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
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), READYMSG, comm, status, mpierr) ; call CHKERR(mpierr)

              if(cnt.le.total_size) then ! we got something to do for a worker -- send him...
                isrc = modulo(cnt-1, Nsrc) +1
                lutindex = (cnt-1) / Nsrc +1
                call mpi_send(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(isrc, 1_mpiint, imp_iinteger, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(present(T), 1_mpiint, imp_logical, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)

              else ! no more work to do... tell the worker to quit
                call mpi_send(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
              endif
              cnt = cnt+1

            case(HAVERESULTSMSG)
              call mpi_recv(lutindex, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(isrc, 1_mpiint, imp_iinteger, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_diff, size(S_diff), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_tol , size(S_tol ), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_dir , size(T_dir ), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_tol , size(T_tol ), imp_irealLUT, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)

              ! Sort coefficients into destination ordering and put em in LUT
              do idst = 1, OPP%diff_streams
                ind = (idst-1) * Nsrc + isrc
                S%c         (ind, lutindex) = S_diff(idst)
                S%stddev_tol(ind, lutindex) = S_tol (idst)
              enddo
              if(present(T)) then
                do idst = 1, OPP%dir_streams
                  ind = (idst-1)*OPP%dir_streams + isrc
                  T%c         (ind, lutindex) = T_dir (idst)
                  T%stddev_tol(ind, lutindex) = T_tol (idst)
                enddo
              endif

              if( mod((lutindex-1)*Nsrc+isrc-1, max(i1, total_size/1000_iintegers)).eq.0 ) & !every .1 percent report status
                  print *,'Calculated LUT...', lutindex, isrc, ((lutindex-1)*Nsrc+isrc-1)*100._irealLUT/total_size,'%'

              call cpu_time(now)
              if( (now-lastsavetime).gt.LUT_dump_interval .or. (now-starttime).gt.LUT_max_create_jobtime ) then !every 30 minutes wall clock time, dump the LUT.
                print *,'Dumping LUT after ',(now-lastsavetime)/60,'minutes'
                if(present(T)) then
                  print *,'Writing table to file... ', T%table_name_c
                  call ncwrite(T%table_name_c, T%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
                  print *,'Writing table to file... ', T%table_name_tol
                  call ncwrite(T%table_name_tol, T%stddev_tol, iierr); call CHKERR(iierr, 'Could not write Table to file')
                endif
                print *,'Writing table to file... ', S%table_name_c
                call ncwrite(S%table_name_c, S%c, iierr); call CHKERR(iierr, 'Could not write Table to file')
                print *,'Writing table to file... ', S%table_name_tol
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr); call CHKERR(iierr, 'Could not write Table to file')
                print *,'done writing!',iierr
                lastsavetime = now ! reset the countdown
                if((now-starttime).gt.LUT_max_create_jobtime) then
                  call CHKERR(int(now-starttime, mpiint),'Maximum duration of create_LUT reached ... exiting')
                endif
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
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
          integer(iintegers) :: isrc
          logical :: ldir

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
                call mpi_recv( lutindex, 1_mpiint, imp_iinteger, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( isrc, 1_mpiint, imp_iinteger, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( ldir, 1_mpiint, imp_logical, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)

                call OPP%LUT_bmc_wrapper(config, lutindex, isrc, ldir, &
                  mpi_comm_self, S_diff, T_dir, S_tol, T_tol)

                if ( maxval(S_tol).gt.stddev_atol+sqrt(epsilon(stddev_atol)) &
                .or. maxval(T_tol).gt.stddev_atol+sqrt(epsilon(stddev_atol))) then
                  print *,'Ttol', T_tol
                  print *,'Stol', S_tol
                  call CHKERR(1_mpiint, 'BOXMC violated stddev constraints!')
                endif

                !print *,'Computed isrc',isrc,'aspect',aspect_zx,'tau',tau_z, w0, g,':', phi, theta
                !print *,myid,'Computed values for ',lutindex, isrc, ldir

                call mpi_send(lutindex , 1_mpiint     , imp_iinteger  , status(MPI_SOURCE) , HAVERESULTSMSG , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(isrc , 1_mpiint     , imp_iinteger  , status(MPI_SOURCE) , HAVERESULTSMSG , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_diff    , size(S_diff) , imp_irealLUT , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_tol     , size(S_tol ) , imp_irealLUT , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(T_dir     , size(T_dir ) , imp_irealLUT , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(T_tol     , size(T_tol ) , imp_irealLUT , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)

                call mpi_send(-i1       , 1_mpiint     , imp_iinteger  , 0_mpiint           , READYMSG       , comm , mpierr); call CHKERR(mpierr)

              case(FINALIZEMSG)
                call mpi_recv(idummy, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_send(-i1, 1_mpiint, imp_iinteger, 0_mpiint, FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
                exit

              end select

            endif !gotmsg
          enddo
      end subroutine
end subroutine createLUT

subroutine prepare_table_space(OPP, config, S, T)
  class(t_optprop_LUT) :: OPP
  type(t_lut_config), intent(in) :: config
  type(t_table),intent(inout) :: S
  type(t_table),intent(inout), optional :: T

  real(irealLUT) :: bytesize, entries
  integer(iintegers) :: errcnt

  if(present(T)) then
    entries = real( &
      OPP%diff_streams*OPP%dir_streams * product(config%dims(:)%N) + &
      OPP%dir_streams**2 * product(config%dims(:)%N) + &
      OPP%diff_streams*OPP%dir_streams * product(config%dims(:)%N) + &
      OPP%dir_streams**2 * product(config%dims(:)%N), irealLUT)
  else
    entries = real( &
      OPP%diff_streams**2 * product(config%dims(:)%N) + &
      OPP%diff_streams**2 * product(config%dims(:)%N), irealLUT)
  endif
  bytesize = real(C_SIZEOF(bytesize), irealLUT) * entries

  print *,'Allocating Space for LUTs ( '//ftoa(bytesize/1024**3)//' Gb)'
  errcnt = 0
  if(present(T)) then
    if(.not.associated(S%c)) allocate(S%c(OPP%diff_streams*OPP%dir_streams, product(config%dims(:)%N)), source=real(nil,irealLUT))
    if(.not.associated(T%c)) allocate(T%c(OPP%dir_streams**2              , product(config%dims(:)%N)), source=real(nil,irealLUT))
    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(OPP%diff_streams*OPP%dir_streams, product(config%dims(:)%N)), source=1e8_irealLUT)
    if(.not.allocated (T%stddev_tol)) allocate(T%stddev_tol(OPP%dir_streams**2,               product(config%dims(:)%N)), source=1e8_irealLUT)
  else
    if(.not.associated(S%c)) allocate(S%c(OPP%diff_streams**2, product(config%dims(:)%N)), source=real(nil,irealLUT))
    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(OPP%diff_streams**2, product(config%dims(:)%N)), source=1e8_irealLUT)
  endif
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
  real(ireals), intent(out) :: sample_pnt
  integer(mpiint), intent(out) :: ierr

  integer(iintegers) :: kdim, nd_indices(size(config%dims))

  call ind_1d_to_nd(config%offsets, index_1d, nd_indices)

  kdim = find_lut_dim_by_name(config, trim(dimname))
  if(kdim.lt.i1) then ! could not find the corresponding dimension
    ierr = 1
    sample_pnt = nil
    return
  endif
  if(nd_indices(kdim).gt.size(config%dims(kdim)%v)) then
    print *,index_1d,'nd_indices', nd_indices
    call CHKERR(1_mpiint, 'wrong indices in kdim')
  endif
  sample_pnt = config%dims(kdim)%v(nd_indices(kdim))
  ierr = 0
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

    real(ireals) :: aspect_zx, aspect_zy, tauz, w0, g, phi, theta
    real(ireals), dimension(2) :: wedge_C
    real(ireals), allocatable :: vertices(:)
    integer(mpiint) :: ierr

    call get_sample_pnt_by_name_and_index(config, 'aspect_zx', index_1d, aspect_zx, ierr); call CHKERR(ierr, 'aspect_zx has to be present')
    call get_sample_pnt_by_name_and_index(config, 'tau', index_1d, tauz, ierr); call CHKERR(ierr, 'tauz has to be present')
    call get_sample_pnt_by_name_and_index(config, 'w0', index_1d, w0, ierr); call CHKERR(ierr, 'w0 has to be present')
    if(dir) then
      call get_sample_pnt_by_name_and_index(config, 'phi', index_1d, phi, ierr); call CHKERR(ierr, 'phi has to be present for direct calculations')
      call get_sample_pnt_by_name_and_index(config, 'theta', index_1d, theta, ierr); call CHKERR(ierr, 'theta has to be present for direct calculations')
    endif

    call get_sample_pnt_by_name_and_index(config, 'g', index_1d, g, ierr)
    if(ierr.ne.0) then
      g = zero ! defaults to isotropic scattering
    endif

    call get_sample_pnt_by_name_and_index(config, 'aspect_zy', index_1d, aspect_zy, ierr)
    if(ierr.ne.0) then
      aspect_zy = aspect_zx ! set dy = dx
    endif

    ! First define Default vertices for Cube Geometry
    select type(OPP)
      class is (t_optprop_LUT_1_2)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_3_6)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_3_10)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_3_16)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_8_10)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_8_12)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_8_16)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_8_18)
        call setup_default_unit_cube_geometry(one, aspect_zx/aspect_zy, aspect_zx, vertices)

      class is (t_optprop_LUT_wedge_5_8)
        call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cx', index_1d, wedge_C(1), ierr); call CHKERR(ierr, 'wedge_coord_Cx has to be present for wedge calculations')
        call get_sample_pnt_by_name_and_index(config, 'wedge_coord_Cy', index_1d, wedge_C(2), ierr); call CHKERR(ierr, 'wedge_coord_Cy has to be present for wedge calculations')
        call setup_default_wedge_geometry([zero, zero], [one, zero], &
          wedge_C, aspect_zx, vertices, &
          sphere_radius=real(wedge_sphere_radius, ireals))

      class default
        call CHKERR(1_mpiint, 'unexpected type for optprop_LUT object!')
    end select

    call bmc_wrapper(OPP, src, vertices, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)

end subroutine

subroutine bmc_wrapper(OPP, src, vertices, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol, inp_atol, inp_rtol)
    class(t_optprop_LUT) :: OPP
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: dir
    integer(mpiint),intent(in) :: comm
    real(ireals), intent(in) :: tauz, w0, g, phi, theta
    real(ireals) :: vertices(:)

    real(irealLUT),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(irealLUT),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)
    real(ireals),intent(in),optional :: inp_atol, inp_rtol
    real(ireals) :: rS_diff(OPP%diff_streams),rT_dir(OPP%dir_streams)
    real(ireals) :: rS_tol (OPP%diff_streams),rT_tol(OPP%dir_streams)

    real(ireals) :: bg(3), dz, atol, rtol

    atol = get_arg(stddev_atol-epsilon(atol)*10, inp_atol)
    rtol = get_arg(stddev_rtol-epsilon(rtol)*10, inp_rtol)

    dz = vertices(size(vertices))

    bg(1) = tauz / dz * (one-w0)
    bg(2) = tauz / dz * w0
    bg(3) = g

    S_diff=nil
    T_dir=nil
    S_tol=nil
    T_tol=nil

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
    call OPP%bmc%get_coeff(comm, bg, src, &
      dir, phi, theta, vertices,   &
      rS_diff, rT_dir, rS_tol, rT_tol, &
      inp_atol=atol, inp_rtol=rtol )
    S_diff = real(rS_diff, irealLUT)
    T_dir  = real(rT_dir, irealLUT)
    S_tol  = real(rS_tol, irealLUT)
    T_tol  = real(rT_tol, irealLUT)
    !print *,'BMC :: dir',T_dir,'diff',S_diff
end subroutine

function lin_index_to_param(idx,rng,N)
    real(irealLUT) :: lin_index_to_param
    real(irealLUT),intent(in) :: idx,rng(2)
    integer(iintegers),intent(in) :: N
    if(N.gt.i1) then
      lin_index_to_param = rng(1) + (idx-1) * ( rng(2)-rng(1) ) / (N-1)
    else
      lin_index_to_param = rng(1)
    endif
end function

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
      lut_dim%v(k) = lin_index_to_param(real(k, irealLUT), vrange, N)
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

    select type(OPP)
      class is (t_optprop_LUT_1_2)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i1, OPP%dirconfig%dims(5), vrange=real([0], irealLUT)) ! azimithally average in 1D
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_8_10)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_8_12)
          OPP%interp_mode = interp_mode_pprts
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
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
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
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_8_16)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_8_18)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_3_10)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_3_16)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_3_6)
          OPP%interp_mode = interp_mode_pprts
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%dirconfig%dims(4), preset=preset_g4)
          call populate_LUT_dim('phi',       i10, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
          call populate_LUT_dim('theta',     i10, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w020,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w020)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('g',         size(preset_g4,kind=iintegers), OPP%diffconfig%dims(4), preset=preset_g4)

      class is (t_optprop_LUT_wedge_5_8)
          OPP%interp_mode = interp_mode_wedge_5_8
          allocate(OPP%dirconfig%dims(7))
          call populate_LUT_dim('tau',       size(preset_tau15,kind=iintegers), OPP%dirconfig%dims(1), preset=preset_tau15)
          call populate_LUT_dim('w0',        size(preset_w010,kind=iintegers), OPP%dirconfig%dims(2), preset=preset_w010)
          call populate_LUT_dim('aspect_zx', size(preset_aspect7,kind=iintegers), OPP%dirconfig%dims(3), preset=preset_aspect7)
          call populate_LUT_dim('wedge_coord_Cx', 5_iintegers, OPP%dirconfig%dims(4), vrange=real([.35,.65], irealLUT))
          call populate_LUT_dim('wedge_coord_Cy', 10_iintegers, OPP%dirconfig%dims(5), vrange=real([.8, .9485], irealLUT))
          call populate_LUT_dim('phi',       71_iintegers, OPP%dirconfig%dims(6), vrange=real([-70,70], irealLUT))
          call populate_LUT_dim('theta',     10_iintegers, OPP%dirconfig%dims(7), vrange=real([0,90], irealLUT))

          !call populate_LUT_dim('tau',       i2, OPP%dirconfig%dims(1), vrange=real([1e-3,1.], irealLUT))
          !call populate_LUT_dim('w0',        i2, OPP%dirconfig%dims(2), vrange=real([.0,.99999], irealLUT))
          !call populate_LUT_dim('aspect_zx', i2, OPP%dirconfig%dims(3), vrange=real([.5,2.], irealLUT))
          !call populate_LUT_dim('wedge_coord_Cx', 5_iintegers, OPP%dirconfig%dims(4), vrange=real([.35,.65], irealLUT))
          !call populate_LUT_dim('wedge_coord_Cy', 5_iintegers, OPP%dirconfig%dims(5), vrange=real([.8, .95], irealLUT))
          !call populate_LUT_dim('phi',       i3, OPP%dirconfig%dims(6), vrange=real([-70,70], irealLUT))
!          call populate_LUT_dim('theta',     i3, OPP%dirconfig%dims(7), vrange=real([0,90], irealLUT))

          allocate(OPP%diffconfig%dims(5))
          call populate_LUT_dim('tau',       size(preset_tau31,kind=iintegers), OPP%diffconfig%dims(1), preset=preset_tau31)
          call populate_LUT_dim('w0',        size(preset_w010,kind=iintegers), OPP%diffconfig%dims(2), preset=preset_w010)
          call populate_LUT_dim('aspect_zx', size(preset_aspect23,kind=iintegers), OPP%diffconfig%dims(3), preset=preset_aspect23)
          call populate_LUT_dim('wedge_coord_Cx', 10_iintegers, OPP%diffconfig%dims(4), vrange=real([.35,.65], irealLUT))
          call populate_LUT_dim('wedge_coord_Cy', 10_iintegers, OPP%diffconfig%dims(5), vrange=real([.8, .95], irealLUT))

          !call populate_LUT_dim('tau',       i2, OPP%diffconfig%dims(1), vrange=real([1e-3,1.], irealLUT))
          !call populate_LUT_dim('w0',        i2, OPP%diffconfig%dims(2), vrange=real([.1,.999], irealLUT))
          !call populate_LUT_dim('aspect_zx', i2, OPP%diffconfig%dims(3), vrange=real([.5,2.], irealLUT))
          !call populate_LUT_dim('wedge_coord_Cx', 2_iintegers, OPP%diffconfig%dims(4), vrange=real([.35,.65], irealLUT))
          !call populate_LUT_dim('wedge_coord_Cy', 2_iintegers, OPP%diffconfig%dims(5), vrange=real([.8, .95], irealLUT))

      class default
        call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop_LUT object!')
    end select

    ! Determine offsets
    allocate(OPP%dirconfig%offsets(size(OPP%dirconfig%dims)))
    call ndarray_offsets(OPP%dirconfig%dims(:)%N, OPP%dirconfig%offsets)

    allocate(OPP%diffconfig%offsets(size(OPP%diffconfig%dims)))
    call ndarray_offsets(OPP%diffconfig%dims(:)%N, OPP%diffconfig%offsets)

    if(ldebug.and..False.) then
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      print *,myid,'set_parameter space dims:', size(OPP%diffconfig%dims), size(OPP%dirconfig%dims)
      do k=1,size(OPP%dirconfig%dims)
        print *,myid,'dim ',trim(OPP%dirconfig%dims(k)%dimname), OPP%dirconfig%offsets(k), OPP%dirconfig%dims(k)%vrange, ':', OPP%dirconfig%dims(k)%v
      enddo
    endif
end subroutine

  subroutine scatter_LUTtables(OPP, comm)
      use m_optprop_parameters, only: luse_memory_map
      integer(mpiint) ,intent(in) :: comm
      class(t_optprop_LUT) :: OPP

      integer(mpiint) :: myid, ierr
      real(irealLUT), pointer :: mmap_ptr(:,:)

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

        mmap_ptr => NULL()
        if(associated(OPP%Sdir%c)) then
          call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr, OPP%Sdir%c)
          deallocate(OPP%Sdir%c)
        else
          call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', mmap_ptr, ierr)
        endif
        OPP%Sdir%c => mmap_ptr

        mmap_ptr => NULL()
        if(associated(OPP%Tdir%c)) then
          call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr, OPP%Tdir%c)
          deallocate(OPP%Tdir%c)
        else
          call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', mmap_ptr, ierr)
        endif
        OPP%Tdir%c => mmap_ptr

      else
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
      print *,myid,'Dimension '//itoa(i)//' '//trim(config%dims(i)%dimname)//' size '//itoa(config%dims(i)%N), '(', config%dims(i)%vrange, ')'
    enddo
  end subroutine

  subroutine check_if_samplepts_in_LUT_bounds(sample_pts, config)
    real(irealLUT),intent(in) :: sample_pts(:)
    type(t_LUT_config), allocatable, intent(in) :: config
    integer(mpiint) :: ierr, kdim

    ierr = 0
    do kdim = 1,size(sample_pts)
      if(sample_pts(kdim).lt.config%dims(kdim)%vrange(1).or.sample_pts(kdim).gt.config%dims(kdim)%vrange(2)) then
        print *,'ERROR value in dimension '//itoa(kdim)//' is outside of LUT range', sample_pts(kdim), 'not in:', config%dims(kdim)%vrange
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
    real(irealLUT) :: pti(size(sample_pts)), norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%dirconfig%dims)) then
        call print_LUT_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%dirconfig%dims)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%dirconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti, kind=iintegers))
      C = OPP%Tdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Tdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
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
          print *,'SUM dir2dir coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),' :: coeff',C( src:size(C):OPP%dir_streams )
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
    real(irealLUT) :: pti(size(sample_pts)), norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%dirconfig%dims)) then
        call print_LUT_config(OPP%dirconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%dirconfig%dims)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%dirconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti, kind=iintegers))
      C = OPP%Sdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Sdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
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
          print *,'SUM dir2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),' :: coeff',C( src:size(C):OPP%dir_streams )
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
    real(irealLUT) :: pti(size(sample_pts)), norm

    if(ldebug_optprop) then
      if(size(sample_pts).ne.size(OPP%diffconfig%dims)) then
        call print_LUT_config(OPP%diffconfig)
        call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in LUT ' &
          //itoa(size(sample_pts, kind=iintegers))//'/'//itoa(size(OPP%diffconfig%dims, kind=iintegers)))
      endif
      call check_if_samplepts_in_LUT_bounds(sample_pts, OPP%diffconfig)
    endif

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%diffconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%diffconfig%offsets, nint(pti, kind=iintegers))
      C = OPP%Sdiff%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Sdiff%c, OPP%diffconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
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
          print *,'SUM diff2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%diff_streams)),' :: coeff',C(src:size(C):OPP%diff_streams)
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
  end subroutine

end module
