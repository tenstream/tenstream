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

module m_netcdfIO
#if defined(__INTEL_COMPILER)
      use ifport
#endif

  use iso_c_binding
  use iso_fortran_env, only: REAL32, REAL64, INT32, INT64

  use netcdf
  use m_data_parameters, only :   &
      default_str_len, &
      ireals,          &
      iintegers, mpiint
  use m_helper_functions, only : CHKWARN, CHKERR, itoa, ftoa, get_arg
  implicit none

  private
  public :: ncwrite, ncload, acquire_file_lock, release_file_lock, &
    get_global_attribute, set_global_attribute

!  integer :: v=11
  integer,parameter :: deflate_lvl=1
!  real(ireals),parameter :: maxwait=600 !in seconds
!  real(ireals),parameter :: waitinterval=.01 ! amount of cpu time to wait before trying anew in seconds
!  integer :: iwait
!  character(default_str_len+10) :: lockfile
  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.

  interface ncwrite
    module procedure ncwrite_1d_r32, ncwrite_2d_r32, ncwrite_3d_r32, ncwrite_4d_r32, ncwrite_5d_r32, ncwrite_6d_r32, ncwrite_7d_r32, &
                     ncwrite_1d_r64, ncwrite_2d_r64, ncwrite_3d_r64, ncwrite_4d_r64, ncwrite_5d_r64, ncwrite_6d_r64, ncwrite_7d_r64
  end interface
  interface ncload
    module procedure &
        ncload_1d_r32, ncload_2d_r32, ncload_3d_r32, ncload_4d_r32, ncload_5d_r32, ncload_6d_r32, ncload_7d_r32, &
        ncload_1d_r64, ncload_2d_r64, ncload_3d_r64, ncload_4d_r64, ncload_5d_r64, ncload_6d_r64, ncload_7d_r64, &
        ncload_2d_r32_ptr, ncload_2d_r64_ptr, &
        ncload_1dint, ncload_2dint
  end interface
  interface get_global_attribute
    module procedure get_global_attribute_str, &
        get_global_attribute_i32, get_global_attribute_i64, &
        get_global_attribute_r32, get_global_attribute_r64
  end interface
  interface set_global_attribute
    module procedure set_global_attribute_str, &
        set_global_attribute_i32, set_global_attribute_i64, &
        set_global_attribute_r32, set_global_attribute_r64
  end interface

#if __HAVE_NC_SET_LOG_LEVEL__
  interface
    function nf90_set_log_level(level) bind (C, name = "nc_set_log_level")
      use iso_c_binding
      implicit none
      integer(c_int) :: nf90_set_log_level
      integer(c_int), intent (in) :: level
    end function nf90_set_log_level
  end interface
#endif

  contains

#if !(__HAVE_NC_SET_LOG_LEVEL__)
  integer(c_int) function nf90_set_log_level(level) result(r)
    use iso_c_binding
    integer(c_int), intent (in) :: level
    r = level ! prevent unused warning
    r = 0
  end function
#endif

    subroutine ncwrite_1d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_2d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_3d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_4d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_5d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_6d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_7d_r32(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL32),intent(in) :: arr(:,:,:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine

    subroutine ncwrite_1d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_2d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_3d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_4d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_5d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_6d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_7d_r64(groups,arr,ierr,arr_shape,startp,countp,stride,map)
        real(REAL64),intent(in) :: arr(:,:,:,:,:,:,:)
        include 'netcdfio_write.inc'
    end subroutine

    subroutine ncload_1dint(groups,arr,ierr,lverbose)
        integer(iintegers),allocatable,intent(inout) :: arr(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2dint(groups,arr,ierr,lverbose)
        integer(iintegers),allocatable,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine

     subroutine ncload_1d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d_r32_ptr(groups,arr,ierr,lverbose)
        real(REAL32),pointer,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1_ptr.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_3d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_4d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_5d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_6d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_7d_r32(groups,arr,ierr,lverbose)
        real(REAL32),allocatable,intent(inout) :: arr(:,:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6),dimsize(7)))
        include 'netcdfio_read_2.inc'
    end subroutine
     subroutine ncload_1d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d_r64_ptr(groups,arr,ierr,lverbose)
        real(REAL64),pointer,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1_ptr.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_3d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_4d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_5d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_6d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_7d_r64(groups,arr,ierr,lverbose)
        real(REAL64),allocatable,intent(inout) :: arr(:,:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6),dimsize(7)))
        include 'netcdfio_read_2.inc'
    end subroutine

    subroutine cpusleep(sec)
        real(ireals) :: sec
        real(ireals) :: t1,t2,dither
        call random_number(dither)
        dither = ( -1.0_ireals + 2*dither ) *.1 + sec !10% dithering on time
        call cpu_time(t1)
        do
          call cpu_time(t2)
          if(t2-t1.gt.dither) return
        enddo
    end subroutine

    subroutine nccheck(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
          if(ldebug) print *,"NetCDF Error ::",trim(nf90_strerror(status)),'::',status
        end if
    end subroutine nccheck

    function get_pid_macro()
    integer(iintegers) :: get_pid_macro
#ifdef _XLF
        get_pid_macro=-1
!        call MGPID(get_pid_macro) 
#else
        get_pid_macro=getpid()
#endif
    end function

    subroutine acquire_file_lock(fname, flock_unit, ierr, lock_fname, blocking, waittime, waitinterval)
      character(len=*), intent(in) :: fname
      integer(mpiint), intent(out) :: flock_unit, ierr

      character(len=*), intent(in), optional :: lock_fname
      logical, intent(in), optional :: blocking          ! if the call blocks until we get the lock or return immediately
      integer(mpiint), intent(in), optional :: waittime  ! if blocking, wait fopr waittime seconds before exiting with error
      real(ireals), intent(in), optional :: waitinterval ! amount of cpu time to wait before trying anew in seconds

      logical :: lblocking
      character(len=default_str_len+5) :: lockfile

      integer(mpiint) :: iwait, maxwait
      real(ireals) :: winterval

      lockfile = trim(get_arg(trim(fname)//'.lock', lock_fname))
      lblocking = get_arg(.True., blocking)
      maxwait = get_arg(120, waittime)
      winterval = get_arg(.5_ireals, waitinterval)

      do iwait=1,int(maxwait/winterval)
        open(newunit=flock_unit,file=lockfile,status='new',err=99)
        write(flock_unit,*) 'file is locked by process: ',get_pid_macro()
        ierr = 0
        return

        99 continue
        if(lblocking) then
          call cpusleep(winterval)
        endif
      enddo
      ierr = iwait
      if(lblocking) then
        call CHKERR(1_mpiint, 'Couldnt lock file '//fname//&
          ' .. waited now for quite a while but we couldnt open the lock: '//lockfile)
      endif
    end subroutine

    subroutine release_file_lock(flock_unit, ierr)
      integer, intent(inout) :: flock_unit
      integer, intent(out) :: ierr
      integer :: i, ios
      logical :: lexist, lnamed, lopened
      real(ireals),parameter :: waitinterval=.1 ! amount of cpu time to wait before trying anew in seconds
      inquire(unit=flock_unit, exist=lexist)
      if(.not. lexist) then
        ierr=2
        return
      else
        call cpusleep(waitinterval)
      endif
      inquire(unit=flock_unit, iostat=ios)
      call CHKERR(ios, 'IOSTAT not 0... is =>'//itoa(ios))

      inquire(unit=flock_unit, named=lnamed)
      if(.not.lnamed) call CHKERR(4_mpiint, 'Release lock file not named')

      inquire(unit=flock_unit, opened=lopened)
      if(.not.lopened) call CHKERR(4_mpiint, 'Release lock file not opened')

      do i=1,10
        close(unit=flock_unit,status='delete',err=99)
        ierr = 0
        return
        99 continue
        call cpusleep(waitinterval)
      enddo
      ierr = 1
      call CHKERR(1_mpiint, 'Error releasing file lock for unit '//itoa(flock_unit))
    end subroutine

    subroutine get_global_attribute_str(fname, attr_name, attr)
      character(len=*) :: fname, attr_name, attr
      integer :: ncid, ierr
      integer :: attrLength
      ierr = nf90_open(trim(fname), nf90_nowrite, ncid); call nccheck(ierr)
      ierr = nf90_inquire_attribute(ncid, nf90_global, attr_name, len=attrLength); call nccheck(ierr)
      if(attrLength.gt.default_str_len) then
        call CHKERR(1_mpiint, 'Cant fit attribute name into default_str_len, its too long: '//itoa(attrLength))
      endif
      ierr = nf90_get_att(ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    end subroutine
    subroutine get_global_attribute_i32(fname, attr_name, attr)
      character(len=*) :: fname, attr_name
      integer(INT32) :: attr
      integer :: ncid, ierr
      ierr = nf90_open(trim(fname), nf90_nowrite, ncid); call nccheck(ierr)
      ierr = nf90_get_att(ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    end subroutine
    subroutine get_global_attribute_i64(fname, attr_name, attr)
      character(len=*) :: fname, attr_name
      integer(INT64) :: attr
      integer :: ncid, ierr
      ierr = nf90_open(trim(fname), nf90_nowrite, ncid); call nccheck(ierr)
      ierr = nf90_get_att(ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    end subroutine
    subroutine get_global_attribute_r32(fname, attr_name, attr)
      character(len=*) :: fname, attr_name
      real(REAL32) :: attr
      integer :: ncid, ierr
      ierr = nf90_open(trim(fname), nf90_nowrite, ncid); call nccheck(ierr)
      ierr = nf90_get_att(ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    end subroutine
    subroutine get_global_attribute_r64(fname, attr_name, attr)
      character(len=*) :: fname, attr_name
      real(REAL64) :: attr
      integer :: ncid, ierr
      ierr = nf90_open(trim(fname), nf90_nowrite, ncid); call nccheck(ierr)
      ierr = nf90_get_att(ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    end subroutine

    subroutine set_global_attribute_str(fname, attr_name, attr)
      character(len=*), intent(in) :: attr
      include 'netcdfio_attr.inc'
    end subroutine
    subroutine set_global_attribute_i32(fname, attr_name, attr)
      integer(INT32), intent(in) :: attr
      include 'netcdfio_attr.inc'
    end subroutine
    subroutine set_global_attribute_i64(fname, attr_name, attr)
      integer(INT64), intent(in) :: attr
      include 'netcdfio_attr.inc'
    end subroutine
    subroutine set_global_attribute_r32(fname, attr_name, attr)
      real(REAL32), intent(in) :: attr
      include 'netcdfio_attr.inc'
    end subroutine
    subroutine set_global_attribute_r64(fname, attr_name, attr)
      real(REAL64), intent(in) :: attr
      include 'netcdfio_attr.inc'
    end subroutine

  end module

  !program main
  !      use arrayIO
  !      implicit none
  !
  !      integer :: r(10,2,6)
  !      integer,allocatable :: o(:,:,:)
  !
  !      character(100) :: fname='test3d'
  !
  !
  !      integer i,j
  !
  !      do i=1,2
  !        do j=1,10
  !          r(j,i,:) = 10*(i-1)+j
  !        enddo
  !      enddo
  !
  !      call write_bin(fname,r)
  !      r = 0
  !      call read_bin(fname,o)
  !
  !      print *,o
  !
  !end program


