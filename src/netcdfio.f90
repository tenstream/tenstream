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

!  use mpi
  use netcdf
  use m_data_parameters, only :   &
      default_str_len, &
      ireals,          &
      iintegers, mpiint
  use m_helper_functions, only : CHKERR, itoa, get_arg
  implicit none

  private
  public :: ncwrite, ncload, acquire_file_lock, release_file_lock

!  integer :: v=11
  integer,parameter :: deflate_lvl=9
!  real(ireals),parameter :: maxwait=600 !in seconds
!  real(ireals),parameter :: waitinterval=.01 ! amount of cpu time to wait before trying anew in seconds
!  integer :: iwait
!  character(default_str_len+10) :: lockfile
  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.

  interface ncwrite
    module procedure ncwrite_1d, ncwrite_2d, ncwrite_3d, ncwrite_4d, ncwrite_5d, ncwrite_6d, ncwrite_7d
  end interface
  interface ncload
    module procedure ncload_1d, ncload_2d, ncload_3d, ncload_4d, ncload_5d, ncload_6d, ncload_7d, &
        ncload_2d_ptr, &
        ncload_1dint, ncload_2dint
  end interface

  contains

    subroutine ncwrite_1d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:)
        real(4)     :: tmp(size(arr))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_2d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_3d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2),size(arr,3))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_4d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:,:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2),size(arr,3),size(arr,4))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_5d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:,:,:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2),size(arr,3),size(arr,4),size(arr,5))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_6d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:,:,:,:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2),size(arr,3),size(arr,4),size(arr,5),size(arr,6))
        include 'netcdfio_write.inc'
    end subroutine
    subroutine ncwrite_7d(groups,arr,ierr)
        real(ireals),intent(in) :: arr(:,:,:,:,:,:,:)
        real(4)     :: tmp(size(arr,1),size(arr,2),size(arr,3),size(arr,4),size(arr,5),size(arr,6),size(arr,7))
        include 'netcdfio_write.inc'
    end subroutine

    subroutine ncload_1dint(groups,arr,ierr)
        integer(iintegers),allocatable,intent(inout) :: arr(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2dint(groups,arr,ierr)
        integer(iintegers),allocatable,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine

     subroutine ncload_1d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d_ptr(groups,arr,ierr)
        real(ireals),pointer,intent(inout) :: arr(:,:)
        include 'netcdfio_read_1_ptr.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_3d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_4d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_5d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_6d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_7d(groups,arr,ierr)
        real(ireals),allocatable,intent(inout) :: arr(:,:,:,:,:,:,:)
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

    subroutine acquire_file_lock(fname, flock_unit, ierr, lock_fname, blocking, waittime)
      character(len=*), intent(in) :: fname
      integer(mpiint), intent(out) :: flock_unit, ierr

      character(len=*), intent(in), optional :: lock_fname
      logical, intent(in), optional :: blocking
      integer(mpiint), intent(in), optional :: waittime

      logical :: lblocking
      character(len=default_str_len+5) :: lockfile

      real(ireals),parameter :: waitinterval=.1 ! amount of cpu time to wait before trying anew in seconds
      integer(mpiint) :: iwait, maxwait

      lockfile = trim(get_arg(trim(fname)//'.lock', lock_fname))
      lblocking = get_arg(.True., blocking)
      maxwait = get_arg(120, waittime)

      do iwait=1,int(maxwait/waitinterval)
        open(newunit=flock_unit,file=lockfile,status='new',err=99)
        write(flock_unit,*) 'file is locked by process: ',get_pid_macro()
        ierr = 0
        return

        99 continue
        if(lblocking) then
          call cpusleep(waitinterval)
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


