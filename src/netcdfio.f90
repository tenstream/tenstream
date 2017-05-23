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
  USE m_data_parameters, ONLY :   &
      default_str_len, &
      ireals,          &
      iintegers
  implicit none

  private
  public :: ncwrite,ncload

  integer :: u=10,v=11
  integer,parameter :: deflate_lvl=9
  real(ireals),parameter :: maxwait=600 !in seconds
  real(ireals),parameter :: waitinterval=.01 ! amount of cpu time to wait before trying anew in seconds
  integer :: iwait
  character(default_str_len+10) :: lockfile
  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.

  interface ncwrite
    module procedure ncwrite_1d,ncwrite_2d,ncwrite_3d,ncwrite_4d,ncwrite_5d,ncwrite_7d
  end interface
  interface ncload
    module procedure ncload_1d,ncload_2d,ncload_3d,ncload_4d,ncload_5d,ncload_7d, ncload_1dint, ncload_2dint
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
        integer(iintegers),allocatable,intent(out) :: arr(:)
        integer(4),allocatable                     :: tmp(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        if(ierr.eq.0) allocate(tmp(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2dint(groups,arr,ierr)
        integer(iintegers),allocatable,intent(out) :: arr(:,:)
        integer(4),allocatable                     :: tmp(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine

     subroutine ncload_1d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:)
        real(4)     ,allocatable             :: tmp(:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1)))
        if(ierr.eq.0) allocate(tmp(dimsize(1)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_2d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:)
        real(4)     ,allocatable             :: tmp(:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_3d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:,:)
        real(4)     ,allocatable             :: tmp(:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2),dimsize(3)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_4d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:,:,:)
        real(4)     ,allocatable             :: tmp(:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_5d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:,:,:,:)
        real(4)     ,allocatable             :: tmp(:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_6d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:,:,:,:,:)
        real(4)     ,allocatable             :: tmp(:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6)))
        include 'netcdfio_read_2.inc'
    end subroutine
    subroutine ncload_7d(groups,arr,ierr)
        real(ireals),allocatable,intent(out) :: arr(:,:,:,:,:,:,:)
        real(4)     ,allocatable             :: tmp(:,:,:,:,:,:,:)
        include 'netcdfio_read_1.inc'
        if(ierr.eq.0) allocate(arr(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6),dimsize(7)))
        if(ierr.eq.0) allocate(tmp(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5),dimsize(6),dimsize(7)))
        include 'netcdfio_read_2.inc'
    end subroutine

    subroutine cpusleep(sec)
        real(ireals) :: sec
        real(ireals) :: t1,t2,dither
        call random_number(dither)
        dither = ( -1.0_ireals + 2*dither ) *.1*sec !10% dithering on time
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


