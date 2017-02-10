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

module m_helper_functions_dp
      use m_data_parameters,only : iintegers,ireal_dp,imp_real_dp,imp_int,imp_logical,mpiint
      use m_helper_functions, only: CHKERR
      use mpi

      implicit none

      private
      public imp_bcast,norm,deg2rad,rmse,mean,approx,rel_approx,delta_scale_optprop,delta_scale,cumsum,inc, &
          mpi_logical_and,mpi_logical_or,imp_allreduce_min,imp_allreduce_max,imp_reduce_sum

      interface imp_bcast
        module procedure imp_bcast_real_1d,imp_bcast_real_2d,imp_bcast_real_3d,imp_bcast_real_5d,imp_bcast_int_1d,imp_bcast_int_2d,imp_bcast_int,imp_bcast_real,imp_bcast_logical
      end interface

      integer(mpiint) :: mpierr
      real(ireal_dp),parameter :: zero=0, one=1, pi=3.141592653589793

    contains
      pure elemental subroutine inc(x,i)
          real(ireal_dp),intent(inout) :: x
          real(ireal_dp),intent(in) :: i
          x=x+i
      end subroutine

      pure function norm(v)
          real(ireal_dp) :: norm
          real(ireal_dp),intent(in) :: v(:)
          norm = sqrt(dot_product(v,v))
      end function

      elemental function deg2rad(deg)
          real(ireal_dp) :: deg2rad
          real(ireal_dp),intent(in) :: deg
          deg2rad = deg *pi/180._ireal_dp
      end function

      pure function rmse(a,b)
          real(ireal_dp) :: rmse(2)
          real(ireal_dp),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireal_dp) :: mean
          real(ireal_dp),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor
          if(present(precision) ) then
            factor = precision
          else
            factor = 10._ireal_dp*epsilon(b)
          endif
          if( a.le.b+factor .and. a.ge.b-factor ) then
            approx = .True.
          else
            approx = .False.
          endif
      end function
      elemental logical function rel_approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor,rel_error
          if(present(precision) ) then
            factor = precision
          else
            factor = 10*epsilon(b)
          endif
          rel_error = abs( (a-b)/ max(epsilon(a), ( (a+b)*.5_ireal_dp ) ) )

          if( rel_error .lt. precision ) then
            rel_approx = .True.
          else
            rel_approx = .False.
          endif
      end function


      function mpi_logical_and(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_and
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_and, 1_mpiint, imp_logical, MPI_LAND, comm, mpierr); call CHKERR(mpierr)
      end function
      function mpi_logical_or(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_or
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_or, 1_mpiint, imp_logical, MPI_LOR, comm, mpierr); call CHKERR(mpierr)
      end function

      subroutine imp_allreduce_min(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MIN,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_allreduce_max(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MAX,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_reduce_sum(comm, v, myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: v
          integer(mpiint),intent(in) :: myid
          integer(mpiint) :: commsize

          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return 

          if(myid.eq.0) then
            call mpi_reduce(MPI_IN_PLACE, v, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          else
            call mpi_reduce(v, MPI_IN_PLACE, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          endif
      end subroutine

      subroutine  imp_bcast_logical(comm, val, sendid, myid)
          integer(mpiint),intent(in) :: comm
          logical,intent(inout) :: val
          integer(mpiint),intent(in) :: sendid,myid

          call mpi_bcast(val, 1_mpiint, imp_logical, sendid, comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int(comm, val,sendid,myid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid,myid

          call mpi_bcast(val,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_1d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_2d(comm, arr,sendid,myid)!
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(2)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real(comm, val,sendid,myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid,myid

          call mpi_bcast(val,1_mpiint,imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_1d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot

          if(sendid.eq.myid) Ntot = size(arr)
          call imp_bcast_int(comm, Ntot, sendid, myid)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_2d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(2)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_3d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(3)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,3_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_5d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:,:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(5)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,5_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3), Ntot(4), Ntot(5) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine

      elemental subroutine delta_scale( kabs,ksca,g,factor ) 
          real(ireal_dp),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: dtau, w0
          dtau = max( kabs+ksca, epsilon(dtau) )
          w0   = ksca/dtau
          g    = g

          if(present(factor)) then
            call delta_scale_optprop( dtau, w0, g, factor)
          else
            call delta_scale_optprop( dtau, w0, g)
          endif

          kabs= dtau * (one-w0)
          ksca= dtau * w0
      end subroutine
      elemental subroutine delta_scale_optprop( dtau, w0, g,factor) 
          real(ireal_dp),intent(inout) :: dtau,w0,g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: f

          g = min( g, one-epsilon(g)*10)
          if(present(factor)) then
            f = factor
          else
            f = g**2
          endif
          dtau = dtau * ( one - w0 * f )
          g    = ( g - f ) / ( one - f )
          w0   = w0 * ( one - f ) / ( one - f * w0 )
      end subroutine

      function cumsum(arr)
          real(ireal_dp),intent(in) :: arr(:)
          real(ireal_dp) :: cumsum(size(arr))
          integer :: i
          cumsum(1) = arr(1)
          do i=2,size(arr)
            cumsum(i) = cumsum(i-1) + arr(i)
          enddo
      end function


      end module
