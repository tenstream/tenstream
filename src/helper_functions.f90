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

module m_helper_functions
      use m_data_parameters,only : iintegers,ireals,pi,zero,one,imp_real,imp_int,imp_logical,mpiint,imp_comm

      use mpi

      implicit none

      private
      public imp_bcast,norm,rad2deg,deg2rad,rmse,mean,approx,rel_approx,delta_scale_optprop,delta_scale,cumsum,inc, &
          mpi_logical_and,mpi_logical_or,imp_allreduce_min,imp_allreduce_max,imp_reduce_sum, search_sorted_bisection, &
          gradient, read_ascii_file_2d, meanvec, swap, imp_allgather_int_inplace, reorder_mpi_comm, CHKERR

      interface imp_bcast
        module procedure imp_bcast_real_1d,imp_bcast_real_2d,imp_bcast_real_3d,imp_bcast_real_5d,imp_bcast_int_1d,imp_bcast_int_2d,imp_bcast_int,imp_bcast_real,imp_bcast_logical
      end interface

      integer(mpiint) :: mpierr

    contains
      pure elemental subroutine swap(x,y)
          real(ireals),intent(inout) :: x,y
          real(ireals) :: tmp
          tmp = x
          x = y
          y = tmp
      end subroutine
      pure elemental subroutine inc(x,i)
          real(ireals),intent(inout) :: x
          real(ireals),intent(in) :: i
          x=x+i
      end subroutine

      subroutine CHKERR(ierr)
        integer(mpiint),intent(in) :: ierr
        if(ierr.ne.0) call mpi_abort(imp_comm, ierr, mpierr)
      end subroutine

      pure function gradient(v)
          real(ireals),intent(in) :: v(:)
          real(ireals) :: gradient(size(v)-1)
          gradient = v(2:size(v))-v(1:size(v)-1)
      end function

      pure function meanvec(v)
          real(ireals),intent(in) :: v(:)
          real(ireals) :: meanvec(size(v)-1)
          meanvec = (v(2:size(v))+v(1:size(v)-1))*.5_ireals
      end function

      pure function norm(v)
          real(ireals) :: norm
          real(ireals),intent(in) :: v(:)
          norm = sqrt(dot_product(v,v))
      end function

      elemental function deg2rad(deg)
          real(ireals) :: deg2rad
          real(ireals),intent(in) :: deg
          deg2rad = deg *pi/180._ireals
      end function
      elemental function rad2deg(rad)
          real(ireals) :: rad2deg
          real(ireals),intent(in) :: rad
          rad2deg = rad /pi*180._ireals
      end function

      pure function rmse(a,b)
          real(ireals) :: rmse(2)
          real(ireals),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireals) :: mean
          real(ireals),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
          real(ireals),intent(in) :: a,b
          real(ireals),intent(in),optional :: precision
          real(ireals) :: factor
          if(present(precision) ) then
            factor = precision
          else
            factor = 10._ireals*epsilon(b)
          endif
          if( a.le.b+factor .and. a.ge.b-factor ) then
            approx = .True.
          else
            approx = .False.
          endif
      end function
      elemental logical function rel_approx(a,b,precision)
          real(ireals),intent(in) :: a,b
          real(ireals),intent(in),optional :: precision
          real(ireals) :: factor,rel_error
          if(present(precision) ) then
            factor = precision
          else
            factor = 10*epsilon(b)
          endif
          rel_error = abs( (a-b)/ max(epsilon(a), ( (a+b)*.5_ireals ) ) )

          if( rel_error .lt. precision ) then
            rel_approx = .True.
          else
            rel_approx = .False.
          endif
      end function


      function mpi_logical_and(comm,lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_and
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_and, 1_mpiint, imp_logical, MPI_LAND, comm, mpierr); call CHKERR(mpierr)
      end function
      function mpi_logical_or(comm,lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_or
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_or, 1_mpiint, imp_logical, MPI_LOR, comm, mpierr); call CHKERR(mpierr)
      end function

      subroutine imp_allreduce_min(comm,v,r)
          integer(mpiint),intent(in) :: comm
          real(ireals),intent(in) :: v
          real(ireals),intent(out) :: r
          call mpi_allreduce(v,r,1_mpiint,imp_real, MPI_MIN,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_allreduce_max(comm,v,r)
          integer(mpiint),intent(in) :: comm
          real(ireals),intent(in) :: v
          real(ireals),intent(out) :: r
          call mpi_allreduce(v,r,1_mpiint,imp_real, MPI_MAX,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_reduce_sum(comm,v)
          real(ireals),intent(inout) :: v
          integer(mpiint),intent(in) :: comm
          integer(mpiint) :: commsize, myid
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(myid.eq.0) then
            call mpi_reduce(MPI_IN_PLACE, v, 1_mpiint, imp_real, MPI_SUM, 0_mpiint, comm, mpierr); call CHKERR(mpierr)
          else
            call mpi_reduce(v, MPI_IN_PLACE, 1_mpiint, imp_real, MPI_SUM, 0_mpiint, comm, mpierr); call CHKERR(mpierr)
          endif
      end subroutine

      subroutine imp_allgather_int_inplace(comm,v)
        integer(mpiint),intent(in) :: comm
        integer(iintegers),intent(inout) :: v(:)
        call mpi_allgather(MPI_IN_PLACE, 0_mpiint, MPI_DATATYPE_NULL, v, 1_mpiint, imp_int, comm, mpierr); call CHKERR(mpierr)
      end subroutine

      subroutine  imp_bcast_logical(comm,val,sendid)
          integer(mpiint),intent(in) :: comm
          logical,intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val, 1_mpiint, imp_logical, sendid, comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int(comm,val,sendid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_1d(comm,arr,sendid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_2d(comm,arr,sendid)!
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(2)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real(comm,val,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireals),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_real,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_1d(comm,arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireals),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_2d(comm,arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireals),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(2)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_3d(comm,arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireals),allocatable,intent(inout) :: arr(:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(3)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,3_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3) ) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_5d(comm,arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireals),allocatable,intent(inout) :: arr(:,:,:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(5)
          integer(mpiint) :: commsize
          call MPI_Comm_size(comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,5_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3), Ntot(4), Ntot(5) ) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine

      elemental subroutine delta_scale( kabs,ksca,g,factor ) 
          real(ireals),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireals),intent(in),optional :: factor
          real(ireals) :: dtau, w0
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
          real(ireals),intent(inout) :: dtau,w0,g
          real(ireals),intent(in),optional :: factor
          real(ireals) :: f

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
          real(ireals),intent(in) :: arr(:)
          real(ireals) :: cumsum(size(arr))
          integer :: i
          cumsum(1) = arr(1)
          do i=2,size(arr)
            cumsum(i) = cumsum(i-1) + arr(i)
          enddo
      end function

      subroutine read_ascii_file_2d(filename, arr, ncolumns, skiplines, ierr)
          character(len=*),intent(in) :: filename
          integer(iintegers),intent(in) :: ncolumns
          integer(iintegers),intent(in),optional :: skiplines

          real(ireals),allocatable,intent(out) :: arr(:,:)

          integer(mpiint) :: ierr

          real :: line(ncolumns)

          integer(iintegers) :: unit, nlines, i, io
          logical :: file_exists=.False.

          ierr=0
          inquire(file=filename, exist=file_exists)

          if(.not.file_exists) then
            print *,'File ',trim(filename), 'does not exist!'
            ierr=1
            return
          endif

          open(newunit=unit, file=filename)
          if(present(skiplines)) then
              do i=1,skiplines
                  read(unit,*)
              enddo
          endif

          nlines = 0
          do
              read(unit, *, iostat=io) line
              !print *,'line',line
              if (io/=0) exit
              nlines = nlines + 1
          end do

          rewind(unit)
          if(present(skiplines)) then
              do i=1,skiplines
                  read(unit,*)
              enddo
          endif

          allocate(arr(nlines,ncolumns))

          do i=1,nlines
              read(unit, *, iostat=io) line
              arr(i,:) = line
          end do

          close(unit)
          print *,'I read ',nlines,'lines'
      end subroutine

      function search_sorted_bisection(arr,val) ! return index+residula i where val is between arr(i) and arr(i+1)
        real(ireals) :: search_sorted_bisection
        real(ireals),intent(in) :: arr(:)
        real(ireals),intent(in) :: val
        real(ireals) :: loc_increment
        integer(iintegers) :: i,j,k

        i=lbound(arr,1)
        j=ubound(arr,1)

        if(arr(i).le.arr(j)) then ! ascending order
          do
            k=(i+j)/2
            if (val < arr(k)) then
              j=k
            else
              i=k
            endif
            if (i+1 >= j) then ! only single or tuple left
              ! i is left bound and j is right bound index
              if(i.eq.j) then
                loc_increment = zero
              else
                loc_increment = (val - arr(i)) / ( arr(j) - arr(i) )
              endif
              search_sorted_bisection= min(max(one*lbound(arr,1), i + loc_increment), one*ubound(arr,1)) ! return `real-numbered` location of val
              return
            endif
          end do
        else !descending order

          do
            k=(i+j)/2
            if (val > arr(k)) then
              j=k
            else
              i=k
            endif
            if (i+1 >= j) then ! only single or tuple left
              ! i is left bound and j is right bound index
              if(i.eq.j) then
                loc_increment = zero
              else
                loc_increment = (val - arr(j)) / ( arr(i) - arr(j) )
              endif
              search_sorted_bisection= min(max(one*lbound(arr,1), j - loc_increment), one*ubound(arr,1)) ! return `real-numbered` location of val
              return
            endif
          end do
        endif
      end function

      subroutine reorder_mpi_comm(icomm, Nrank_x, Nrank_y, new_comm)
        integer(mpiint), intent(in) :: icomm
        integer(mpiint), intent(out) :: new_comm
        integer(iintegers) :: Nrank_x, Nrank_y

        ! This is the code snippet from Petsc FAQ to change from PETSC (C) domain splitting to MPI(Fortran) domain splitting
        ! the numbers of processors per direction are (int) x_procs, y_procs, z_procs respectively
        ! (no parallelization in direction 'dir' means dir_procs = 1)

        integer(iintegers) :: x,y
        integer(mpiint) :: orig_id, petsc_id ! id according to fortran decomposition

        call MPI_COMM_RANK( icomm, orig_id, mpierr ); call CHKERR(mpierr)

        ! calculate coordinates of cpus in MPI ordering:
        x = int(orig_id) / Nrank_y
        y = modulo(orig_id ,Nrank_y)

        ! set new rank according to PETSc ordering:
        petsc_id = y*Nrank_x + x

        ! create communicator with new ranks according to PETSc ordering:
        call MPI_Comm_split(icomm, 1_iintegers, petsc_id, new_comm, mpierr)

        !print *,'Reordering communicator'
        !print *,'setup_petsc_comm: MPI_COMM_WORLD',orig_id,'calc_id',petsc_id
      end subroutine
  end module
