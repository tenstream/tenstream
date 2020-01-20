module test_convolution

#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, i1, i2
  use m_tenstream_options, only : read_commandline_options
  use m_pprts_base, only : t_solver_3_10
  use m_pprts, only : init_pprts, destroy_pprts
  use m_petsc_helpers, only : petscVecToF90, petscGlobalVecToZero, f90VecToPetsc, dmda_convolve_ediff_srfc
  use m_helper_functions, only : CHKERR

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver
  type(tVec) :: gvec, lvec

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
    call read_commandline_options(this%getMpiCommunicator())
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call destroy_pprts(solver, lfinalizepetsc=.True.)
  end subroutine teardown

  @test(npes =[2])
  subroutine dmda_convolve_2(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, ierr

    integer(iintegers),parameter :: nxp=8,nyp=8,nv=2
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: dz=dx
    real(ireals) :: dz1d(nv)
    real(ireals), parameter :: eps=sqrt(epsilon(eps))

    real(ireals),allocatable,dimension(:,:,:,:) :: local_arr
    real(ireals) :: targt

    integer(iintegers) :: j,k

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    associate(C=>solver%C_diff)
      allocate(local_arr(C%dof, C%zm , C%xm,  C%ym ))
      local_arr = modulo(myid,2)

      call dmda_convolve_ediff_srfc(C%da, i1, local_arr(:,C%zs+1, :, :))
      do k=0,numnodes-1
        if(myid.eq.k) then
          do j=1,C%ym
            print *, myid, j, local_arr(1,C%zs+1,:,j)
          enddo
        endif
        call mpi_barrier(comm, ierr)
      enddo

      if(modulo(myid,2).eq.0) then
        targt = 1._ireals/3._ireals
      else
        targt = 2._ireals/3._ireals
      endif
      @mpiassertEqual(targt, local_arr(1,C%zs+1,1,1), eps)

      local_arr = modulo(myid,2)
      call dmda_convolve_ediff_srfc(C%da, i2, local_arr(:,C%zm, :, :))
      do k=0,numnodes-1
        if(myid.eq.k) then
          do j=1,C%ym
            print *, myid, j, local_arr(1,C%zm,:,j)
          enddo
        endif
        call mpi_barrier(comm, ierr)
      enddo
      if(modulo(myid,2).eq.0) then
        targt = 2._ireals/5._ireals
      else
        targt = 3._ireals/5._ireals
      endif
      @mpiassertEqual(targt, local_arr(1,C%zm,1,1), eps)
      if(modulo(myid,2).eq.0) then
        targt = 1._ireals/5._ireals
      else
        targt = 4._ireals/5._ireals
      endif
      @mpiassertEqual(targt, local_arr(1,C%zm,1,2), eps)

      local_arr = modulo(myid,2)
      call dmda_convolve_ediff_srfc(C%da, i1, local_arr(:,C%zm, :, :))
      call dmda_convolve_ediff_srfc(C%da, i1, local_arr(:,C%zm, :, :))
      do k=0,numnodes-1
        if(myid.eq.k) then
          do j=1,C%ym
            print *, myid, j, local_arr(1,C%zm,:,j)
          enddo
        endif
        call mpi_barrier(comm, ierr)
      enddo

      if(modulo(myid,2).eq.0) then
        targt = 1._ireals/3._ireals
      else
        targt = 2._ireals/3._ireals
      endif
      @mpiassertEqual(targt, local_arr(1,C%zm,1,1), eps)
      if(modulo(myid,2).eq.0) then
        targt = 1._ireals/9._ireals
      else
        targt = 8._ireals/9._ireals
      endif
      @mpiassertEqual(targt, local_arr(1,C%zm,1,2), eps)
    end associate
  end subroutine

  @test(npes =[1])
  subroutine dmda_convolve_1(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, ierr

    integer(iintegers),parameter :: nxp=5,nyp=5,nv=2
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: dz=dx
    real(ireals) :: dz1d(nv)
    real(ireals), parameter :: eps=sqrt(epsilon(eps))

    real(ireals),allocatable,dimension(:,:,:,:) :: local_arr

    integer(iintegers) :: j,k

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    associate(C=>solver%C_diff)
      allocate(local_arr(C%dof, C%zm , C%xm,  C%ym ))
      local_arr = zero
      local_arr(i1,:,3,3) = one
      local_arr(i2,:,3,3) = one*10

      call dmda_convolve_ediff_srfc(C%da, i1, local_arr(:,C%zs+1, :, :))
      do k=0,numnodes-1
        if(myid.eq.k) then
          do j=1,C%ym
            print *, myid, j, local_arr(1,C%zs+1,:,j)
          enddo
        endif
        call mpi_barrier(comm, ierr)
      enddo
      @assertEqual(one/9._ireals, local_arr(1,C%zs+1,3,3), eps)
      @assertEqual(10*one/9._ireals, local_arr(2,C%zs+1,3,3), eps)


      call dmda_convolve_ediff_srfc(C%da, i2, local_arr(:,C%zm, :, :))
      do k=0,numnodes-1
        if(myid.eq.k) then
          do j=1,C%ym
            print *, myid, j, local_arr(1,C%zm,:,j)
          enddo
        endif
        call mpi_barrier(comm, ierr)
      enddo
      @assertEqual(one/25._ireals, local_arr(1,C%zm,3,3), eps)
      @assertEqual(10*one/25._ireals, local_arr(2,C%zm,3,3), eps)

    end associate
  end subroutine

  subroutine assert_equivalence(a,b)
    real(ireals),dimension(:,:,:,:) :: a,b
    integer(iintegers) :: i,j,k,d
    do j=1,ubound(a,4)
      do i=1,ubound(a,3)
        do k=1,ubound(a,2)
          do d=1,ubound(a,1)
            @assertEqual(a(d,k,i,j),b(d,k,i,j))
          enddo
        enddo
      enddo
    enddo
  end subroutine

end module
