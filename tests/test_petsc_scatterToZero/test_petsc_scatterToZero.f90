module test_petsc_scatterToZero

#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : iintegers, ireals, mpiint, zero, one
  use m_pprts, only : init_pprts, t_solver_8_10, destroy_pprts
  use m_petsc_helpers, only : petscVecToF90, petscGlobalVecToZero, f90VecToPetsc
  use m_helper_functions, only : CHKERR

  use pfunit_mod

  implicit none

contains

  @test(npes =[2])
  subroutine petsc_scatterToZero(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, ierr

    integer(iintegers),parameter :: nxp=9,nyp=9,nv=2
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: dz=dx
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:,:) :: local_arr
    real(ireals),allocatable,dimension(:,:,:,:) :: local_arr_2
    real(ireals),allocatable,dimension(:,:,:,:) :: global_arr_on_rank0

    integer(iintegers) :: i,j,k,d

    type(t_solver_8_10) :: solver

    type(tVec) :: gvec, lvec

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    allocate(local_arr(solver%C_one%dof, solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    local_arr = myid

    ! Gen new global vector
    call DMGetGlobalVector(solver%C_one%da, gvec, ierr); call CHKERR(ierr)

    ! copy local fortran array into global vec
    print *,myid,'copy local fortran array into global vec'
    call f90VecToPetsc(local_arr, solver%C_one%da, gvec)

    ! Check that we can get an local array from a global petsc vec
    call petscVecToF90(gvec, solver%C_one%da, local_arr_2)
    call assert_equivalence(local_arr, local_arr_2)

    deallocate(local_arr_2)

    ! do that feat again with a 3D array
    call f90VecToPetsc(local_arr(1,:,:,:), solver%C_one%da, gvec)
    call petscVecToF90(gvec, solver%C_one%da, local_arr_2)
    call assert_equivalence(local_arr, local_arr_2)

    ! copy global vec to rank 0, lVec(full size)
    call petscGlobalVecToZero(gVec, solver%C_one%da, lVec)

    ! Put data from petsc vec into fortran array
    if(myid.eq.0) then
      call petscVecToF90(lvec, solver%C_one%da, global_arr_on_rank0, opt_l_only_on_rank0=.True.)
      print *,myid,' size(global_arr_on_rank0)', size(global_arr_on_rank0)
      do j=1,ubound(global_arr_on_rank0,4)
        do i=1,ubound(global_arr_on_rank0,3)
          do k=1,ubound(global_arr_on_rank0,2)
            do d=1,ubound(global_arr_on_rank0,1)
              if(d.le.solver%C_one%dof .and. k.le.solver%C_one%zm .and. i.le.solver%C_one%xm .and. j.le.solver%C_one%ym) then
                @assertEqual(zero, global_arr_on_rank0(d,k,i,j))
              else
                @assertEqual(one, global_arr_on_rank0(d,k,i,j))
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    ! Debug Output
    call PetscObjectSetName(gvec, 'VecGlobal', ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_gvec', ierr); call CHKERR(ierr)
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call PetscObjectSetName(lvec, 'VecLocal', ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(lvec, PETSC_NULL_VEC, '-show_lvec', ierr); call CHKERR(ierr)
    ! Debug Output

    call DMRestoreGlobalVector(solver%C_one%da, gvec, ierr); call CHKERR(ierr)
    call VecDestroy(lvec, ierr); call CHKERR(ierr)
    call destroy_pprts(solver, .True.)
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
