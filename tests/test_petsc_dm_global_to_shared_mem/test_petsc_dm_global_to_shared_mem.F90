module test_petsc_dm_global_to_shared_mem

#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : iintegers, ireals, mpiint, zero, one, &
    & init_mpi_data_parameters
  use m_pprts_base, only : t_solver_3_10, destroy_pprts
  use m_pprts, only : init_pprts
  use m_petsc_helpers, only : &
    gen_shared_subcomm, gen_shared_scatter_ctx, &
    getvecpointer, restorevecpointer
  use m_helper_functions, only : itoa, CHKERR, imp_allreduce_sum

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    logical :: lpetsc_is_initialized
    integer(mpiint) :: mpierr
    call PetscInitialized(lpetsc_is_initialized, mpierr)
    if(lpetsc_is_initialized) call PetscFinalize(mpierr)
  end subroutine teardown


  @test(npes =[3])
  subroutine test_gen_shared_subcomm(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, numnodes, myid
    integer(mpiint) :: subcomm, subnumnodes, submyid
    integer(mpiint) :: ierr

    integer(iintegers) :: k

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    if(myid.eq.0) print *,'Main Comm:'
    do k=0,numnodes-1
      if(myid.eq.k) print *,'Hi, from main comm: this is',myid,'out of ',numnodes
      call mpi_barrier(comm,ierr); call CHKERR(ierr)
    enddo
    call mpi_barrier(comm,ierr); call CHKERR(ierr)

    call gen_shared_subcomm(comm, subcomm, ierr); call CHKERR(ierr)

    call mpi_comm_rank(subcomm,submyid,ierr); call CHKERR(ierr)
    call mpi_comm_size(subcomm,subnumnodes,ierr); call CHKERR(ierr)


    if(myid.eq.0) print *,'Sub Comm:'
    do k=0,subnumnodes-1
      if(submyid.eq.k) print *,'Hi, from sub comm: this is',myid,'out of ',numnodes,' sub:', submyid, 'of', subnumnodes
      call mpi_barrier(subcomm,ierr); call CHKERR(ierr)
    enddo
    call mpi_barrier(comm,ierr); call CHKERR(ierr)
  end subroutine


  @test(npes =[1,2,3])
  subroutine test_scatter_to_shared_subcomm(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, numnodes, myid
    integer(mpiint) :: subcomm, subnumnodes, submyid
    integer(mpiint) :: ierr

    integer(iintegers),parameter :: nxp=3,nyp=3,nv=2
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: sundir(3) = 0
    real(ireals),parameter :: dz=dx
    real(ireals) :: dz1d(nv)

    type(t_solver_3_10) :: solver
    type(tVec) :: gvec, svec, gvec_target
    real(ireals), pointer :: x1d(:)=>NULL(), xv(:,:,:,:)=>NULL()

    integer(iintegers) :: N, k, num_shared_masters
    type(tVecScatter) :: ctx
    real(ireals) :: diff

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call gen_shared_subcomm(comm, subcomm, ierr); call CHKERR(ierr)

    call mpi_comm_rank(subcomm,submyid,ierr); call CHKERR(ierr)
    call mpi_comm_size(subcomm,subnumnodes,ierr); call CHKERR(ierr)

    call init_pprts(comm, nv, nxp, nyp, dx,dy, sundir, solver, dz1d)

    associate(C=>solver%C_one)
      ! Gen new global vector
      call DMGetGlobalVector(C%da, gvec, ierr); call CHKERR(ierr)

      call getVecPointer(C%da, gvec, x1d, xv)
      x1d(:) = 100._ireals + real(myid, ireals)
      call restoreVecPointer(C%da, gvec, x1d, xv)
      call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_gvec', ierr); call CHKERR(ierr)

      if(submyid.eq.0) then
        call VecGetSize(gvec, N, ierr); call CHKERR(ierr)
      else
        N = 0
      endif
      call VecCreateSeq(PETSC_COMM_SELF, N, svec, ierr); call CHKERR(ierr)
      call VecSetFromOptions(svec, ierr); call CHKERR(ierr)
      call PetscObjectSetName(svec, 'shared_mem_vec_rank_'//itoa(myid), ierr);call CHKERR(ierr)

      call gen_shared_scatter_ctx(gvec, svec, ctx, ierr); call CHKERR(ierr)

      call VecScatterBegin(ctx, gvec, svec, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
      call VecScatterEnd  (ctx, gvec, svec, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

      do k=0, numnodes-1
        if(k.eq.myid) then
          call PetscObjectViewFromOptions(svec, PETSC_NULL_VEC, '-show_svec', ierr); call CHKERR(ierr)
        endif
        call mpi_barrier(comm, ierr); call CHKERR(ierr)
      enddo
      call mpi_barrier(comm, ierr); call CHKERR(ierr)


      ! And read back the average values
      call VecDuplicate(gvec, gvec_target, ierr); call CHKERR(ierr)
      call VecCopy(gvec, gvec_target, ierr); call CHKERR(ierr)
      call VecSet(gvec, zero, ierr); call CHKERR(ierr)

      call VecScatterBegin(ctx, svec, gvec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
      call VecScatterEnd  (ctx, svec, gvec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

      if(submyid.eq.0) then
        k = 1
      else
        k = 0
      endif
      call imp_allreduce_sum(comm, k, num_shared_masters)
      call VecScale(gvec, one/real(num_shared_masters, ireals), ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_gvec', ierr); call CHKERR(ierr)

      call VecAXPY(gvec, -one, gvec_target, ierr); call CHKERR(ierr)
      call VecNorm(gvec, NORM_1, diff, ierr); call CHKERR(ierr)

      @assertEqual(zero, diff, 10*epsilon(diff), "difference should be zero between vec before sending it to shared comm and after averaging")
    end associate

    ! Tidy up
    call VecScatterDestroy(ctx, ierr); call CHKERR(ierr)
    call DMRestoreGlobalVector(solver%C_one%da, gvec, ierr); call CHKERR(ierr)
    call destroy_pprts(solver, lfinalizepetsc=.True.)

  end subroutine

end module
