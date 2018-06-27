@test(npes =[9])
subroutine test_tenstream_ex1(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint
    use m_pprts_base, only : t_coord, t_solver_3_10
    use m_pprts, only : init_pprts, destroy_pprts
    use m_helper_functions, only : reorder_mpi_comm, chkerr
    use m_tenstream_options, only: read_commandline_options
    use pfunit_mod
#include "petsc/finclude/petsc.h"
    use petsc

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers),parameter :: nxp=9,nyp=9,nv=3
    real(ireals),parameter :: dx=10000,dy=dx
    real(ireals),parameter :: phi0=-1, theta0=-1
    real(ireals),parameter :: dz=dx/1000
    real(ireals) :: dz1d(nv)

    integer(mpiint),parameter :: Nrank_x=3, Nrank_y=3

    integer(iintegers) :: neighbors_orig(4), neighbors_reorder(4)
    integer(mpiint) :: myid, mpierr, orig_id, numnodes

    type(t_solver_3_10) :: solver

    MPI_Comm :: comm, reorder_comm

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()

    call mpi_comm_rank( comm, orig_id, mpierr)

    call init_mpi_data_parameters(comm)

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, "-twostr_only", mpierr); call CHKERR(mpierr)
    call init_pprts(comm, nv, nxp, nyp, dx, dy, phi0, theta0, solver, dz1d=dz1d)
    call mpi_comm_rank( comm, myid, mpierr)
    print *,'I am originally', orig_id, 'my rank is now', myid, &
      ' and Neighbors are', solver%C_one%neighbors([10,4,16,22])

    neighbors_orig = solver%C_one%neighbors([10,4,16,22])
    call destroy_pprts(solver, .True.)

    call mpi_barrier(comm, mpierr)

    if (myid.eq.0) print *,'Reordering Communicator now!'
    call reorder_mpi_comm(comm, Nrank_x, Nrank_y, reorder_comm)

    call init_mpi_data_parameters(reorder_comm)
    call init_pprts(reorder_comm, nv, nxp, nyp, dx, dy, phi0, theta0, solver, dz1d=dz1d)
    call mpi_comm_rank( comm, myid, mpierr)
    print *,'I am originally', orig_id, 'my rank is now', myid, &
      ' and Neighbors are', solver%C_one%neighbors([10,4,16,22])

    neighbors_reorder = solver%C_one%neighbors([10,4,16,22])
    call destroy_pprts(solver, .True.)


    ! Now check if the neighbors fit:
    ! if we write it down, it should look like
    !
    ! original:
    !
    ! 2  5  8
    ! 1  4  7
    ! 0  3  6
    !
    ! reordered:
    !
    ! 6  7  8
    ! 3  4  5
    ! 0  1  2
    !
    ! And the ordering of petsc numbering is clockwise, starting at the bottom for original neighbors
    ! and anti-clockwise starting with left neighbor on the reordered grid

    if (any(orig_id.eq.[0,4,8])) &
      @assertEqual(neighbors_orig, neighbors_reorder, 'neighbors have to stay the same in the middle and for the corners')

    if (orig_id.eq.1) then
      @assertEqual([0,7,2,4], neighbors_orig)
      @assertEqual([5,0,4,6], neighbors_reorder)
    endif

    if (orig_id.eq.2) then
      @assertEqual([1,8,0,5], neighbors_orig)
      @assertEqual([8,3,7,0], neighbors_reorder)
    endif

    if (orig_id.eq.3) then
      @assertEqual([5,0,4,6], neighbors_orig)
      @assertEqual([0,7,2,4], neighbors_reorder)
    endif

end subroutine
