@test(npes=[9, 6])
subroutine test_tenstream_ex1(this)

  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint
  use m_pprts_base, only: &
    & allocate_pprts_solver_from_commandline, &
    & destroy_pprts, &
    & t_coord, &
    & t_solver
  use m_pprts, only: init_pprts
  use m_helper_functions, only: reorder_mpi_comm, chkerr
  use pfunit_mod
!!#ifdef HAVE_PETSC
!!#include "petsc/finclude/petsc.h"
!!  use petsc
!!#endif

  implicit none

  class(MpiTestMethod), intent(inout) :: this

  integer(iintegers), parameter :: nv = 3
  real(ireals), parameter :: dx = 10000, dy = dx
  real(ireals), parameter :: sundir(3) = 0
  real(ireals), parameter :: dz = dx / 1000
  real(ireals) :: dz1d(nv)

  integer(iintegers) :: nxp, nyp
  integer(mpiint) :: Nrank_x, Nrank_y

  integer(iintegers) :: neighbors_orig(4), neighbors_reorder(4)
  integer(mpiint) :: myid, mpierr, orig_id, numnodes

  class(t_solver), allocatable :: solver

  integer(mpiint) :: test_comm, comm, reorder_comm

  dz1d = dz

  test_comm = this%getMpiCommunicator()
  numnodes = this%getNumProcesses()

  call MPI_Comm_dup(test_comm, comm, mpierr); call CHKERR(mpierr)

  call mpi_comm_rank(comm, orig_id, mpierr)

  ! Select domain size and processor grid based on rank count.
  if (numnodes .eq. 9) then
    nxp = 9; nyp = 9
    Nrank_x = 3; Nrank_y = 3
  else  ! 6 ranks: nxp=2, nyp=3 (2x3 grid, tests asymmetric decomposition)
    nxp = 6; nyp = 6
    Nrank_x = 2; Nrank_y = 3
  end if

  call init_mpi_data_parameters(comm)
  call allocate_pprts_solver_from_commandline(solver, '3_10', mpierr); call CHKERR(mpierr)

  call init_pprts(comm, nv, nxp, nyp, dx, dy, sundir, solver, dz1d=dz1d)
  call mpi_comm_rank(comm, myid, mpierr)
  print *, 'I am originally', orig_id, 'my rank is ', myid, &
    ' and Neighbors are', solver%C_one%neighbors([10, 4, 16, 22])

  neighbors_orig = solver%C_one%neighbors([10, 4, 16, 22])
  call destroy_pprts(solver, .true.)

  call mpi_barrier(comm, mpierr)

  if (myid .eq. 0) print *, 'Reordering Communicator now!'
  call reorder_mpi_comm(comm, Nrank_x, Nrank_y, reorder_comm)

  call init_mpi_data_parameters(reorder_comm)
  call init_pprts(reorder_comm, nv, nxp, nyp, dx, dy, sundir, solver, dz1d=dz1d)
  call mpi_comm_rank(reorder_comm, myid, mpierr)
  print *, 'I am originally', orig_id, 'my rank is now', myid, &
    ' and Neighbors are', solver%C_one%neighbors([10, 4, 16, 22])

  neighbors_reorder = solver%C_one%neighbors([10, 4, 16, 22])
  call destroy_pprts(solver, .true.)

  ! Now check if the neighbors fit.
  ! Column-major rank numbering: rank = xi + yi*nxp (xi fastest, matches PETSc DMDA).
  ! neighbors([10, 4, 16, 22]) = [west, south, east, north].
  ! And the ordering of petsc numbering is clockwise, starting at the bottom for original neighbors
  ! and anti-clockwise starting with left neighbor on the reordered grid

  if (numnodes .eq. 9) then
    ! 3x3 grid:
    !   original (rank = xi + yi*3):      reordered (rank = xi + yi*3 after reorder):
    !   2  5  8                            6  7  8
    !   1  4  7                            3  4  5
    !   0  3  6                            0  1  2
    !
    ! Ranks 0, 4, 8 stay at the same spatial position after reordering.
    if (any(orig_id .eq. [0, 4, 8])) &
      @assertEqual(neighbors_orig, neighbors_reorder, '9-rank: diagonal corners keep same neighbors')

    if (orig_id .eq. 0) then
      @assertEqual([2, 6, 1, 3], neighbors_orig)
      @assertEqual([2, 6, 1, 3], neighbors_reorder)
    end if

    if (orig_id .eq. 1) then
      @assertEqual([0, 7, 2, 4], neighbors_orig)
      @assertEqual([5, 0, 4, 6], neighbors_reorder)
    end if
    if (orig_id .eq. 2) then
      @assertEqual([1, 8, 0, 5], neighbors_orig)
      @assertEqual([8, 3, 7, 0], neighbors_reorder)
    end if
    if (orig_id .eq. 3) then
      @assertEqual([5, 0, 4, 6], neighbors_orig)
      @assertEqual([0, 7, 2, 4], neighbors_reorder)
    end if
  else
    ! 6 ranks, 2x3 grid (Nrank_x=2, Nrank_y=3):
    !   original (rank = xi*3 + yi):   reordered (rank = xi + yi*2):
    !   2  5    (yi=2)                  4  5
    !   1  4    (yi=1)                  2  3
    !   0  3    (yi=0)                  0  1
    !   xi: 0  1
    !
    if (any(orig_id .eq. [0, 5])) &
      @assertEqual(neighbors_orig, neighbors_reorder, '6-rank: diagonal corners keep same neighbors')

    if (orig_id .eq. 0) then
      @assertEqual([1, 4, 1, 2], neighbors_orig)
      @assertEqual([1, 4, 1, 2], neighbors_reorder)
    end if
    if (orig_id .eq. 1) then
      @assertEqual([0, 5, 0, 3], neighbors_orig)
      @assertEqual([3, 0, 3, 4], neighbors_reorder)
    end if
    if (orig_id .eq. 2) then
      @assertEqual([3, 0, 3, 4], neighbors_orig)
      @assertEqual([5, 2, 5, 0], neighbors_reorder)
    end if
    if (orig_id .eq. 3) then
      @assertEqual([2, 1, 2, 5], neighbors_orig)
      @assertEqual([0, 5, 0, 3], neighbors_reorder)
    end if
    if (orig_id .eq. 4) then
      @assertEqual([5, 2, 5, 0], neighbors_orig)
      @assertEqual([2, 1, 2, 5], neighbors_reorder)
    end if
  end if

end subroutine
