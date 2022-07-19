module test_petsc_scatterToZero

#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only: iintegers, ireals, mpiint, zero, one
  use m_pprts_base, only: t_solver_3_10, destroy_pprts
  use m_pprts, only: init_pprts
  use m_petsc_helpers, only: petscVecToF90, petscGlobalVecToZero, f90VecToPetsc
  use m_helper_functions, only: CHKERR

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver
  type(tVec) :: gvec, lvec

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    ! Tidy up
    call DMRestoreGlobalVector(solver%C_one%da, gvec, ierr); call CHKERR(ierr)
    call VecDestroy(lvec, ierr); call CHKERR(ierr)
    call destroy_pprts(solver, lfinalizepetsc=.true.)
  end subroutine teardown

  @test(npes=[2, 3])
  subroutine petsc_scatterToZero(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, ierr

    integer(iintegers), parameter :: nxp = 9, nyp = 9, nv = 2
    real(ireals), parameter :: dx = 100, dy = dx
    real(ireals), parameter :: sundir(3) = 0
    real(ireals), parameter :: dz = dx
    real(ireals) :: dz1d(nv)

    real(ireals), allocatable, dimension(:, :, :, :) :: local_arr
    real(ireals), allocatable, dimension(:, :, :, :) :: local_arr_2
    real(ireals), allocatable, dimension(:, :, :, :) :: global_arr_on_rank0

    integer(iintegers) :: i, j, k, d

    dz1d = dz

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_pprts(comm, nv, nxp, nyp, dx, dy, sundir, solver, dz1d)

    associate (C => solver%C_one)
      allocate (local_arr(C%dof, C%zm, C%xm, C%ym))

      local_arr = myid

      ! Gen new global vector
      call DMGetGlobalVector(C%da, gvec, ierr); call CHKERR(ierr)

      ! copy local fortran array into global vec
      call f90VecToPetsc(local_arr, C%da, gvec)

      ! Check that we can get an local array from a global petsc vec
      call petscVecToF90(gvec, C%da, local_arr_2)
      call assert_equivalence(local_arr, local_arr_2, 'get local arr from global petsc vec')

      deallocate (local_arr_2)

      ! do that feat again with a 3D array
      call f90VecToPetsc(local_arr(:, :, :, :), C%da, gvec)
      call petscVecToF90(gvec, C%da, local_arr_2)
      call assert_equivalence(local_arr, local_arr_2, 'get local arr from global petsc vec in 3D')

      ! Debug Output
      call PetscObjectSetName(gvec, 'VecGlobal', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_gvec', ierr); call CHKERR(ierr)

      ! copy global vec to rank 0, lVec(full size)
      call petscGlobalVecToZero(gVec, C%da, lVec)

      ! Debug Output
      call PetscObjectSetName(lvec, 'VecLocal', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(lvec, PETSC_NULL_VEC, '-show_lvec', ierr); call CHKERR(ierr)

      ! Put data from petsc vec into fortran array
      if (myid .eq. 0) then
        call petscVecToF90(lvec, C%da, global_arr_on_rank0, only_on_rank0=.true.)
        do j = 1, ubound(global_arr_on_rank0, 4)
          do i = 1, ubound(global_arr_on_rank0, 3)
            do k = 1, ubound(global_arr_on_rank0, 2)
              do d = 1, ubound(global_arr_on_rank0, 1)
                @assertEqual(real((j - 1) / C%ym), global_arr_on_rank0(d, k, i, j), '')
              end do
            end do
          end do
        end do
      end if
    end associate
  end subroutine

  subroutine assert_equivalence(a, b, msg)
    real(ireals), dimension(:, :, :, :) :: a, b
    character(len=*), intent(in) :: msg
    integer(iintegers) :: i, j, k, d
    do j = 1, ubound(a, 4)
      do i = 1, ubound(a, 3)
        do k = 1, ubound(a, 2)
          do d = 1, ubound(a, 1)
            @assertEqual(a(d, k, i, j), b(d, k, i, j), trim(msg))
          end do
        end do
      end do
    end do
  end subroutine

end module
