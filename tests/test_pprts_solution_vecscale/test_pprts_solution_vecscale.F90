module test_pprts_solution_vecscale

  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, irealLUT, zero, one, pi, mpiint

#include "petsc/finclude/petsc.h"
  use petsc

  use m_pprts_base, only: t_solver, t_solver_3_10, t_solver_8_10, t_solver_8_16, &
                          prepare_solution, print_solution, destroy_pprts
  use m_pprts, only: init_pprts, scale_flx
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: CHKERR

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver_3_10
  type(t_solver_8_10) :: solver_8_10
  type(t_solver_8_16) :: solver_8_16

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
    if (solver_3_10%linitialized) call destroy_pprts(solver_3_10, lfinalizepetsc=.false.)
    if (solver_8_10%linitialized) call destroy_pprts(solver_8_10, lfinalizepetsc=.false.)
    if (solver_8_16%linitialized) call destroy_pprts(solver_8_16, lfinalizepetsc=.false.)
    call PetscFinalize(ierr); call CHKERR(ierr)
  end subroutine teardown

  @test(npes=[2, 1])
  subroutine test_pprts_solution_vecscale_back_and_forth(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: nxp = 5, nyp = 5, nv = 5, uid = 0
    real(ireals), parameter :: dx = 100, dy = dx
    real(ireals), parameter :: sundir(3) = 0
    real(ireals), parameter :: dz = dx
    real(ireals) :: dz1d(nv)

    dz1d = dz
    dz1d(1) = 10 * dx

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call this_test(solver_3_10)
    call this_test(solver_8_10)
    call this_test(solver_8_16)

  contains
    subroutine this_test(solver)
      class(t_solver), intent(inout) :: solver
      real(ireals) :: target_dirnorm, target_diffnorm
      real(ireals) :: dirnorm, diffnorm
      integer(mpiint) :: ierr
      call init_pprts(comm, nv, nxp, nyp, dx, dy, sundir, solver, dz1d)

      associate (S => solver%solutions(uid))
        call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
                              lsolar=.true., lthermal=.false., solution=S, uid=uid)

        call VecSet(S%edir, one, ierr); call CHKERR(ierr)
        S%lWm2_dir = .true.
        call VecSet(S%ediff, one, ierr); call CHKERR(ierr)
        S%lWm2_diff = .true.

        call VecNorm(S%edir, NORM_2, target_dirnorm, ierr); call CHKERR(ierr)
        call VecNorm(S%ediff, NORM_2, target_diffnorm, ierr); call CHKERR(ierr)

        call print_solution(S)

        call scale_flx(solver, S, lWm2=.false.)

        call print_solution(S)

        call scale_flx(solver, S, lWm2=.true.)

        call print_solution(S)

        call VecNorm(S%edir, NORM_2, dirnorm, ierr); call CHKERR(ierr)
        call VecNorm(S%ediff, NORM_2, diffnorm, ierr); call CHKERR(ierr)

        @assertEqual(target_dirnorm, dirnorm, epsilon(dirnorm), 'direct radiation vec norm changed when scaling back and forth')
        @assertEqual(target_diffnorm, diffnorm, epsilon(diffnorm), 'diffuse radiation vec norm changed when scaling back and forth')
      end associate
    end subroutine
  end subroutine
end module
