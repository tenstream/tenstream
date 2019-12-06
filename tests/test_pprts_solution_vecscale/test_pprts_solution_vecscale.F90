module test_pprts_solution_vecscale

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, irealLUT, zero, one, pi, mpiint

#include "petsc/finclude/petsc.h"
  use petsc

  use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_8_10, t_solver_8_16, &
    prepare_solution, print_solution
  use m_pprts, only : init_pprts, scale_flx, destroy_pprts
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
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    if(solver_3_10%linitialized) call destroy_pprts(solver_3_10, lfinalizepetsc=.False.)
    if(solver_8_10%linitialized) call destroy_pprts(solver_8_10, lfinalizepetsc=.False.)
    if(solver_8_16%linitialized) call destroy_pprts(solver_8_16, lfinalizepetsc=.False.)
    call PetscFinalize(ierr) ;call CHKERR(ierr)
  end subroutine teardown

  @test(npes =[2,1])
  subroutine test_pprts_solution_vecscale_back_and_forth(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers),parameter :: nxp=5,nyp=5,nv=5,uid=0
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=10, theta0=60
    real(ireals),parameter :: albedo=0., dz=dx
    real(ireals),parameter :: incSolar=1000
    real(ireals),parameter :: atolerance=.1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir0,fdn0,fup0,fdiv0
    real(ireals),allocatable,dimension(:,:,:) :: fdir1,fdn1,fup1,fdiv1
    real(ireals),allocatable,dimension(:,:,:) :: fdir2,fdn2,fup2,fdiv2
    real(ireals),allocatable,dimension(:,:,:) :: fdir3,fdn3,fup3,fdiv3

    integer(iintegers) :: j
    integer(iintegers) :: cx, cy      ! global indices of cloud

    dz1d = dz
    dz1d(1) = 10*dx

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call this_test(solver_3_10)
    call this_test(solver_8_10)
    call this_test(solver_8_16)

    contains
      subroutine this_test(solver)
      class(t_solver), intent(inout) :: solver
        real(ireals) :: target_dirnorm, target_diffnorm
        real(ireals) :: dirnorm, diffnorm
        integer(mpiint) :: ierr
        call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

        associate( S => solver%solutions(uid) )
          call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
            lsolar=.True., solution=S, uid=uid)

          call VecSet(S%edir , one, ierr); call CHKERR(ierr)
          S%lWm2_dir  = .True.
          call VecSet(S%ediff, one, ierr); call CHKERR(ierr)
          S%lWm2_diff = .True.

          call VecNorm(S%edir , NORM_2, target_dirnorm , ierr); call CHKERR(ierr)
          call VecNorm(S%ediff, NORM_2, target_diffnorm, ierr); call CHKERR(ierr)

          call print_solution(S)

          call scale_flx(solver, S, lWm2=.False. )

          call print_solution(S)

          call scale_flx(solver, S, lWm2=.True. )

          call print_solution(S)

          call VecNorm(S%edir , NORM_2, dirnorm , ierr); call CHKERR(ierr)
          call VecNorm(S%ediff, NORM_2, diffnorm, ierr); call CHKERR(ierr)

          @assertEqual(target_dirnorm, dirnorm, epsilon(dirnorm), 'direct radiation vec norm changed when scaling back and forth')
          @assertEqual(target_diffnorm, diffnorm, epsilon(diffnorm), 'diffuse radiation vec norm changed when scaling back and forth')
        end associate
    end subroutine
  end subroutine
end module
