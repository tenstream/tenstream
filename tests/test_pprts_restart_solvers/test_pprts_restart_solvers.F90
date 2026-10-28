module test_pprts_restart_solvers
  ! Explicit diffuse solvers for a cold solve followed by a warm restart, as in LES runs:
  ! the first call converges with the default (tight) tolerance, the second one starts from the previous solution
  ! after the cloud moved by one column and solves for -solar_diff_ksp_rtol.
  !
  ! Compares SOR (-solar_diff_ksp_type richardson), BiCGStab (-solar_diff_ksp_type bcgs) and
  ! BiCGStab for the first solve only (-solar_diff_initial_ksp_type bcgs) and prints the cost of both calls.
  ! The domain size can be set with -test_restart_nxy and -test_restart_nz, e.g. to measure the cost on a larger scene.

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one

  use mpi, only: MPI_Wtime

  use pfunit_mod

  use m_helper_functions, only: &
    & CHKERR, &
    & get_petsc_opt, &
    & insert_petsc_opt, &
    & spherical_2_cartesian, &
    & toStr

  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result_toZero

  implicit none

  real(ireals), parameter :: dx = 100, dy = dx, dz = 100
  real(ireals), parameter :: theta0 = 40, phi0 = 30
  real(ireals), parameter :: incSolar = 1000
  real(ireals), parameter :: albedo = .2_ireals
  real(ireals), parameter :: rtol_restart = 1e-2_ireals
  integer(iintegers), parameter :: pc_sor_its = 2

  ! results of one configuration: iterations and wall time of both calls, fluxes on rank 0
  type t_result
    integer(iintegers) :: niter(2)
    real(ireals) :: time(2)
    real(ireals), allocatable, dimension(:, :, :) :: edn1, eup1, edn2, eup2
  end type

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    call insert_petsc_opt('-solar_diff_ksp_type richardson -solar_diff_initial_ksp_type richardson', ierr); call CHKERR(ierr)
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes=[2, 1])
  subroutine restart_sor_vs_bcgs(this)
    class(MpiTestMethod), intent(inout) :: this
    type(t_result) :: sor, bcgs, initial_bcgs
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    call solve_twice(comm, 'richardson', 'richardson', sor)
    call solve_twice(comm, 'bcgs', 'bcgs', bcgs)
    call solve_twice(comm, 'richardson', 'bcgs', initial_bcgs)

    if (myid .eq. 0) then
      call print_result('sor', sor, .false., .false.)
      call print_result('bcgs', bcgs, .true., .true.)
      call print_result('initial_bcgs', initial_bcgs, .true., .false.)
    end if

    ! the first solves are converged tightly, i.e. all solvers have to agree closely
    if (myid .eq. 0) then
      @assertLessThan(rel_diff(bcgs%edn1, sor%edn1), 1e-3_ireals, 'bcgs edn, first call')
      @assertLessThan(rel_diff(bcgs%eup1, sor%eup1), 1e-3_ireals, 'bcgs eup, first call')
      @assertLessThan(rel_diff(initial_bcgs%edn1, sor%edn1), 1e-3_ireals, 'initial_bcgs edn, first call')
      @assertLessThan(rel_diff(initial_bcgs%eup1, sor%eup1), 1e-3_ireals, 'initial_bcgs eup, first call')
      ! the restarts only for rtol_restart
      @assertLessThan(rel_diff(bcgs%edn2, sor%edn2), 2e-2_ireals, 'bcgs edn, restart')
      @assertLessThan(rel_diff(bcgs%eup2, sor%eup2), 2e-2_ireals, 'bcgs eup, restart')
      @assertLessThan(rel_diff(initial_bcgs%edn2, sor%edn2), 2e-2_ireals, 'initial_bcgs edn, restart')
      @assertLessThan(rel_diff(initial_bcgs%eup2, sor%eup2), 2e-2_ireals, 'initial_bcgs eup, restart')
    end if

    ! -solar_diff_initial_ksp_type bcgs has to do the very same first solve as -solar_diff_ksp_type bcgs
    @assertEqual(bcgs%niter(1), initial_bcgs%niter(1), 'initial_ksp_type: first solve is not done with bcgs')
  end subroutine

  ! cold solve, move the cloud by one column and solve again starting from the previous solution
  subroutine solve_twice(comm, ksp_type, initial_ksp_type, res)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: ksp_type, initial_ksp_type
    type(t_result), intent(out) :: res

    class(t_solver), allocatable :: solver
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, edn, eup, abso, edir
    real(ireals), allocatable :: dz1d(:)
    integer(iintegers) :: Nxy, Nz, icall
    real(kind(1d0)) :: tstart
    logical :: lflg
    integer(mpiint) :: ierr

    Nxy = 16
    call get_petsc_opt('', '-test_restart_nxy', Nxy, lflg, ierr); call CHKERR(ierr)
    Nz = 12
    call get_petsc_opt('', '-test_restart_nz', Nz, lflg, ierr); call CHKERR(ierr)

    ! always set all of the options, they stay in the options database.
    ! The test is about the explicit solvers, override e.g. PETSC_OPTIONS of the CI implicit job
    call insert_petsc_opt('-solar_dir_explicit -solar_diff_explicit', ierr); call CHKERR(ierr)
    call insert_petsc_opt('-solar_diff_ksp_type '//ksp_type, ierr); call CHKERR(ierr)
    call insert_petsc_opt('-solar_diff_initial_ksp_type '//initial_ksp_type, ierr); call CHKERR(ierr)
    call insert_petsc_opt('-solar_diff_ksp_rtol '//toStr(rtol_restart), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-solar_diff_ksp_atol 0', ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)
    allocate (dz1d(Nz), source=dz)
    call init_pprts(comm, Nz, Nxy, Nxy, dx, dy, spherical_2_cartesian(phi0, theta0), solver, dz1d=dz1d)

    associate (C => solver%C_one)
      allocate (kabs(C%zm, C%xm, C%ym), ksca(C%zm, C%xm, C%ym), g(C%zm, C%xm, C%ym))
    end associate

    do icall = 1, 2
      call cloud(solver, Nxy, Nz, shift=icall - 1, kabs=kabs, ksca=ksca, g=g)
      call set_optical_properties(solver, albedo, kabs, ksca, g)
      tstart = MPI_Wtime()
      call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=incSolar, opt_solution_uid=0_iintegers)
      res%time(icall) = real(MPI_Wtime() - tstart, ireals)
      res%niter(icall) = solver%solutions(0)%Niter_diff
      call pprts_get_result_toZero(solver, edn, eup, abso, edir, opt_solution_uid=0_iintegers)
      if (icall .eq. 1) then
        if (allocated(edn)) allocate (res%edn1, source=edn)
        if (allocated(eup)) allocate (res%eup1, source=eup)
      else
        if (allocated(edn)) allocate (res%edn2, source=edn)
        if (allocated(eup)) allocate (res%eup2, source=eup)
      end if
    end do

    call destroy_pprts(solver, lfinalizepetsc=.false.)
  end subroutine

  ! smooth, optically thick cloud of about 3/4 of the domain size, centered at the global column (Nxy/2 + shift, Nxy/2)
  subroutine cloud(solver, Nxy, Nz, shift, kabs, ksca, g)
    class(t_solver), intent(in) :: solver
    integer(iintegers), intent(in) :: Nxy, Nz, shift
    real(ireals), intent(out), dimension(:, :, :) :: kabs, ksca, g
    real(ireals), parameter :: ksca_cloud = 5._ireals / dz, ksca_clear = 1e-5_ireals
    real(ireals) :: r2, w
    integer(iintegers) :: i, j, k, gi, gj
    associate (C => solver%C_one)
      do j = 1, C%ym
        do i = 1, C%xm
          gi = C%xs + i
          gj = C%ys + j
          do k = 1, C%zm
            r2 = (real(gi - Nxy / 2 - shift, ireals) / (Nxy*.375_ireals))**2 &
               & + (real(gj - Nxy / 2, ireals) / (Nxy*.375_ireals))**2 &
               & + (real(k - Nz / 2, ireals) / (Nz*.25_ireals))**2
            w = max(zero, one - r2)
            ksca(k, i, j) = ksca_clear + ksca_cloud * w
            kabs(k, i, j) = 1e-6_ireals + 1e-4_ireals * ksca_cloud * w
            g(k, i, j) = .85_ireals * w
          end do
        end do
      end do
    end associate
  end subroutine

  ! max absolute difference relative to the max of the reference
  pure function rel_diff(a, ref) result(d)
    real(ireals), intent(in) :: a(:, :, :), ref(:, :, :)
    real(ireals) :: d
    d = maxval(abs(a - ref)) / max(maxval(abs(ref)), tiny(d))
  end function

  ! number of SOR sweeps: one per SOR iteration, BiCGStab needs 2*pc_sor_its per iteration plus pc_sor_its to start
  subroutine print_result(name, res, lbcgs_first, lbcgs_restart)
    character(len=*), intent(in) :: name
    type(t_result), intent(in) :: res
    logical, intent(in) :: lbcgs_first, lbcgs_restart
    print *, 'restart test '//name//' :: first call '//toStr(res%niter(1))//' iterations, '// &
      & toStr(sweeps(res%niter(1), lbcgs_first))//' sweeps, '//toStr(res%time(1))//' s'// &
      & ' :: restart '//toStr(res%niter(2))//' iterations, '// &
      & toStr(sweeps(res%niter(2), lbcgs_restart))//' sweeps, '//toStr(res%time(2))//' s'
  contains
    integer(iintegers) function sweeps(niter, lbcgs)
      integer(iintegers), intent(in) :: niter
      logical, intent(in) :: lbcgs
      if (lbcgs) then
        sweeps = pc_sor_its * (1 + 2 * niter)
      else
        sweeps = niter
      end if
    end function
  end subroutine
end module
