module m_examples_pprts_box_cld

  use m_data_parameters, only: init_mpi_data_parameters, finalize_mpi, &
    & iintegers, ireals, mpiint, zero, pi, imp_ireals

  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, &
                     pprts_get_result
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts

  use m_helper_functions, only: CHKERR, get_mem_footprint, spherical_2_cartesian, is_inrange, get_petsc_opt

  use mpi, only: MPI_COMM_WORLD, MPI_SUM, MPI_IN_PLACE

  implicit none

  integer(mpiint) :: myid, ierr

contains
  subroutine ex_pprts_box_cld()
    implicit none

    integer(iintegers), parameter :: nxp = 17, nyp = 33, nv = 20
    real(ireals), parameter :: dx = 200, dy = dx
    real(ireals), parameter :: dz = 100
    real(ireals), parameter :: phi0 = 0, theta0 = 70
    real(ireals), parameter :: incSolar = 1364
    integer(iintegers) :: cld_width, cld_k(2), nk
    real(ireals) :: cld_tau, cld_w0, cld_g
    real(ireals) :: dz1d(nv)

    real(ireals) :: sundir(3), albedo
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, plck
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    integer(iintegers) :: xl, xu, yl, yu, i, j, k
    real(ireals), allocatable :: prof(:, :) ! (layer, [domain sum, cloud column sum, ncol domain, ncol cloud])

    class(t_solver), allocatable :: solver
    logical :: lthermal, lsolar, lflg

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(MPI_COMM_WORLD)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)
    dz1d = dz

    sundir = spherical_2_cartesian(phi0, theta0)

    call init_pprts(MPI_COMM_WORLD, nv, nxp, nyp, dx, dy, sundir, solver, dz1d=dz1d)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

    allocate (kabs(solver%C_one_atm%zm, solver%C_one_atm%xm, solver%C_one_atm%ym))
    allocate (ksca(solver%C_one_atm%zm, solver%C_one_atm%xm, solver%C_one_atm%ym))
    allocate (g(solver%C_one_atm%zm, solver%C_one_atm%xm, solver%C_one_atm%ym))
    allocate (plck(solver%C_one_atm1%zm, solver%C_one_atm1%xm, solver%C_one_atm1%ym))

    kabs = .001_ireals / (dz * nv)
    ksca = .005_ireals / (dz * nv)
    g = zero
    plck = 100._ireals / pi

    ! cloud box: horizontal half width [px], vertical layer range, total optical thickness, ssa and asymmetry
    cld_width = 5
    call get_petsc_opt('', '-cld_width', cld_width, lflg, ierr); call CHKERR(ierr)
    cld_k = nv / 2
    nk = 2
    call get_petsc_opt('', '-cld_k', cld_k, nk, lflg, ierr); call CHKERR(ierr)
    if (lflg .and. nk .ne. 2) call CHKERR(1_mpiint, '-cld_k needs two values: first and last cloud layer')
    cld_tau = 10 + 1.001_ireals / real(nv, ireals)
    call get_petsc_opt('', '-cld_tau', cld_tau, lflg, ierr); call CHKERR(ierr)
    cld_w0 = 10 / cld_tau
    call get_petsc_opt('', '-cld_w0', cld_w0, lflg, ierr); call CHKERR(ierr)
    cld_g = .9
    call get_petsc_opt('', '-cld_g', cld_g, lflg, ierr); call CHKERR(ierr)

    xl = int((nxp + 1._ireals) / 2._ireals) - cld_width
    xu = int((nxp + 1._ireals) / 2._ireals) + cld_width
    yl = int((nyp + 1._ireals) / 2._ireals) - cld_width
    yu = int((nyp + 1._ireals) / 2._ireals) + cld_width
    if (myid .eq. 0) print *, 'cld between x:', xl, xu, 'y:', yl, yu

    do j = 1_iintegers + solver%C_one%ys, solver%C_one%ye + 1
      do i = 1_iintegers + solver%C_one%xs, solver%C_one%xe + 1
        if (is_inrange(i, xl, xu) .and. is_inrange(j, yl, yu)) then
          kabs(cld_k(1):cld_k(2), i - solver%C_one%xs, j - solver%C_one%ys) = &
            & cld_tau * (1 - cld_w0) / (dz * real(cld_k(2) - cld_k(1) + 1, ireals))
          ksca(cld_k(1):cld_k(2), i - solver%C_one%xs, j - solver%C_one%ys) = &
            & cld_tau * cld_w0 / (dz * real(cld_k(2) - cld_k(1) + 1, ireals))
          g(cld_k(1):cld_k(2), i - solver%C_one%xs, j - solver%C_one%ys) = cld_g
        end if
      end do
    end do

    albedo = .1
    call get_petsc_opt('', '-Ag', albedo, lflg, ierr); call CHKERR(ierr)

    call set_optical_properties(solver, albedo, kabs, ksca, g, plck)

    lsolar = .true.
    lthermal = .true.
    call get_petsc_opt('', '-solar', lsolar, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-thermal', lthermal, lflg, ierr); call CHKERR(ierr)

    call solve_pprts(solver, lthermal=lthermal, lsolar=lsolar, edirTOA=incSolar)

    allocate (fdir(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate (fdn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate (fup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate (fdiv(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    call pprts_get_result(solver, fdn, fup, fdiv, fdir)

    ! domain mean and cloud column mean of absorption [W/m3] per layer and of TOA eup
    allocate (prof(0:size(fdiv, 1), 4), source=0._ireals)
    do j = 1, size(fdiv, 3)
      do i = 1, size(fdiv, 2)
        prof(1:, 1) = prof(1:, 1) + fdiv(:, i, j)
        prof(0, 1) = prof(0, 1) + fup(1, i, j)
        prof(:, 3) = prof(:, 3) + 1
        if (is_inrange(i + solver%C_one%xs, xl, xu) .and. is_inrange(j + solver%C_one%ys, yl, yu)) then
          prof(1:, 2) = prof(1:, 2) + fdiv(:, i, j)
          prof(0, 2) = prof(0, 2) + fup(1, i, j)
          prof(:, 4) = prof(:, 4) + 1
        end if
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, prof, int(size(prof), mpiint), imp_ireals, MPI_SUM, MPI_COMM_WORLD, ierr); call CHKERR(ierr)
    if (myid .eq. 0) then
      print *, 'TOA eup domain mean, cloud column mean', prof(0, 1) / prof(0, 3), prof(0, 2) / max(1._ireals, prof(0, 4))
      do k = 1, size(fdiv, 1)
        print *, 'abso_profile', k, prof(k, 1) / prof(k, 3), prof(k, 2) / max(1._ireals, prof(k, 4))
      end do
    end if

    if (myid .eq. 0) then
      print *, 'kabs:', kabs(:, 1, 1)
      print *, 'fdir:', fdir(:, 1, 1)
      print *, 'fdn:', fdn(:, 1, 1)
      print *, 'fup:', fup(:, 1, 1)
      print *, 'fdiv:', fdiv(:, 1, 1)
    end if
    print *, myid, 'Memory:', get_mem_footprint(MPI_COMM_WORLD)

    call destroy_pprts(solver, .true.)
    call finalize_mpi(MPI_COMM_WORLD, .true., .true.)
  end subroutine
end module
