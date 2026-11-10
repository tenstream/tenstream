!> @brief tests for the pprts MonteCarlo solver mcdmda
!> @details the checks only use properties that hold exactly for each MonteCarlo realization (not only in the mean),
!>          e.g. Beer-Lambert for the domain mean direct radiation without scattering or energy conservation,
!>          so that we can use tight tolerances with few photons.
module test_mcdmda

  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, pi

#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#endif

  use m_pprts_base, only: t_solver_mcdmda, destroy_pprts
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, set_angles, pprts_get_result
  use m_helper_functions, only: CHKERR, deg2rad, imp_allreduce_mean, spherical_2_cartesian, toStr

  use pfunit_mod

  implicit none

  integer(iintegers), parameter :: Nx = 4, Ny = 4, Nlay = 3
  real(ireals), parameter :: dx = 100, dy = dx, dz = 100
  real(ireals), parameter :: incSolar = 1000

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    continue
  end subroutine teardown

  !> @brief run mcdmda on a horizontally homogeneous atmosphere with Nlay layers
  !> @details planck is always given to the solver, also for solar computations,
  !>          i.e. solar results must not be affected by thermal emission
  subroutine run_mcdmda(comm, lsolar, theta0, albedo, dtau, w0, g, B, fdir, fdn, fup, fdiv)
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lsolar
    real(ireals), intent(in) :: theta0, albedo, dtau, w0, g, B ! dtau is the total optical thickness of the column
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: fdir, fdn, fup, fdiv

    type(t_solver_mcdmda) :: solver
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, gg, planck
    real(ireals), allocatable, dimension(:, :) :: planck_srfc
    real(ireals) :: dz1d(Nlay), sundir(3)

    dz1d = dz
    sundir = spherical_2_cartesian(0._ireals, theta0)

    call init_pprts(comm, Nlay, Nx, Ny, dx, dy, sundir, solver, dz1d)

    associate (Catm => solver%C_one_atm, Catm1 => solver%C_one_atm1)
      allocate (kabs(Catm%zm, Catm%xm, Catm%ym), source=dtau / (dz * Nlay) * (1 - w0))
      allocate (ksca(Catm%zm, Catm%xm, Catm%ym), source=dtau / (dz * Nlay) * w0)
      allocate (gg(Catm%zm, Catm%xm, Catm%ym), source=g)
      allocate (planck(Catm1%zm, Catm1%xm, Catm1%ym), source=B)
      allocate (planck_srfc(Catm1%xm, Catm1%ym), source=B)
    end associate

    call set_optical_properties(solver, albedo, kabs, ksca, gg, planck=planck, planck_srfc=planck_srfc, ldelta_scaling=.false.)
    call set_angles(solver, sundir)

    call solve_pprts(solver, lthermal=.not. lsolar, lsolar=lsolar, edirTOA=incSolar)
    call pprts_get_result(solver, fdn, fup, fdiv, fdir)

    call destroy_pprts(solver, lfinalizepetsc=.true.)
  end subroutine

  real(ireals) function dmean(comm, arr) result(m)
    integer(mpiint), intent(in) :: comm
    real(ireals), intent(in) :: arr(:, :)
    call imp_allreduce_mean(comm, arr, m)
  end function

  !> @brief solar, purely absorbing atmosphere but with a planck field present (regression test for thermal emission in solar runs)
  @test(npes=[2, 1])
  subroutine test_mcdmda_solar_absorbing_with_planck(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: theta0 = 60, dtau = 1, B = 100
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    real(ireals) :: mu, trgt, eps
    integer(mpiint) :: comm
    integer(iintegers) :: k

    comm = this%getMpiCommunicator()
    mu = cos(deg2rad(theta0))
    eps = 1e-4_ireals * incSolar

    call run_mcdmda(comm, .true., theta0, 0._ireals, dtau, 0._ireals, 0._ireals, B, fdir, fdn, fup, fdiv)

    do k = 1, Nlay + 1
      trgt = incSolar * mu * exp(-dtau * real(k - 1, ireals) / real(Nlay, ireals) / mu)
      @mpiassertEqual(trgt, dmean(comm, fdir(k, :, :)), eps, 'edir @ level '//toStr(k))
      @mpiassertEqual(0._ireals, maxval(abs(fdn(k, :, :))), eps, 'edn has to be zero without scattering @ level '//toStr(k))
      @mpiassertEqual(0._ireals, maxval(abs(fup(k, :, :))), eps, 'eup has to be zero without scattering @ level '//toStr(k))
    end do

    do k = 1, Nlay
      trgt = (dmean(comm, fdir(k, :, :)) - dmean(comm, fdir(k + 1, :, :))) / dz
      @mpiassertEqual(trgt, dmean(comm, fdiv(k, :, :)), eps / dz, 'absorption @ layer '//toStr(k))
    end do
  end subroutine

  !> @brief solar, conservative scattering, i.e. all energy that is not reflected to space is absorbed at the surface
  @test(npes=[2, 1])
  subroutine test_mcdmda_solar_conservative_scattering(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: theta0 = 30, dtau = 5, albedo = .3, g = .85, B = 100
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    real(ireals) :: mu, net_toa, net_srfc, eps, refl
    integer(mpiint) :: comm
    integer(iintegers) :: k

    comm = this%getMpiCommunicator()
    mu = cos(deg2rad(theta0))
    eps = 1e-4_ireals * incSolar

    call run_mcdmda(comm, .true., theta0, albedo, dtau, 1._ireals, g, B, fdir, fdn, fup, fdiv)

    @mpiassertEqual(incSolar * mu, dmean(comm, fdir(1, :, :)), eps, 'TOA edir')
    @mpiassertEqual(0._ireals, dmean(comm, fdn(1, :, :)), eps, 'TOA edn')

    net_toa = dmean(comm, fdir(1, :, :) + fdn(1, :, :) - fup(1, :, :))
    net_srfc = dmean(comm, fdir(Nlay + 1, :, :) + fdn(Nlay + 1, :, :) - fup(Nlay + 1, :, :))
    @mpiassertEqual(net_toa, net_srfc, eps, 'net flux TOA vs srfc (no absorption in the atmosphere)')

    @mpiassertEqual(albedo * dmean(comm, fdir(Nlay + 1, :, :) + fdn(Nlay + 1, :, :)), dmean(comm, fup(Nlay + 1, :, :)), eps, 'srfc eup')

    do k = 1, Nlay
      @mpiassertEqual(0._ireals, dmean(comm, fdiv(k, :, :)), eps / dz, 'absorption @ layer '//toStr(k))
    end do

    ! plausibility of the reflectance, i.e. we are not testing a trivial solution
    refl = dmean(comm, fup(1, :, :)) / (incSolar * mu)
    @mpiassertTrue(refl .gt. .1_ireals .and. refl .lt. .9_ireals, 'TOA reflectance implausible '//toStr(refl))
  end subroutine

  !> @brief thermal, isothermal purely absorbing atmosphere with surface at the same temperature, black surface
  !> @details every upward photon carries exactly the blackbody radiance, i.e. eup = pi * B everywhere
  @test(npes=[2, 1])
  subroutine test_mcdmda_thermal_isothermal(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: dtau = 2, B = 100
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    real(ireals) :: eps
    integer(mpiint) :: comm
    integer(iintegers) :: k

    comm = this%getMpiCommunicator()
    eps = 1e-4_ireals * pi * B

    call run_mcdmda(comm, .false., 0._ireals, 0._ireals, dtau, 0._ireals, 0._ireals, B, fdir, fdn, fup, fdiv)

    @mpiassertEqual(0._ireals, dmean(comm, fdn(1, :, :)), eps, 'TOA edn')
    do k = 1, Nlay + 1
      @mpiassertEqual(pi * B, minval(fup(k, :, :)), eps, 'min eup @ level '//toStr(k))
      @mpiassertEqual(pi * B, maxval(fup(k, :, :)), eps, 'max eup @ level '//toStr(k))
      @mpiassertTrue(minval(fdn(k, :, :)) .ge. -eps .and. maxval(fdn(k, :, :)) .le. pi * B + eps, 'edn out of bounds @ '//toStr(k))
    end do
    ! edn has to increase towards the surface
    @mpiassertTrue(dmean(comm, fdn(Nlay + 1, :, :)) .gt. dmean(comm, fdn(2, :, :)), 'edn does not increase downwards')
  end subroutine
end module
