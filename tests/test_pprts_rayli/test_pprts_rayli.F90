module test_pprts_rayli
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & default_str_len
  use pfunit_mod
  use m_helper_functions, only: deg2rad, imp_allreduce_mean, CHKERR
  use m_pprts_external_solvers, only: destroy_rayli_info
  use m_examples_pprts_ex1, only: pprts_ex1
  use m_example_pprts_rrtm_lw_sw, only: ex_pprts_rrtm_lw_sw

  implicit none

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call destroy_rayli_info()
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes=[4, 1])
  subroutine test_pprts_rayli_clear_sky(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 1
    real(ireals), parameter :: dx = 1, dy = dx, dz = 1
    real(ireals), parameter :: phi0 = 180, theta0 = 60
    real(ireals), parameter :: albedo = 0
    real(ireals), parameter :: incSolar = 1
    real(ireals), parameter :: Bplck = 0, Bplck_srfc = 0
    logical, parameter :: lthermal = .false., lsolar = .true.

    real(ireals), parameter :: dtau_clearsky = 1, w0_clearsky = 0, g_clearsky = 0
    integer(iintegers), parameter :: cld_layer_idx(2) = [0, -1]
    real(ireals), parameter :: dtau_cloud = 1, w0_cloud = 0, g_cloud = 0
    logical, parameter :: lverbose = .true.

    real(ireals) :: trgt
    real(ireals) :: val
    real(ireals), parameter :: eps = sqrt(epsilon(eps))

    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    integer(iintegers) :: i

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call pprts_ex1( &
        & comm, &
        & lthermal, &
        & lsolar, &
        & Nx, Ny, Nlay, &
        & dx, dy, &
        & phi0, theta0, &
        & albedo, dz, &
        & incSolar, &
        & Bplck, &
        & Bplck_srfc, &
        & dtau_clearsky, w0_clearsky, g_clearsky, &
        & cld_layer_idx, &
        & dtau_cloud, w0_cloud, g_cloud, &
        & fdir, fdn, fup, fdiv, &
        & lverbose)

    trgt = incSolar * cos(deg2rad(theta0))
    call imp_allreduce_mean(comm, fdir(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA edir')

    trgt = incSolar * cos(deg2rad(theta0)) * exp(-dtau_clearsky / cos(deg2rad(theta0)))
    call imp_allreduce_mean(comm, fdir(Nlay + 1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC edir')

    ! check that flx divergence is same as absorption output
    do i = 1, size(fdiv, 1)
      call imp_allreduce_mean(comm, &
        & ((fdir(i, :, :) + fdn(i, :, :) - fup(i, :, :)) - (fdir(i + 1, :, :) + fdn(i + 1, :, :) - fup(i + 1, :, :))) / dz, &
        & trgt)
      call imp_allreduce_mean(comm, fdiv(i, :, :), val)
      @mpiassertEqual(trgt, val, eps, 'flx div consistency')

      trgt = incSolar * cos(deg2rad(theta0)) * (1._ireals - exp(-dtau_clearsky / cos(deg2rad(theta0))))
      @mpiassertEqual(trgt, val, eps, 'absorption wrong')
    end do

  end subroutine

  @test(npes=[4, 1])
  subroutine test_pprts_rayli_single_cld_lay(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 4
    real(ireals), parameter :: dx = 1, dy = 1, dz = 5
    real(ireals), parameter :: phi0 = 180, theta0 = 0
    real(ireals), parameter :: albedo = 0.1
    real(ireals), parameter :: incSolar = 1
    real(ireals), parameter :: Bplck = 0, Bplck_srfc = 0
    logical, parameter :: lthermal = .false., lsolar = .true.

    real(ireals), parameter :: dtau_clearsky = 1, w0_clearsky = 0, g_clearsky = 0
    integer(iintegers), parameter :: cld_layer_idx(2) = [2, 3]
    real(ireals), parameter :: dtau_cloud = 1, w0_cloud = 0.2, g_cloud = 0.2
    logical, parameter :: lverbose = .true.

    real(ireals), parameter :: eps = 1e-2
    real(ireals) :: trgt, val
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv
    integer(iintegers) :: i

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call pprts_ex1( &
        & comm, &
        & lthermal, &
        & lsolar, &
        & Nx, Ny, Nlay, &
        & dx, dy, &
        & phi0, theta0, &
        & albedo, dz, &
        & incSolar, &
        & Bplck, &
        & Bplck_srfc, &
        & dtau_clearsky, w0_clearsky, g_clearsky, &
        & cld_layer_idx, &
        & dtau_cloud, w0_cloud, g_cloud, &
        & fdir, fdn, fup, fdiv, &
        & lverbose)

    trgt = incSolar * cos(deg2rad(theta0))
    call imp_allreduce_mean(comm, fdir(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA edir')

    trgt = incSolar * cos(deg2rad(theta0)) * exp(-(dtau_clearsky + dtau_cloud) / cos(deg2rad(theta0)))
    call imp_allreduce_mean(comm, fdir(Nlay + 1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC edir')

    ! check that flx divergence is same as absorption output
    do i = 1, size(fdiv, 1)
      call imp_allreduce_mean(comm, &
        & ((fdir(i, :, :) + fdn(i, :, :) - fup(i, :, :)) - (fdir(i + 1, :, :) + fdn(i + 1, :, :) - fup(i + 1, :, :))) / dz, &
        & trgt)
      call imp_allreduce_mean(comm, fdiv(i, :, :), val)
      @mpiassertEqual(trgt, val, eps, 'flx div consistency')
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_pprts_rayli_rrtmg_solar(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 4
    real(ireals), parameter :: dx = 1e6, dy = 1e6
    real(ireals), parameter :: phi0 = 180, theta0 = 0
    real(ireals), parameter :: albedo_th = 0.1
    real(ireals), parameter :: albedo_sol = 0.25
    logical, parameter :: lthermal = .false., lsolar = .true.
    character(len=default_str_len), parameter :: atm_filename = 'atm.dat'

    integer(iintegers) :: N
    real(ireals), parameter :: incSolar = 1368.222
    real(ireals), parameter :: eps = 10
    real(ireals) :: trgt, val
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, "-rayli_photons 1e2", ierr); call CHKERR(ierr)

    call ex_pprts_rrtm_lw_sw( &
      & comm, Nx, Ny, Nlay, &
      & dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & fdir, fdn, fup, fdiv, &
      & vlwc=1e-2_ireals)

    N = ubound(fdn, 1)

    trgt = incSolar * cos(deg2rad(theta0))
    call imp_allreduce_mean(comm, fdir(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA edir')

    trgt = 0
    call imp_allreduce_mean(comm, fdn(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA edn')

    trgt = 340.2
    call imp_allreduce_mean(comm, fup(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA eup')

    trgt = 612.1
    call imp_allreduce_mean(comm, fdir(N, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC edir')

    trgt = 423.3
    call imp_allreduce_mean(comm, fdn(N, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC edn')

    trgt = (612.1 + 423.3) * albedo_sol
    call imp_allreduce_mean(comm, fup(N, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC eup')
  end subroutine

  @test(npes=[1])
  subroutine test_pprts_rayli_rrtmg_thermal(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 4
    real(ireals), parameter :: dx = 1e6, dy = 1e6
    real(ireals), parameter :: phi0 = 180, theta0 = 0
    real(ireals), parameter :: albedo_th = 0.1
    real(ireals), parameter :: albedo_sol = 0.25
    logical, parameter :: lthermal = .true., lsolar = .false.
    character(len=default_str_len), parameter :: atm_filename = 'atm.dat'

    integer(iintegers) :: N
    real(ireals), parameter :: eps = 10
    real(ireals) :: trgt, val
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, "-rayli_photons 1e2", ierr); call CHKERR(ierr)

    call ex_pprts_rrtm_lw_sw( &
      & comm, Nx, Ny, Nlay, &
      & dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & fdir, fdn, fup, fdiv, &
      & vlwc=1e-2_ireals)

    N = ubound(fdn, 1)

    trgt = 0
    call imp_allreduce_mean(comm, fdn(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA edn')

    trgt = 247.2
    call imp_allreduce_mean(comm, fup(1, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'TOA eup')

    trgt = 342.4
    call imp_allreduce_mean(comm, fdn(N, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC edn')

    trgt = 385.3
    call imp_allreduce_mean(comm, fup(N, :, :), val)
    @mpiassertEqual(trgt, val, eps, 'SRFC eup')
  end subroutine
end module
