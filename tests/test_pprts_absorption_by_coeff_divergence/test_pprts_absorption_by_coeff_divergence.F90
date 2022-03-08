module test_pprts_absorption_by_coeff_divergence
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & pi

  use pfunit_mod
  use m_helper_functions, only: deg2rad, imp_allreduce_mean, CHKERR
  use m_pprts_external_solvers, only: destroy_rayli_info
  use m_examples_pprts_ex1, only: pprts_ex1

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
  subroutine test_pprts_thermal(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 4
    real(ireals), parameter :: dx = 1, dy = 1, dz = 1
    real(ireals), parameter :: phi0 = 180, theta0 = 0
    real(ireals), parameter :: albedo = 0.1
    real(ireals), parameter :: incSolar = 1
    real(ireals), parameter :: Bplck = 100 / pi, Bplck_srfc = 100
    logical, parameter :: lthermal = .true., lsolar = .false.

    real(ireals), parameter :: dtau_clearsky = 1, w0_clearsky = 0, g_clearsky = 0
    integer(iintegers), parameter :: cld_layer_idx(2) = [2, 3]
    real(ireals), parameter :: dtau_cloud = 1, w0_cloud = 0.2, g_cloud = 0.2
    logical, parameter :: lverbose = .true.

    real(ireals), parameter :: eps = 1e-2
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv, fdiv_coeff
    integer(mpiint) :: ierr

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, '-absorption_by_coeff_divergence no', ierr); call CHKERR(ierr)
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

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, '-absorption_by_coeff_divergence yes', ierr); call CHKERR(ierr)
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
      & fdir, fdn, fup, fdiv_coeff, &
      & lverbose)

    print *, 'fdiv by divergence', fdiv
    print *, 'fdiv by coeff_div ', fdiv_coeff

    @mpiassertEqual(fdiv, fdiv_coeff, maxval(fdiv) * eps, 'absorption should be same for both methods')

  end subroutine
  @test(npes=[4, 1])
  subroutine test_pprts_solar(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 4
    real(ireals), parameter :: dx = 1, dy = 1, dz = 1
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
    real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv, fdiv_coeff
    integer(mpiint) :: ierr

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, '-absorption_by_coeff_divergence no', ierr); call CHKERR(ierr)
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

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, '-absorption_by_coeff_divergence yes', ierr); call CHKERR(ierr)
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
      & fdir, fdn, fup, fdiv_coeff, &
      & lverbose)

    print *, 'fdiv by divergence', fdiv
    print *, 'fdiv by coeff_div ', fdiv_coeff

    @mpiassertEqual(fdiv, fdiv_coeff, maxval(fdiv) * eps, 'absorption should be same for both methods')

  end subroutine

end module
