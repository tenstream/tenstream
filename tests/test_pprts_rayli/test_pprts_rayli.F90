module test_pprts_rayli

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint
  use pfunit_mod
  use m_helper_functions, only: deg2rad
  use m_examples_pprts_ex1, only: pprts_ex1

  implicit none

contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine teardown

  @test(npes = [1])
  subroutine test_pprts_rayli_check_direct_radiation(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: Nx = 3, Ny = 3, Nlay = 10
    real(ireals), parameter :: dx = 100, dy = 100, dz = 500
    real(ireals), parameter :: phi0 = 180, theta0 = 0
    real(ireals), parameter :: albedo = 0.1
    real(ireals), parameter :: incSolar = 1

    real(ireals), parameter :: dtau_clearsky = 1, w0_clearsky = 0, g_clearsky = 0
    integer(iintegers), parameter :: cld_layer_idx(2) = [0,-1]
    real(ireals), parameter :: dtau_cloud = 1, w0_cloud = 0, g_cloud = 0
    logical, parameter :: lverbose=.True.

    real(ireals) :: trgt
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call pprts_ex1( &
        & comm, &
        & Nx, Ny, Nlay, &
        & dx, dy, &
        & phi0, theta0, &
        & albedo, dz, &
        & incSolar, &
        & dtau_clearsky, w0_clearsky, g_clearsky, &
        & cld_layer_idx, &
        & dtau_cloud, w0_cloud, g_cloud, &
        & fdir,fdn,fup,fdiv, &
        & lverbose )

    trgt = incSolar * cos(deg2rad(theta0))
    @assertEqual(trgt, fdir(1,:,:), trgt * 0.03, 'TOA edir')

    trgt = incSolar * cos(deg2rad(theta0)) * exp(-dtau_clearsky/cos(deg2rad(theta0)))
    @assertEqual(trgt, fdir(Nlay + 1,:,:), trgt * 0.03, 'Should be equal')

  end subroutine
end module
