module test_plexrt_rayli

  use m_examples_plex_ex_fish, only: plex_ex_fish

  use m_tenstream_options, only: read_commandline_options

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & meanval, &
    & spherical_2_cartesian

  use pfunit_mod

  implicit none

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
    call read_commandline_options(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes=[1])
  subroutine plexrt_regular_sw_sza0_abso_only(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 2, Ny = 2, Nz = 2
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0., 0.)
    w0 = .0
    Ag = .1

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), meanval(edir(1, :)), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-dtau/abs(sundir(3))), meanval(edir(Nz,:)), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(meanval(edir(Nz, :) + edn(Nz, :)) * Ag, meanval(eup(Nz, :)), eps, 'eup should be (edir+edn)*Ag')
    @assertEqual(meanval(eup(Nz, :) * exp(-dtau)), meanval(eup(1, :)), eps, 'eup TOA')
  end subroutine

  @test(npes=[1])
  subroutine plexrt_regular_sw_sza0(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 2, Ny = 2, Nz = 2
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(270., 40.)
    w0 = 1
    Ag = .1

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), meanval(edir(1, :)), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-dtau/abs(sundir(3))), meanval(edir(Nz,:)), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(meanval(edir(Nz, :) + edn(Nz, :)) * Ag, meanval(eup(Nz, :)), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

end module
