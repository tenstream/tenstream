module test_plexrt_fish

  use m_examples_plex_ex_fish, only: plex_ex_fish

  use m_tenstream_options, only: read_commandline_options

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi

  use m_helper_functions, only: &
    & CHKERR, &
    & imp_allreduce_mean, &
    & spherical_2_cartesian

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm = this%getMpiCommunicator()
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm = this%getMpiCommunicator()
    call finalize_mpi(comm, lfinalize_mpi=.false., lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes=[4, 2, 1])
  subroutine plexrt_regular_sw_sza0(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 4, Ny = 5, Nz = 3
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0., 0.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz, :) + edn(Nz, :)) * Ag, eup(Nz, :), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine plexrt_regular_sw_sza20(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 6, Ny = 3, Nz = 3
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 2e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90., 20.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz, :) + edn(Nz, :)) * Ag, eup(Nz, :), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine plexrt_regular_sw_sza40_no_scatter_no_Ag(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 6, Ny = 3, Nz = 4
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 5e-2

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90., 40.)
    w0 = 0
    Ag = 0

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(zero, eup, eps, 'no scatter, no albedo, eup should be 0')
    @assertEqual(zero, edn, eps, 'no scatter, no albedo, edn should be 0')
  end subroutine

  @test(npes=[2, 1])
  subroutine plexrt_regular_lw_emis_dtau_high(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 2, Ny = 3, Nz = 3
    logical, parameter :: lthermal = .true., lsolar = .false., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1000 * (Nz - 1), g = 0, Bplck = 1
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    integer(iintegers) :: k
    real(ireals) :: m

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0., 0.)
    w0 = .0
    Ag = .0

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @mpiassertEqual(0._ireals, edir, 'thermal computation, edir should be zero')
    @mpiassertEqual(0._ireals, edn(1, :), 'edn TOA should be zero')
    @mpiassertEqual(Bplck*pi, edn(2:Nz,:), Bplck*pi*eps, 'for high optical absorption thickness Edn emissivities should be plain Planck')
    @mpiassertEqual(Bplck * pi, eup, Bplck * pi * eps, 'for high optical absorption thickness Eup emissivities should be plain Planck')
    @mpiassertEqual(-Bplck*pi, abso(1,:) , Bplck*pi*eps, 'for high optical absorption thickness top most layer should emit plain Planck')
    @mpiassertEqual(0._ireals, abso(2:Nz - 1, :), eps, 'for high optical absorption thickness absorption in atm should be saturated')
    do k = 2, size(abso, 1)
      call imp_allreduce_mean(comm, abso(k, :), m)
      @mpiassertEqual(0._ireals, m, eps, 'for high optical absorption thickness absorption in atm should be zero')
    end do
  end subroutine

  @test(npes=[2, 1])
  subroutine plexrt_regular_lw_emis_dtau_low(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 2, Ny = 3, Nz = 10
    logical, parameter :: lthermal = .true., lsolar = .false., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 0 * (Nz - 1), g = 0, Bplck = 100
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .true.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    integer(iintegers) :: k
    real(ireals) :: m

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0., 0.)
    w0 = .0
    Ag = .5

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @mpiassertEqual(0._ireals, edir, 'thermal computation, edir should be zero')
    @mpiassertEqual(0._ireals, edn, eps, 'for low optical absorption thickness Edn emissivities should be zero')
    do k = 1, size(eup, 1)
      call imp_allreduce_mean(comm, eup(k, :), m)
      @mpiassertEqual((1._ireals - Ag)*Bplck*pi, m, Bplck*pi*eps, 'for low optical absorption thickness Eup emissivities should be 1-Albedo * Srfc_Planck')

      call imp_allreduce_mean(comm, abso(k, :), m)
      @mpiassertEqual(0._ireals, m, eps, 'for low optical absorption thickness absorption in atm should be zero')
    end do
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine plexrt_fish_sw_sza0(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 4, Ny = 5, Nz = 3
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .false.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0., 0.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz, :) + edn(Nz, :)) * Ag, eup(Nz, :), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine plexrt_fish_sw_sza20(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 6, Ny = 3, Nz = 3
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 2e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .false.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90., 20.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz, :) + edn(Nz, :)) * Ag, eup(Nz, :), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine plexrt_fish_sw_sza40_no_scatter_no_Ag(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx = 6, Ny = 3, Nz = 3
    logical, parameter :: lthermal = .false., lsolar = .true., lverbose = .true.
    real(ireals), parameter :: dz = 1, dtau = 1, g = 0, Bplck = -1
    real(ireals), parameter :: eps = 5e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh = .false.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    comm = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90., 40.)
    w0 = 0
    Ag = 0

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, g, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1, :), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(zero, eup, eps, 'no scatter, no albedo, eup should be 0')
    @assertEqual(zero, edn, eps, 'no scatter, no albedo, edn should be 0')
  end subroutine

end module
