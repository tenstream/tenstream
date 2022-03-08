module test_plexrt_rrtmg_lw_sw
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & ireals, iintegers, mpiint, default_str_len

  use m_helper_functions, only: chkerr, spherical_2_cartesian, meanval

  use m_examples_plex_rrtmg_fish, only: ex_plex_rrtmg_fish

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

  @test(npes=[1, 2])
  subroutine plexrt_rrtmg_sw(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .true.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .false., lsolar = .true.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 2, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = 1 ! TOA idx
    trgt = -sundir(3) * S0
    @assertEqual(trgt, edir(k, :), trgt * eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k, :), trgt * eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 346.0154 ! eup @ TOA determined with 2str
    trgt = 342.7360 ! eup @ TOA determined with rayli
    trgt = 343.4313 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    trgt = 4.07007942e-03 ! abso @ TOA determined with rayli
    trgt = 4.07351227e-03 ! abso @ TOA determined with plexrt
    @assertEqual(trgt, abso(k, :), trgt * eps, 'abso @ TOA')

    k = size(edir, dim=1) ! surface idx
    @assertEqual((edir(k, :) + edn(k, :)) * Ag_sol, eup(k, :), eps, 'eup should be (edir+edn)*Ag')

    trgt = 1050.883 ! edir @ srfc determined with 2str
    trgt = 1052.004 ! edir @ srfc determined with rayli
    trgt = 1051.049 ! edir @ srfc determined with plexrt
    @assertEqual(trgt, edir(k, :), trgt * eps, 'direct radiation @ surface')

    trgt = 72.29650 ! edn @ srfc determined with 2str
    trgt = 75.53327 ! edn @ srfc determined with rayli
    trgt = 73.24891 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse radiation @ surface')

    k = size(abso, dim=1) ! surface idx
    trgt = 3.69826481e-02 ! abso @ srfc determined with rayli
    trgt = 3.79099995e-02 ! abso @ srfc determined with plexrt
    @assertEqual(trgt, abso(k, :), trgt * eps, 'abso @ TOA')

  end subroutine

  @test(npes=[1, 2])
  subroutine plexrt_rrtmg_lw(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .true.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .true., lsolar = .false.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 2, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = 1 ! TOA idx
    @assertEqual(0._ireals, edn(k, :), tiny(edn), 'diffuse down radiation @ TOA should be zero')

    trgt = 251.9889  ! eup @ TOA determined with 2str
    trgt = 251.3705  ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    k = size(edn, dim=1) ! surface idx

    trgt = 292.0275  ! edn @ srfc determined with 2str
    trgt = 290.8192  ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse down radiation @ surface')

    trgt = 387.195 ! eup @ srfc determined with 2str
    trgt = 386.140 ! eup @ srfc determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse up radiation @ surface')
  end subroutine

  @test(npes=[1, 2])
  subroutine plexrt_rrtmg_fish_sw(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .false.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .false., lsolar = .true.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 3, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = 1 ! TOA idx
    trgt = -sundir(3) * S0
    @assertEqual(trgt, edir(k, :), trgt * eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k, :), trgt * eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 346.0154 ! eup @ TOA determined with 2str -- regular mesh
    trgt = 345.7406 ! eup @ TOA determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    k = size(edir, dim=1) ! surface idx
    @assertEqual((edir(k, :) + edn(k, :)) * Ag_sol, eup(k, :), eps, 'eup should be (edir+edn)*Ag')

    trgt = 1050.883 ! edir @ srfc determined with 2str -- regular mesh
    @assertEqual(trgt, edir(k, :), trgt * eps, 'direct radiation @ surface')

    trgt = 72.29650 ! edn @ srfc determined with 2str -- regular mesh
    trgt = 73.63184 ! edn @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse radiation @ surface')
  end subroutine

  @test(npes=[1, 2])
  subroutine plexrt_rrtmg_fish_lw(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .false.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .true., lsolar = .false.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 3, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = 1 ! TOA idx
    @assertEqual(0._ireals, edn(k, :), tiny(edn), 'diffuse down radiation @ TOA should be zero')

    trgt = 251.99    ! eup @ TOA determined with 2str -- regular mesh
    trgt = 251.42    ! eup @ TOA determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    k = size(edn, dim=1) ! surface idx

    trgt = 292.03    ! edn @ srfc determined with 2str -- regular mesh
    trgt = 291.70    ! edn @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse down radiation @ surface')

    trgt = 387.20    ! eup @ srfc determined with 2str -- regular mesh
    trgt = 386.184   ! eup @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse up radiation @ surface')
  end subroutine

  @test(npes=[2, 1])
  subroutine plexrt_rrtmg_sw_cld(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .true.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .false., lsolar = .true.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 2, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 1
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(iintegers) :: k
    integer(mpiint) :: comm

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = 1 ! TOA idx
    trgt = -sundir(3) * S0
    @assertEqual(trgt, edir(k, :), trgt * eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k, :), trgt * eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 655.0020 ! eup @ TOA determined with 2str
    trgt = 658.2601 ! eup @ TOA determined with rayli
    trgt = 613.3604 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    k = size(edir, dim=1, kind=iintegers) ! surface idx
    @assertEqual((edir(k, :) + edn(k, :)) * Ag_sol, eup(k, :), eps, 'eup should be (edir+edn)*Ag')

    trgt = 8.018153 ! edir @ srfc determined with 2str
    trgt = 7.824601 ! edir @ srfc determined with rayli
    trgt = 10.21104 ! edir @ srfc determined with plexrt
    @assertEqual(trgt, edir(k, :), trgt * eps, 'direct radiation @ surface')

    trgt = 537.3950 ! edn @ srfc determined with 2str
    trgt = 533.3429 ! edn @ srfc determined with rayli
    trgt = 592.9749 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse radiation @ surface')

    k = size(edir, dim=1, kind=iintegers) - Nz + 1 ! above cloud idx
    trgt = 1060.844 ! edir determined with 2str
    trgt = 1060.979 ! edir determined with rayli
    trgt = 1060.844 ! edir determined with plexrt
    @assertEqual(trgt, edir(k, :), trgt * eps, 'direct radiation @ cld top')

    trgt = 94.81436 ! edn determined with 2str
    trgt = 104.6741 ! edn determined with rayli
    trgt = 93.82275 ! edn determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'edn radiation @ cld top')

    trgt = 688.9160 ! eup determined with 2str
    trgt = 702.0040 ! eup determined with rayli
    trgt = 640.7009 ! eup determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'eup radiation @ cld top')

    trgt = 0.838109791 ! abso determined with 2str
    trgt = 0.837467730 ! abso determined with rayli
    trgt = 0.875359833 ! abso determined with plexrt
    @assertEqual(trgt, abso(k, :), trgt * eps, 'abso @ cld layer')

    k = size(edir, dim=1, kind=iintegers) - Nz + 2 ! below cloud idx
    trgt = 8.061157 ! edir determined with 2str
    trgt = 7.866495 ! edir determined with rayli
    trgt = 10.26560 ! edir determined with plexrt
    @assertEqual(trgt, edir(k, :), trgt * eps, 'direct radiation @ cld bot')

    trgt = 540.7412 ! edn determined with 2str
    trgt = 536.6765 ! edn determined with rayli
    trgt = 597.6715 ! edn determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'edn radiation @ cld bot')

    trgt = 163.6163 ! eup determined with 2str
    trgt = 162.3884 ! eup determined with rayli
    trgt = 181.1754 ! eup determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'eup radiation @ cld bot')
  end subroutine

  @test(npes=[2, 1])
  subroutine plexrt_rrtmg_lw_cld(this)
    class(MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh = .true.
    logical, parameter :: lverbose = .true.
    logical, parameter :: lthermal = .true., lsolar = .false.
    character(len=*), parameter :: atm_filename = "atm.dat"
    integer(iintegers), parameter :: Nx = 4, Ny = 2, Nz = 4
    real(ireals), parameter :: dx = 100, dz = 100
    real(ireals), parameter :: Ag_th = .05, Ag_sol = .3
    real(ireals), parameter :: phi = 0, theta = 0
    real(ireals), parameter :: lwc = 1
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals), parameter :: eps = 1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(iintegers) :: k
    integer(mpiint) :: comm

    comm = this%getMpiCommunicator()

    sundir = spherical_2_cartesian(phi, theta)

    call ex_plex_rrtmg_fish(&
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    k = size(edn, dim=1, kind=iintegers) - Nz + 1_iintegers ! above cloud idx
    print *, 'abso', abso(k, :)
    print *, 'edn ', edn(k, :)
    print *, 'eup ', eup(k, :)

    k = 1 ! TOA idx
    @assertEqual(0._ireals, edn(k, :), 'diffuse down radiation @ TOA should be zero')

    trgt = 251.96   ! eup @ TOA determined with 2str
    trgt = 251.46   ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse upward radiation @ TOA')

    k = size(edn, dim=1, kind=iintegers) ! surface idx

    trgt = 386.82   ! edn @ srfc determined with 2str
    trgt = 384.13   ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'diffuse dn radiation @ surface')

    trgt = 390.95   ! eup @ srfc determined with 2str
    trgt = 390.81   ! eup @ srfc determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'diffuse up radiation @ surface')

    k = size(edn, dim=1, kind=iintegers) - Nz + 1_iintegers ! above cloud idx
    trgt = 280.63   ! edn determined with 2str
    trgt = 280.63   ! edn determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'edn radiation @ cld top')

    trgt = 380.92   ! eup determined with 2str
    trgt = 379.19   ! eup determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'eup radiation @ cld top')

    trgt = -0.98677    ! abso determined with 2str
    trgt = -0.95925    ! abso determined with plexrt
    @assertEqual(trgt, meanval(abso(k, :)), abs(trgt) * eps, 'abso @ cld layer')
    @assertEqual(trgt, abso(k, :), abs(trgt) * sqrt(eps), 'abso @ cld layer should also be horizontally homogeneous')

    k = size(edn, dim=1, kind=iintegers) - Nz + 2_iintegers ! below cloud idx
    trgt = 383.90   ! edn determined with 2str
    trgt = 381.74   ! edn determined with plexrt
    @assertEqual(trgt, edn(k, :), trgt * eps, 'edn radiation @ cld bot')

    trgt = 388.176  ! eup determined with 2str
    trgt = 386.954  ! eup determined with plexrt
    @assertEqual(trgt, eup(k, :), trgt * eps, 'eup radiation @ cld bot')
  end subroutine
end module
