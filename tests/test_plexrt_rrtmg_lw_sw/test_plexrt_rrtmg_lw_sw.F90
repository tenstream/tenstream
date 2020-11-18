module test_plexrt_rrtmg_lw_sw
  use m_tenstream_options, only : read_commandline_options
  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & ireals, iintegers, mpiint, default_str_len

  use m_helper_functions, only: chkerr, spherical_2_cartesian

  use m_examples_plex_rrtmg_fish, only: ex_plex_rrtmg_fish

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm     = this%getMpiCommunicator()
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm     = this%getMpiCommunicator()
    call finalize_mpi(comm, lfinalize_mpi=.False., lfinalize_petsc=.True.)
  end subroutine teardown

  !@test(npes =[1,2])
  subroutine plexrt_rrtmg_sw(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.True.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.False., lsolar=.True.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=2, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm     = this%getMpiCommunicator()

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
    trgt=-sundir(3)*S0
    @assertEqual(trgt, edir(k,:), trgt*eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k,:), trgt*eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 346.0154 ! eup @ TOA determined with 2str
    trgt = 345.7406 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edir,dim=1) ! surface idx
    @assertEqual((edir(k,:)+edn(k,:))*Ag_sol, eup(k,:), eps, 'eup should be (edir+edn)*Ag')

    trgt = 1050.883 ! edir @ srfc determined with 2str
    @assertEqual(trgt, edir(k,:), trgt*eps, 'direct radiation @ surface')

    trgt = 72.29650 ! edn @ srfc determined with 2str
    trgt = 72.28502 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse radiation @ surface')
  end subroutine

  !@test(npes =[1,2])
  subroutine plexrt_rrtmg_lw(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.True.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.True., lsolar=.False.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=2, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm     = this%getMpiCommunicator()

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
    @assertEqual(0._ireals, edn(k,:), tiny(edn), 'diffuse down radiation @ TOA should be zero')

    trgt = 254.99942 ! eup @ TOA determined with 2str
    trgt = 254.85485 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edn,dim=1) ! surface idx

    trgt = 288.03708 ! edn @ srfc determined with 2str
    trgt = 288.51068 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse down radiation @ surface')

    trgt = 386.99545 ! eup @ srfc determined with 2str
    trgt = 386.02422 ! eup @ srfc determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse up radiation @ surface')
  end subroutine

  !@test(npes =[1,2])
  subroutine plexrt_rrtmg_fish_sw(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.False.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.False., lsolar=.True.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=3, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm     = this%getMpiCommunicator()

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
    trgt=-sundir(3)*S0
    @assertEqual(trgt, edir(k,:), trgt*eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k,:), trgt*eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 346.0154 ! eup @ TOA determined with 2str -- regular mesh
    trgt = 345.7406 ! eup @ TOA determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edir,dim=1) ! surface idx
    @assertEqual((edir(k,:)+edn(k,:))*Ag_sol, eup(k,:), eps, 'eup should be (edir+edn)*Ag')

    trgt = 1050.883 ! edir @ srfc determined with 2str -- regular mesh
    @assertEqual(trgt, edir(k,:), trgt*eps, 'direct radiation @ surface')

    trgt = 72.29650 ! edn @ srfc determined with 2str -- regular mesh
    trgt = 72.28502 ! edn @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse radiation @ surface')
  end subroutine

  !@test(npes =[1,2])
  subroutine plexrt_rrtmg_fish_lw(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.False.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.True., lsolar=.False.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=3, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 0
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(mpiint) :: comm, k

    comm     = this%getMpiCommunicator()

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
    @assertEqual(0._ireals, edn(k,:), tiny(edn), 'diffuse down radiation @ TOA should be zero')

    trgt = 254.99942 ! eup @ TOA determined with 2str -- regular mesh
    trgt = 254.85485 ! eup @ TOA determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edn,dim=1) ! surface idx

    trgt = 288.03708 ! edn @ srfc determined with 2str -- regular mesh
    trgt = 288.51068 ! edn @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse down radiation @ surface')

    trgt = 386.99545 ! eup @ srfc determined with 2str -- regular mesh
    trgt = 386.02422 ! eup @ srfc determined with plexrt -- regular mesh
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse up radiation @ surface')
  end subroutine

  !@test(npes =[1,2])
  subroutine plexrt_rrtmg_sw_cld(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.True.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.False., lsolar=.True.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=2, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 1
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals), parameter :: S0 = 1368.22229
    real(ireals) :: sundir(3), trgt
    integer(iintegers) :: k
    integer(mpiint) :: comm

    comm     = this%getMpiCommunicator()

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
    trgt=-sundir(3)*S0
    @assertEqual(trgt, edir(k,:), trgt*eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(0._ireals, edn(k,:), trgt*eps, 'diffuse down radiation @ TOA should be zero')

    trgt = 655.0020 ! eup @ TOA determined with 2str
    trgt = 614.6086 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edir,dim=1,kind=iintegers) ! surface idx
    @assertEqual((edir(k,:)+edn(k,:))*Ag_sol, eup(k,:), eps, 'eup should be (edir+edn)*Ag')

    trgt = 8.018153 ! edir @ srfc determined with 2str
    trgt = 10.21104 ! edir @ srfc determined with plexrt
    @assertEqual(trgt, edir(k,:), trgt*eps, 'direct radiation @ surface')

    trgt = 537.3950 ! edn @ srfc determined with 2str
    trgt = 595.8433 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse radiation @ surface')

    k = size(edir,dim=1,kind=iintegers) - Nz +1 ! above cloud idx
    trgt = 1060.844 ! edir determined with 2str
    trgt = 1060.844 ! edir determined with plexrt
    @assertEqual(trgt, edir(k,:), trgt*eps, 'direct radiation @ cld top')

    trgt = 94.81436 ! edn determined with 2str
    trgt = 91.66837 ! edn determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'edn radiation @ cld top')

    trgt = 688.9160 ! eup determined with 2str
    trgt = 642.0978 ! eup determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'eup radiation @ cld top')

    trgt = 0.838109791 ! abso determined with 2str
    trgt = 0.840916336 ! abso determined with plexrt
    @assertEqual(trgt, abso(k,:), trgt*eps, 'abso @ cld layer')

    k = size(edir,dim=1,kind=iintegers) - Nz +2 ! below cloud idx
    trgt = 8.061157 ! edir determined with 2str
    trgt = 10.26560 ! edir determined with plexrt
    @assertEqual(trgt, edir(k,:), trgt*eps, 'direct radiation @ cld bot')

    trgt = 540.7412 ! edn determined with 2str
    trgt = 600.2974 ! edn determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'edn radiation @ cld bot')

    trgt = 163.6163 ! eup determined with 2str
    trgt = 181.9779 ! eup determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'eup radiation @ cld bot')
  end subroutine

  @test(npes =[1])
  subroutine plexrt_rrtmg_lw_cld(this)
  class (MpiTestMethod), intent(inout) :: this
    logical, parameter :: lregular_mesh=.True.
    logical, parameter :: lverbose=.True.
    logical, parameter :: lthermal=.True., lsolar=.False.
    character(len=*), parameter :: atm_filename="atm.dat"
    integer(iintegers), parameter :: Nx=4, Ny=2, Nz=4
    real(ireals), parameter :: dx=100, dz=100
    real(ireals), parameter :: Ag_th=.05, Ag_sol=.3
    real(ireals), parameter :: phi=0, theta=0
    real(ireals), parameter :: lwc = 1
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    real(ireals), parameter :: eps=1e-2 ! better than 1 percent
    real(ireals) :: sundir(3), trgt
    integer(iintegers) :: k
    integer(mpiint) :: comm

    comm     = this%getMpiCommunicator()

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

    k = size(edn,dim=1,kind=iintegers) - Nz + 1_iintegers ! above cloud idx
    print *,'abso', abso(k,:)

    k = 1 ! TOA idx
    @assertEqual(0._ireals, edn(k,:), 'diffuse down radiation @ TOA should be zero')

    trgt = 256.2525 ! eup @ TOA determined with 2str
    trgt = 255.7656 ! eup @ TOA determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse upward radiation @ TOA')

    k = size(edn,dim=1,kind=iintegers) ! surface idx

    trgt = 386.6776 ! edn @ srfc determined with 2str
    trgt = 384.7783 ! edn @ srfc determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'diffuse dn radiation @ surface')

    trgt = 390.9326 ! eup @ srfc determined with 2str
    trgt = 390.8375 ! eup @ srfc determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'diffuse up radiation @ surface')

    k = size(edn,dim=1,kind=iintegers) - Nz + 1_iintegers ! above cloud idx
    trgt = 275.4420 ! edn determined with 2str
    trgt = 272.3709 ! edn determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'edn radiation @ cld top')

    trgt = 380.9236 ! eup determined with 2str
    trgt = 382.0963 ! eup determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'eup radiation @ cld top')

    trgt = -1.03882027 ! abso determined with 2str
    trgt = -1.05032933 ! abso determined with plexrt
    @assertEqual(trgt, abso(k,:), trgt*eps, 'abso @ cld layer')

    k = size(edn,dim=1,kind=iintegers) - Nz + 2_iintegers ! below cloud idx
    trgt = 383.8705 ! edn determined with 2str
    trgt = 381.3875 ! edn determined with plexrt
    @assertEqual(trgt, edn(k,:), trgt*eps, 'edn radiation @ cld bot')

    trgt = 388.2633 ! eup determined with 2str
    trgt = 388.9045 ! eup determined with plexrt
    @assertEqual(trgt, eup(k,:), trgt*eps, 'eup radiation @ cld bot')
  end subroutine
end module
