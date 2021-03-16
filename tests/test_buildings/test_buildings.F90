module test_buildings
#include "petsc/finclude/petsc.h"
  use petsc
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only : &
    init_mpi_data_parameters,   &
    finalize_mpi,               &
    iintegers, ireals, mpiint,  &
    zero, one, pi, default_str_len

  use m_helper_functions, only : &
    & CHKERR, &
    & toStr, cstr, &
    & deg2rad, &
    & spherical_2_cartesian, &
    & rotate_angle_z, &
    & meanval, &
    & is_inrange, &
    & ind_1d_to_nd

  use m_pprts_base, only: &
    & t_solver, &
    & destroy_pprts, &
    & allocate_pprts_solver_from_commandline

  use m_pprts, only : init_pprts, &
    & set_optical_properties, solve_pprts, &
    & pprts_get_result, set_angles, &
    & gather_all_toZero

  use m_buildings, only: &
    & t_pprts_buildings, &
    & check_buildings_consistency, &
    & faceidx_by_cell_plus_offset, &
    & init_buildings,    &
    & clone_buildings,   &
    & destroy_buildings, &
    & PPRTS_TOP_FACE,    &
    & PPRTS_BOT_FACE,    &
    & PPRTS_LEFT_FACE,   &
    & PPRTS_RIGHT_FACE,  &
    & PPRTS_REAR_FACE,   &
    & PPRTS_FRONT_FACE

  use m_examples_pprts_buildings, only: ex_pprts_buildings
  use m_examples_pprts_rrtm_buildings, only: ex_pprts_rrtm_buildings

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.False., &
      & lfinalize_petsc=.True.)
  end subroutine teardown

  @test(npes =[1])
  subroutine test_init_buildings(this)
    class (MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    integer(iintegers), parameter :: Nfaces=6
    integer(iintegers), parameter :: da_sizes(4) = [ integer(iintegers) :: 6, 1, 1, 1 ]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings (B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated (B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    enddo
    call check_iface_data_consistency(B1)

    call destroy_buildings(B1, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B1))

    contains
      subroutine check_iface_data_consistency(B)
        type(t_pprts_buildings), intent(in) :: B
        integer(iintegers) :: i
        do i = 1, Nfaces
          @assertEqual(i, B%iface(i))
          @assertEqual(real(i, ireals), B%albedo(i))
        enddo
      end subroutine
  end subroutine

  @test(npes =[1])
  subroutine test_clone_buildings(this)
    class (MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    type(t_pprts_buildings), allocatable :: B2
    type(t_pprts_buildings), allocatable :: B3
    integer(iintegers), parameter :: Nfaces=6
    integer(iintegers), parameter :: da_sizes(4) = [ integer(iintegers) :: 6, 1, 1, 1 ]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings (B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated (B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    enddo
    call check_iface_data_consistency(B1)

    call clone_buildings(B1, B2, .True., ierr); call CHKERR(ierr)
    @assertEqual(2_iintegers, B1%iface_data%ref_count, 'origin ref_count did not increase')

    @assertTrue(allocated (B2%albedo))
    @assertTrue(associated(B2%iface))
    @assertTrue(associated(B2%iface_data%data))
    @assertEqual(Nfaces, size(B2%iface_data%data, kind=iintegers))
    @assertEqual(2_iintegers, B2%iface_data%ref_count)
    @assertEqual(B1%iface, B2%iface)
    @assertTrue(associated(B1%iface, B2%iface), 'cloned iface pointer dont point to same target')
    call check_iface_data_consistency(B2)

    call clone_buildings(B1, B3, .False., ierr); call CHKERR(ierr)
    @assertEqual(3_iintegers, B1%iface_data%ref_count, 'origin ref_count did not increase')

    @assertFalse(allocated(B3%albedo))
    @assertTrue(associated(B3%iface))
    @assertTrue(associated(B3%iface_data%data))
    @assertEqual(Nfaces, size(B3%iface_data%data, kind=iintegers))
    @assertEqual(3_iintegers, B3%iface_data%ref_count)
    @assertEqual(B1%iface, B3%iface)
    @assertTrue(associated(B1%iface, B3%iface), 'cloned iface pointer dont point to same target')
    call check_iface_data_consistency(B3)

    call destroy_buildings(B1, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B1))

    call check_iface_data_consistency(B2)
    call check_iface_data_consistency(B3)

    call destroy_buildings(B2, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B2))
    call destroy_buildings(B3, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B3))

    contains
      subroutine check_iface_data_consistency(B)
        type(t_pprts_buildings), intent(in) :: B
        integer(iintegers) :: i
        do i = 1, Nfaces
          @assertEqual(i, B%iface(i), "iface does not match")
          if(allocated(B%albedo)) then
            @assertEqual(real(i, ireals), B%albedo(i), "albedo does not match")
          endif
        enddo
      end subroutine
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_overhead_sun(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0
    box_planck = 0
    phi0       = 0
    theta0     = 0
    Ag         = 0
    dtau       = 0
    w0         = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      @assertEqual(0, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(S0/dz, gabso(glob_box_k, glob_box_i, glob_box_j), atol, 'box should absorb all incoming solar radiation')
    endif
    if(size(buildings%edir).gt.0) then
      @assertEqual(S0       , buildings%edir(1), atol)
      @assertEqual(0._ireals, buildings%edir(2:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_albedo(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0.5_ireals
    box_planck = 0
    phi0       = 0
    theta0     = 0
    Ag         = 1.0_ireals
    dtau       = 0
    w0         = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      @assertEqual( 0, gedir(Nlay+1      , glob_box_i, glob_box_j), atol, 'edir beneath building should be zero')
      @assertEqual( 0, gedir(glob_box_k+1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(S0, gedir(glob_box_k  , glob_box_i, glob_box_j), atol, 'edir at top of building should be S0')

      @assertEqual(S0*box_albedo, geup (glob_box_k  , glob_box_i, glob_box_j), atol, 'eup at top of building should be S0*B_Ag')

    endif
    if(size(buildings%edir).gt.0) then
      @assertEqual(S0       , buildings%edir(1), atol, 'edir at top of building should be S0')
      @assertEqual(0._ireals, buildings%edir(2:6), atol, 'if zenith angle is 0, all sides of the building except the top should be 0 edir')

      @assertEqual((buildings%edir+buildings%incoming)*box_albedo, buildings%outgoing, atol, 'total incoming times building albedo should give outgoing')
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_emission(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal   = .True.
    lsolar     = .False.
    box_albedo = 0._ireals
    box_planck = 1
    phi0       = 0
    theta0     = 180
    Ag         = 0
    dtau       = 0
    w0         = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      @assertEqual(box_planck*pi, geup(glob_box_k  , glob_box_i, glob_box_j), atol, 'eup at top of building should be emission')
      @assertEqual(box_planck*pi, gedn(glob_box_k+1, glob_box_i, glob_box_j), atol, 'edn at top of building should be emission')

    endif
    if(size(buildings%edir).gt.0) then
      @assertEqual(box_planck*pi, buildings%outgoing, atol, 'emission on buildings should be planck')
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_fortyfive_sun_azi0(this)
  class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0
    box_planck = 0
    phi0       = 0
    theta0     = 45
    Ag         = 0
    dtau       = 0
    w0         = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j  ), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j-1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j-1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j-2), atol, 'edir building shadow should go south and should be zero')
    endif
    if(size(buildings%edir).gt.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0*cos(deg2rad(theta0)) ! front face bc phi==0 is north)
      @assertEqual(trgt, buildings%edir(6), atol)

      @assertEqual(0._ireals              , buildings%edir(2:5), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_fortyfive_sun_azi180(this)
  class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0
    box_planck = 0
    phi0       = 180
    theta0     = 45
    Ag         = 0
    dtau       = 0
    w0         = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j  ), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j+1), atol, 'edir building shadow should go north and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j+1), atol, 'edir building shadow should go north and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j+2), atol, 'edir building shadow should go north and should be zero')
    endif
    if(size(buildings%edir).gt.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0*cos(deg2rad(theta0)) ! back face bc phi==180 is south)
      @assertEqual(trgt, buildings%edir(5), atol)

      @assertEqual(0._ireals, buildings%edir(2:4), atol)
      @assertEqual(0._ireals, buildings%edir(6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_fortyfive_sun_azi90(this)
  class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0
    box_planck = 0
    theta0     = 45
    Ag         = 0
    dtau       = 0
    w0         = 0
    phi0       = 90

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j  ), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i-1, glob_box_j), atol, 'edir building shadow should go west and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i-1, glob_box_j), atol, 'edir building shadow should go west and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i-2, glob_box_j), atol, 'edir building shadow should go west and should be zero')
    endif
    if(size(buildings%edir).gt.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0*cos(deg2rad(theta0)) ! right face bc phi==90 is east sun)
      @assertEqual(trgt, buildings%edir(4), atol)

      @assertEqual(0._ireals, buildings%edir(2:3), atol)
      @assertEqual(0._ireals, buildings%edir(5:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif
  end subroutine

  @test(npes =[4,2,1])
  subroutine test_buildings_example_fortyfive_sun_azi270(this)
  class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=3, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=3, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0
    box_planck = 0
    theta0     = 45
    Ag         = 0
    dtau       = 0
    w0         = 0
    phi0       = 270

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings                                 )

    if(myid.eq.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j  ), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i+1, glob_box_j), atol, 'edir building shadow should go east and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i+1, glob_box_j), atol, 'edir building shadow should go east and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i+2, glob_box_j), atol, 'edir building shadow should go east and should be zero')
    endif
    if(size(buildings%edir).gt.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0*cos(deg2rad(theta0)) ! leftt face bc phi==270 is west sun)
      @assertEqual(trgt, buildings%edir(3), atol)

      @assertEqual(0._ireals, buildings%edir(2), atol)
      @assertEqual(0._ireals, buildings%edir(4:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif
  end subroutine


  !@test(npes=[4,2,1])
  subroutine test_pprts_rrtmg_buildings_example_solar(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm, myid, ierr

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=10
    real(ireals), parameter :: dx=100, dy=100
    character(len=*), parameter :: atm_filename='afglus_100m.dat'
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: rtol=1e-5_ireals
    logical :: lsolar, lthermal
    real(ireals) :: buildings_albedo, buildings_temp
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag_solar, Ag_thermal
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso
    type(t_pprts_buildings), allocatable :: buildings_solar, buildings_thermal

    integer(iintegers) :: local_dims(6)
    integer(iintegers) :: ke, iface, idx(4)
    integer(iintegers) :: k, i, j

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lsolar   = .True.
    lthermal = .False.
    buildings_albedo = 0
    buildings_temp   = 0
    phi0   = 0
    theta0 = 0
    Ag_solar   = 0
    Ag_thermal = 0

    do i=1,2
      call ex_pprts_rrtm_buildings(           &
        & comm, lverbose,                     &
        & lthermal, lsolar,                   &
        & Nx, Ny, Nlay,                       &
        & buildings_albedo, buildings_temp,   &
        & dx, dy,                             &
        & atm_filename,                       &
        & phi0, theta0,                       &
        & Ag_solar, Ag_thermal,               &
        & gedir, gedn, geup, gabso,           &
        & buildings_solar, buildings_thermal, &
        & local_dims )
    enddo

    print *,myid, 'shape edir', shape(gedir)

    associate( Bs => buildings_solar, Bt => buildings_thermal )
      @assertFalse(allocated(Bt%edir), 'edir in thermal should never be allocated')
      @assertFalse(allocated(Bt%incoming), 'buildings_thermal should not have results allocated because lthermal=.F. ')
      @assertFalse(allocated(Bt%outgoing), 'buildings_thermal should not have results allocated because lthermal=.F. ')

      ke = size(gabso,dim=1)

      do iface = 1, size(Bs%incoming)
        @assertEqual(0._ireals, Bs%outgoing(iface), rtol, 'outgoing fluxes should be zero because buildings_albedo = 0')

        call ind_1d_to_nd(Bs%da_offsets, Bs%iface(iface), idx)
        associate(d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
          k = local_dims(1)+lk
          i = local_dims(3)+li
          j = local_dims(5)+lj

          print *,myid, 'Building iface', iface, 'idx =>', idx, 'edir', Bs%edir(iface), 'inc', Bs%incoming(iface), 'out', Bs%outgoing(iface)

          if(k.eq.ke .and. i.eq.3 .and. j.eq.3) then ! center pyramid box
            @assertEqual(0._ireals, Bs%edir    (iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
            @assertEqual(0._ireals, Bs%incoming(iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
            @assertEqual(0._ireals, Bs%outgoing(iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
          endif

          if(d.eq.PPRTS_TOP_FACE .and. k.eq.ke-1 .and. i.eq.3 .and. j.eq.3) then ! center pyramid box above
            @assertEqual(gedir(k,i,j), Bs%edir(iface), rtol*gedir(k,i,j), 'edir flux of lifted center building should be same as edir in atmosphere ')
          endif

          if(d.eq.PPRTS_TOP_FACE .and. k.eq.ke .and. i.ne.3 .and. j.ne.3) then ! lower cells without the center one
            @assertEqual(gedir(k,i,j), Bs%edir(iface), rtol*gedir(k,i,j), 'edir flux of surrounding buildings should be same as edir in atmosphere ')
          endif

          if(d.ne.PPRTS_TOP_FACE) then
            @assertEqual(0._ireals, Bs%edir(iface), rtol, 'edir on all but the top face should be zero because sza=0')
          endif
        end associate
      enddo

    end associate
    call destroy_buildings(buildings_solar,  ierr); call CHKERR(ierr)
    call destroy_buildings(buildings_thermal,ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[4,2,1])
  subroutine test_pprts_rrtmg_buildings_example_thermal(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm, myid, ierr

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=10
    real(ireals), parameter :: dx=100, dy=100
    character(len=*), parameter :: atm_filename='afglus_100m.dat'
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: rtol=1e-2_ireals, sb=5.67e-8
    logical :: lsolar, lthermal
    real(ireals) :: buildings_albedo, buildings_temp
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag_solar, Ag_thermal
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso
    type(t_pprts_buildings), allocatable :: buildings_solar, buildings_thermal

    real(ireals) :: trgt
    integer(iintegers) :: local_dims(6)
    integer(iintegers) :: ke, iface, idx(4)
    integer(iintegers) :: k, i, j

    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lsolar   = .False.
    lthermal = .True.
    buildings_albedo = 0
    buildings_temp   = 288._ireals
    phi0   = 0
    theta0 = -1
    Ag_solar   = 0
    Ag_thermal = 0

    do i=1,2
      call ex_pprts_rrtm_buildings(           &
        & comm, lverbose,                     &
        & lthermal, lsolar,                   &
        & Nx, Ny, Nlay,                       &
        & buildings_albedo, buildings_temp,   &
        & dx, dy,                             &
        & atm_filename,                       &
        & phi0, theta0,                       &
        & Ag_solar, Ag_thermal,               &
        & gedir, gedn, geup, gabso,           &
        & buildings_solar, buildings_thermal, &
        & local_dims )
    enddo

    if(myid.eq.0.and.lverbose) then
      do k=ubound(gedn,1)-Nlay, ubound(gedn,1)
        print *,'j(center) k', k, 'edn', gedn(k,:,Ny/2)
      enddo
      do k=ubound(gedn,1)-Nlay, ubound(gedn,1)
        print *,'j(center) k', k, 'eup', geup(k,:,Ny/2)
      enddo
    endif
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    associate( Bs => buildings_solar, Bt => buildings_thermal )
      @mpiassertFalse(allocated(Bt%edir),     'edir in thermal should never be allocated')
      @mpiassertFalse(allocated(Bs%edir),     'buildings_solar should not have results allocated because lsolar=.F.')
      @mpiassertFalse(allocated(Bs%incoming), 'buildings_solar should not have results allocated because lsolar=.F.')
      @mpiassertFalse(allocated(Bs%outgoing), 'buildings_solar should not have results allocated because lsolar=.F.')

      ke = size(gabso,dim=1)

      if(lverbose) then
        do iface = 1, size(Bt%outgoing)
          call ind_1d_to_nd(Bt%da_offsets, Bt%iface(iface), idx)
          associate(d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
            k = local_dims(1)+lk
            i = local_dims(3)+li
            j = local_dims(5)+lj

            print *,myid, 'Building iface', iface, 'idx =>', idx, 'inc', Bt%incoming(iface), 'out', Bt%outgoing(iface)
          end associate
        enddo
      endif

      trgt = sb*buildings_temp**4
      do iface = 1, size(Bt%outgoing)
        @assertEqual(trgt, Bt%outgoing(iface), trgt*rtol, 'outgoing fluxes should be stefan boltzmann emission because buildings_albedo=0:')
      enddo

      do iface = 1, size(Bt%incoming)

        call ind_1d_to_nd(Bt%da_offsets, Bt%iface(iface), idx)
        associate(d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
          k = local_dims(1)+lk
          i = local_dims(3)+li
          j = local_dims(5)+lj

          if(k.eq.ke .and. i.eq.3 .and. j.eq.3) then ! center pyramid box
            @assertEqual(Bt%outgoing(iface), Bt%incoming(iface), rtol, 'all in fluxes of center building should be same as outgoing bc it is encircled by others')
          endif
        end associate
      enddo

    end associate
    call destroy_buildings(buildings_solar,  ierr); call CHKERR(ierr)
    call destroy_buildings(buildings_thermal,ierr); call CHKERR(ierr)
  end subroutine

  !@test(npes =[4,2,1])
  subroutine test_buildings_virtual_faces(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx=6, Ny=6, Nlay=4, icollapse=1
    integer(iintegers), parameter :: glob_box_i=3, glob_box_j=4, glob_box_k=2
    real(ireals), parameter :: dx=100, dy=100, dz=100, S0=1
    logical, parameter :: lverbose=.True.
    real(ireals), parameter :: atol=1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    real(ireals) :: dz1d(Nlay)
    real(ireals) :: sundir(3), trgt
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,plck
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver
    integer(iintegers) :: Nbuildings, Nvirtual
    logical :: lhave_box, lhave_virtual_box

    integer(iintegers) :: k, i
    integer(iintegers) :: box_k, box_i, box_j
    integer(mpiint) :: numnodes, ierr


    comm = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    lthermal = .False.
    lsolar   = .True.
    box_albedo = 0.
    box_planck = 0
    phi0       = 0
    theta0     = 45
    Ag         = 0
    dtau       = 0
    w0         = 0

    dz1d = dz

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(comm, Nlay, Nx, Ny, dx,dy, sundir, solver, dz1d, collapseindex=icollapse)

    associate(Ca => solver%C_one_atm, C1 => solver%C_one)
      allocate(kabs(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(ksca(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(g   (Ca%zm  , Ca%xm, Ca%ym ))

      if(lthermal) then
        allocate(plck(Ca%zm+1, Ca%xm, Ca%ym ))
        plck(:,:,:) = 0
        plck(Ca%zm+1,:,:) = 0
      endif

      kabs = dtau*(one-w0)/dz/real(Nlay, ireals)
      ksca = dtau*w0/dz/real(Nlay, ireals)
      g    = zero

      box_k = glob_box_k - C1%zs
      box_i = glob_box_i - C1%xs
      box_j = glob_box_j - C1%ys

      if(lverbose) print *, myid, 'Have box:', &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze+1), &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe+1), &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye+1)

      if( &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze+1).and. &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe+1).and. &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye+1)      ) then

        lhave_box = .True.
        Nbuildings = 12
      else
        lhave_box = .False.
        Nbuildings = 0
      endif

      if( &
        !& .False. .and. &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze+1).and. &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe+1).and. &
        & is_inrange(glob_box_j-1, C1%ys+1, C1%ye+1)      ) then

        if(lverbose) print *, myid, 'Have virtual box'
        Nvirtual = 3
        lhave_virtual_box = .True.
      else
        Nvirtual = 0
        lhave_virtual_box = .False.
      endif
      Nbuildings = Nbuildings + Nvirtual

      call init_buildings(buildings, &
        & [integer(iintegers) :: 6, C1%zm, C1%xm,  C1%ym], &
        & Nbuildings, &
        & ierr); call CHKERR(ierr)

      if(lthermal) allocate(buildings%planck(Nbuildings))
      if(lhave_box) then
        ! upper box
        buildings%iface( 1) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_TOP_FACE)
        buildings%iface( 2) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_BOT_FACE)
        buildings%iface( 3) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_LEFT_FACE)
        buildings%iface( 4) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_RIGHT_FACE)
        buildings%iface( 5) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_REAR_FACE)
        buildings%iface( 6) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j, PPRTS_FRONT_FACE)

        ! lower box
        buildings%iface( 7) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_TOP_FACE)
        buildings%iface( 8) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_BOT_FACE)
        buildings%iface( 9) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_LEFT_FACE)
        buildings%iface(10) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_RIGHT_FACE)
        buildings%iface(11) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_REAR_FACE)
        buildings%iface(12) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k+1, box_i, box_j, PPRTS_FRONT_FACE)


        buildings%albedo( 1) = box_albedo
        buildings%albedo( 2) = box_albedo
        buildings%albedo( 3) = box_albedo
        buildings%albedo( 4) = box_albedo
        buildings%albedo( 5) = box_albedo
        buildings%albedo( 6) = box_albedo

        ! lower boxbox_albedo
        buildings%albedo( 7) = box_albedo
        buildings%albedo( 8) = box_albedo
        buildings%albedo( 9) = box_albedo
        buildings%albedo(10) = box_albedo
        buildings%albedo(11) = box_albedo
        buildings%albedo(12) = box_albedo


        if(lthermal) then
          buildings%planck(1:12) = box_planck
        endif
      endif
      if(lhave_virtual_box) then
        ! virtual face to get the side flux but not change anything in the RT
        k = size(buildings%iface) - Nvirtual+1
        buildings%iface(k+0) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j-1, PPRTS_TOP_FACE)
        buildings%iface(k+1) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j-1, PPRTS_BOT_FACE)
        buildings%iface(k+2) = faceidx_by_cell_plus_offset(buildings%da_offsets, box_k, box_i, box_j-1, PPRTS_FRONT_FACE)

        buildings%albedo(k:k+Nvirtual-1) = -1
        if(lthermal) buildings%planck(k:k+Nvirtual-1) = -1
      endif

      call check_buildings_consistency(buildings, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)

      if(lthermal) then
        call set_optical_properties(solver, Ag, kabs, ksca, g, plck)
      else
        call set_optical_properties(solver, Ag, kabs, ksca, g)
      endif
      call set_angles(solver, sundir)

      call solve_pprts(solver, &
        & lthermal=lthermal, &
        & lsolar=lsolar, &
        & edirTOA=S0, &
        & opt_buildings=buildings)

      call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_buildings=buildings)

      if(lsolar) then
        call gather_all_toZero(solver%C_one_atm1, fdir, gedir)
      endif
      call gather_all_toZero(solver%C_one_atm1, fdn, gedn)
      call gather_all_toZero(solver%C_one_atm1, fup, geup)
      call gather_all_toZero(solver%C_one_atm, fdiv, gabso)

      if(lverbose .and. myid.eq.0) then
        print *,'y-slice: i='//toStr(box_i)
        i = box_i
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, i, :) ), 'red')
          enddo
        endif ! lsolar
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, i, :) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, i, :) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, i, :) ), 'purple')
        enddo

        print *,''
        print *,'x-slice: j='//toStr(box_j)
        i = box_j
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, :, i) ), 'red')
          enddo
        endif
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, :, i) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, :, i) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, :, i) ), 'purple')
        enddo

        print *,''
        if(lsolar) then
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, 'mean ', &
              & 'edir', meanval(gedir(k,:,:)), &
              & 'edn' , meanval(gedn(k,:,:)), &
              & 'eup' , meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edir', meanval(gedir(k,:,:)), &
            & 'edn' , meanval(gedn(k,:,:)), &
            & 'eup' , meanval(geup(k,:,:))
        else
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, &
              & 'mean ', &
              & 'edn', meanval(gedn(k,:,:)), &
              & 'eup', meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edn', meanval(gedn(k,:,:)), &
            & 'eup', meanval(geup(k,:,:))
        endif
      endif
    end associate
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    if(lverbose) then
      do k = 0, numnodes-1
        if(k.eq.myid) then
          if(allocated(buildings%edir)) then
            print *,cstr(' * Buildings on rank '//toStr(k)//':', 'red')
            do i=1, size(buildings%iface)
              print *, 'building_face', i, 'edir', buildings%edir(i), &
                & 'in/out', buildings%incoming(i), buildings%outgoing(i)
            enddo
          else
            do i=1, size(buildings%iface)
              print *, 'building_face', i, &
                & 'in/out', buildings%incoming(i), buildings%outgoing(i)
            enddo
          endif
        endif
        call mpi_barrier(comm, ierr); call CHKERR(ierr)
      enddo
    endif

    call destroy_pprts(solver, .False.)

    if(myid.eq.0) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay+1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j  ), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k+1, glob_box_i, glob_box_j-1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j-1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j-2), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+3, glob_box_i, glob_box_j-2), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+3, glob_box_i, glob_box_j-3), atol, 'edir building shadow should go south and should be zero')

      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j  ), atol, 'edir at bot of lower building should be zero')
      @assertEqual(0, gedir(glob_box_k+2, glob_box_i, glob_box_j-1), atol, 'lower building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+3, glob_box_i, glob_box_j-1), atol, 'lower building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k+3, glob_box_i, glob_box_j-2), atol, 'lower building shadow should go south and should be zero')
    endif

    if(lhave_box) then
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0*cos(deg2rad(theta0)) ! front face bc phi==0 is north)
      @assertEqual(trgt, buildings%edir(6), atol)

      @assertEqual(0._ireals, buildings%edir(2:5), atol)

      ! second building block below, only have on sunlit face on the side
      @assertEqual(0._ireals, buildings%edir(7:11), atol)
      @assertEqual(trgt, buildings%edir(12), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    endif

    if(lhave_virtual_box) then
      k = size(buildings%edir) - Nvirtual+1
      trgt = S0*cos(deg2rad(theta0))
      @assertEqual(trgt     , buildings%edir(k+0), atol)
      @assertEqual(0._ireals, buildings%edir(k+1), atol)
      @assertEqual(0._ireals, buildings%edir(k+2), atol)
    endif
  end subroutine
end module
