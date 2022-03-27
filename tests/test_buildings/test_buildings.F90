module test_buildings
#include "petsc/finclude/petsc.h"
  use petsc
  use iso_fortran_env, only: real32, real64
  use m_data_parameters, only: &
    init_mpi_data_parameters, &
    finalize_mpi, &
    iintegers, ireals, mpiint, &
    zero, one, pi, default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & deg2rad, &
    & ind_1d_to_nd

  ! main entry point for solver, and desctructor
!  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

!  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

!  use m_pprts_base, only : t_solver_3_10
  use m_buildings, only: &
    & t_pprts_buildings, &
    & init_buildings,    &
    & clone_buildings,   &
    & destroy_buildings, &
    & PPRTS_TOP_FACE

  use m_examples_pprts_buildings, only: ex_pprts_buildings
  use m_examples_pprts_rrtm_buildings, only: ex_pprts_rrtm_buildings

  use pfunit_mod

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
    ! Tidy up
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes=[1])
  subroutine test_init_buildings(this)
    class(MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    integer(iintegers), parameter :: Nfaces = 6
    integer(iintegers), parameter :: da_sizes(4) = [integer(iintegers) :: 6, 1, 1, 1]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings(B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated(B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    end do
    call check_iface_data_consistency(B1)

    call destroy_buildings(B1, ierr); call CHKERR(ierr)
    @assertFalse(allocated(B1))

  contains
    subroutine check_iface_data_consistency(B)
      type(t_pprts_buildings), intent(in) :: B
      integer(iintegers) :: i
      do i = 1, Nfaces
        @assertEqual(i, B%iface(i))
        @assertEqual(real(i, ireals), B%albedo(i))
      end do
    end subroutine
  end subroutine

  @test(npes=[1])
  subroutine test_clone_buildings(this)
    class(MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    type(t_pprts_buildings), allocatable :: B2
    type(t_pprts_buildings), allocatable :: B3
    integer(iintegers), parameter :: Nfaces = 6
    integer(iintegers), parameter :: da_sizes(4) = [integer(iintegers) :: 6, 1, 1, 1]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings(B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated(B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    end do
    call check_iface_data_consistency(B1)

    call clone_buildings(B1, B2, .true., ierr); call CHKERR(ierr)
    @assertEqual(2_iintegers, B1%iface_data%ref_count, 'origin ref_count did not increase')

    @assertTrue(allocated(B2%albedo))
    @assertTrue(associated(B2%iface))
    @assertTrue(associated(B2%iface_data%data))
    @assertEqual(Nfaces, size(B2%iface_data%data, kind=iintegers))
    @assertEqual(2_iintegers, B2%iface_data%ref_count)
    @assertEqual(B1%iface, B2%iface)
    @assertTrue(associated(B1%iface, B2%iface), 'cloned iface pointer dont point to same target')
    call check_iface_data_consistency(B2)

    call clone_buildings(B1, B3, .false., ierr); call CHKERR(ierr)
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
    @assertFalse(allocated(B1))

    call check_iface_data_consistency(B2)
    call check_iface_data_consistency(B3)

    call destroy_buildings(B2, ierr); call CHKERR(ierr)
    @assertFalse(allocated(B2))
    call destroy_buildings(B3, ierr); call CHKERR(ierr)
    @assertFalse(allocated(B3))

  contains
    subroutine check_iface_data_consistency(B)
      type(t_pprts_buildings), intent(in) :: B
      integer(iintegers) :: i
      do i = 1, Nfaces
        @assertEqual(i, B%iface(i), "iface does not match")
        if (allocated(B%albedo)) then
          @assertEqual(real(i, ireals), B%albedo(i), "albedo does not match")
        end if
      end do
    end subroutine
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_overhead_sun(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0
    box_planck = 0
    phi0 = 0
    theta0 = 0
    Ag = 0
    dtau = 0
    w0 = 0

    call ex_pprts_buildings(comm, lverbose,         &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,  &
      & glob_box_i, glob_box_j, glob_box_k,         &
      & box_Ni, box_nj, box_Nk,                     &
      & box_albedo, box_planck,                     &
      & dx, dy, dz,                                 &
      & S0, phi0, theta0,                           &
      & Ag, dtau, w0,                               &
      & gedir, gedn, geup, gabso,                   &
      & buildings)

    if (myid .eq. 0) then
      @assertEqual(0, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(S0 / dz, gabso(glob_box_k, glob_box_i, glob_box_j), atol, 'box should absorb all incoming solar radiation')
    end if
    if (size(buildings%edir) .gt. 0) then
      @assertEqual(S0, buildings%edir(1), atol)
      @assertEqual(0._ireals, buildings%edir(2:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_albedo(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0.5_ireals
    box_planck = 0
    phi0 = 0
    theta0 = 0
    Ag = 1.0_ireals
    dtau = 0
    w0 = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      @assertEqual(0, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(S0, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be S0')

      @assertEqual(S0 * box_albedo, geup(glob_box_k, glob_box_i, glob_box_j), atol, 'eup at top of building should be S0*B_Ag')

    end if
    if (size(buildings%edir) .gt. 0) then
      @assertEqual(S0, buildings%edir(1), atol, 'edir at top of building should be S0')
      @assertEqual(0._ireals, buildings%edir(2:6), atol, 'if zenith angle is 0, all sides of the building except the top should be 0 edir')

      @assertEqual((buildings%edir+buildings%incoming)*box_albedo, buildings%outgoing, atol, 'total incoming times building albedo should give outgoing')
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_emission(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .true.
    lsolar = .false.
    box_albedo = 0._ireals
    box_planck = 1
    phi0 = 0
    theta0 = 180
    Ag = 0
    dtau = 0
    w0 = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      @assertEqual(box_planck * pi, geup(glob_box_k, glob_box_i, glob_box_j), atol, 'eup at top of building should be emission')
      @assertEqual(box_planck * pi, gedn(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edn at top of building should be emission')
    end if
    @assertEqual(box_planck * pi, buildings%outgoing, atol, 'emission on buildings should be planck')
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_fortyfive_sun_azi0(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0
    box_planck = 0
    phi0 = 0
    theta0 = 45
    Ag = 0
    dtau = 0
    w0 = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j - 1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i, glob_box_j - 1), atol, 'edir building shadow should go south and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i, glob_box_j - 2), atol, 'edir building shadow should go south and should be zero')
    end if
    if (size(buildings%edir) .gt. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0 * cos(deg2rad(theta0)) ! front face bc phi==0 is north)
      @assertEqual(trgt, buildings%edir(6), atol)

      @assertEqual(0._ireals, buildings%edir(2:5), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_fortyfive_sun_azi180(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0
    box_planck = 0
    phi0 = 180
    theta0 = 45
    Ag = 0
    dtau = 0
    w0 = 0

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j + 1), atol, 'edir building shadow should go north and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i, glob_box_j + 1), atol, 'edir building shadow should go north and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i, glob_box_j + 2), atol, 'edir building shadow should go north and should be zero')
    end if
    if (size(buildings%edir) .gt. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0 * cos(deg2rad(theta0)) ! back face bc phi==180 is south)
      @assertEqual(trgt, buildings%edir(5), atol)

      @assertEqual(0._ireals, buildings%edir(2:4), atol)
      @assertEqual(0._ireals, buildings%edir(6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_fortyfive_sun_azi90(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0
    box_planck = 0
    theta0 = 45
    Ag = 0
    dtau = 0
    w0 = 0
    phi0 = 90

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i - 1, glob_box_j), atol, 'edir building shadow should go west and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i - 1, glob_box_j), atol, 'edir building shadow should go west and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i - 2, glob_box_j), atol, 'edir building shadow should go west and should be zero')
    end if
    if (size(buildings%edir) .gt. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0 * cos(deg2rad(theta0)) ! right face bc phi==90 is east sun)
      @assertEqual(trgt, buildings%edir(4), atol)

      @assertEqual(0._ireals, buildings%edir(2:3), atol)
      @assertEqual(0._ireals, buildings%edir(5:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_buildings_example_fortyfive_sun_azi270(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 3, icollapse = 1
    integer(iintegers), parameter :: glob_box_i = 3, glob_box_j = 3, glob_box_k = 2
    integer(iintegers), parameter :: box_Ni = 1, box_nj = 1, box_Nk = 1
    real(ireals), parameter :: dx = 100, dy = 100, dz = 100, S0 = 1
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: atol = 1e-5_ireals

    logical :: lsolar, lthermal
    real(ireals) :: box_albedo, box_planck
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag, dtau, w0
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals) :: trgt

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lthermal = .false.
    lsolar = .true.
    box_albedo = 0
    box_planck = 0
    theta0 = 45
    Ag = 0
    dtau = 0
    w0 = 0
    phi0 = 270

    call ex_pprts_buildings(comm, lverbose,       &
      & lthermal, lsolar, Nx, Ny, Nlay, icollapse,&
      & glob_box_i, glob_box_j, glob_box_k,       &
      & box_Ni, box_nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings)

    if (myid .eq. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, gedir(Nlay + 1, glob_box_i, glob_box_j), atol, 'edir beneath building should be clear sky value')
      @assertEqual(trgt, gedir(glob_box_k, glob_box_i, glob_box_j), atol, 'edir at top of building should be the one from clear sky')

      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i, glob_box_j), atol, 'edir at bot of building should be zero')
      @assertEqual(0, gedir(glob_box_k + 1, glob_box_i + 1, glob_box_j), atol, 'edir building shadow should go east and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i + 1, glob_box_j), atol, 'edir building shadow should go east and should be zero')
      @assertEqual(0, gedir(glob_box_k + 2, glob_box_i + 2, glob_box_j), atol, 'edir building shadow should go east and should be zero')
    end if
    if (size(buildings%edir) .gt. 0) then
      trgt = S0 * cos(deg2rad(theta0))
      @assertEqual(trgt, buildings%edir(1), atol)
      trgt = S0 * cos(deg2rad(theta0)) ! leftt face bc phi==270 is west sun)
      @assertEqual(trgt, buildings%edir(3), atol)

      @assertEqual(0._ireals, buildings%edir(2), atol)
      @assertEqual(0._ireals, buildings%edir(4:6), atol)

      @assertEqual(0._ireals, buildings%incoming, atol)
      @assertEqual(0._ireals, buildings%outgoing, atol)
    end if
  end subroutine

  @test(npes=[4, 1])
  subroutine test_pprts_rrtmg_buildings_example_solar(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm, myid, ierr

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 10
    real(ireals), parameter :: dx = 100, dy = 100
    character(len=*), parameter :: atm_filename = 'afglus_100m.dat'
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: rtol = 1e-5_ireals
    logical :: lsolar, lthermal
    real(ireals) :: buildings_albedo, buildings_temp
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag_solar, Ag_thermal
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso
    type(t_pprts_buildings), allocatable :: buildings_solar, buildings_thermal

    integer(iintegers) :: local_dims(6)
    integer(iintegers) :: ke, iface, idx(4)
    integer(iintegers) :: k, i, j

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lsolar = .true.
    lthermal = .false.
    buildings_albedo = 0
    buildings_temp = 0
    phi0 = 0
    theta0 = 0
    Ag_solar = 0
    Ag_thermal = 0

    do i = 1, 2
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
        & local_dims)
    end do

    print *, myid, 'shape edir', shape(gedir)

    associate (Bs => buildings_solar, Bt => buildings_thermal)
      @assertFalse(allocated(Bt%edir), 'edir in thermal should never be allocated')
      @assertFalse(allocated(Bt%incoming), 'buildings_thermal should not have results allocated because lthermal=.F. ')
      @assertFalse(allocated(Bt%outgoing), 'buildings_thermal should not have results allocated because lthermal=.F. ')

      ke = size(gabso, dim=1)

      do iface = 1, size(Bs%incoming)
        @assertEqual(0._ireals, Bs%outgoing(iface), rtol, 'outgoing fluxes should be zero because buildings_albedo = 0')

        call ind_1d_to_nd(Bs%da_offsets, Bs%iface(iface), idx)
        associate (d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
          k = local_dims(1) + lk
          i = local_dims(3) + li
          j = local_dims(5) + lj

          print *, myid, 'Building iface', iface, 'idx =>', idx, &
            & 'edir', Bs%edir(iface), 'inc', Bs%incoming(iface), 'out', Bs%outgoing(iface)

          if (k .eq. ke .and. i .eq. 3 .and. j .eq. 3) then ! center pyramid box
            @assertEqual(0._ireals, Bs%edir(iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
            @assertEqual(0._ireals, Bs%incoming(iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
            @assertEqual(0._ireals, Bs%outgoing(iface), rtol, 'all fluxes of center building should be zero bc it is encircled by others')
          end if

          if (d .eq. PPRTS_TOP_FACE .and. k .eq. ke - 1 .and. i .eq. 3 .and. j .eq. 3) then ! center pyramid box above
            @assertEqual(gedir(k,i,j), Bs%edir(iface), rtol*gedir(k,i,j), 'edir flux of lifted center building should be same as edir in atmosphere ')
          end if

          if (d .eq. PPRTS_TOP_FACE .and. k .eq. ke .and. i .ne. 3 .and. j .ne. 3) then ! lower cells without the center one
            @assertEqual(gedir(k,i,j), Bs%edir(iface), rtol*gedir(k,i,j), 'edir flux of surrounding buildings should be same as edir in atmosphere ')
          end if

          if (d .ne. PPRTS_TOP_FACE) then
            @assertEqual(0._ireals, Bs%edir(iface), rtol, 'edir on all but the top face should be zero because sza=0')
          end if
        end associate
      end do

    end associate
    call destroy_buildings(buildings_solar, ierr); call CHKERR(ierr)
    call destroy_buildings(buildings_thermal, ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[4, 1])
  subroutine test_pprts_rrtmg_buildings_example_thermal(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm, myid, ierr

    integer(iintegers), parameter :: Nx = 6, Ny = 6, Nlay = 10
    real(ireals), parameter :: dx = 100, dy = 100
    character(len=*), parameter :: atm_filename = 'afglus_100m.dat'
    logical, parameter :: lverbose = .true.
    real(ireals), parameter :: rtol = 1e-2_ireals, sb = 5.67e-8
    logical :: lsolar, lthermal
    real(ireals) :: buildings_albedo, buildings_temp
    real(ireals) :: phi0, theta0
    real(ireals) :: Ag_solar, Ag_thermal
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso
    type(t_pprts_buildings), allocatable :: buildings_solar, buildings_thermal

    real(ireals) :: trgt
    integer(iintegers) :: local_dims(6)
    integer(iintegers) :: ke, iface, idx(4)
    integer(iintegers) :: k, i, j

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    lsolar = .false.
    lthermal = .true.
    buildings_albedo = 0
    buildings_temp = 288._ireals
    phi0 = 0
    theta0 = -1
    Ag_solar = 0
    Ag_thermal = 0

    do i = 1, 2
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
        & local_dims)
    end do

    if (myid .eq. 0 .and. lverbose) then
      do k = ubound(gedn, 1) - Nlay, ubound(gedn, 1)
        print *, 'j(center) k', k, 'edn', gedn(k, :, Ny / 2)
      end do
      do k = ubound(gedn, 1) - Nlay, ubound(gedn, 1)
        print *, 'j(center) k', k, 'eup', geup(k, :, Ny / 2)
      end do
    end if
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    associate (Bs => buildings_solar, Bt => buildings_thermal)
      @mpiassertFalse(allocated(Bt%edir), 'edir in thermal should never be allocated')
      @mpiassertFalse(allocated(Bs%edir), 'buildings_solar should not have results allocated because lsolar=.F.')
      @mpiassertFalse(allocated(Bs%incoming), 'buildings_solar should not have results allocated because lsolar=.F.')
      @mpiassertFalse(allocated(Bs%outgoing), 'buildings_solar should not have results allocated because lsolar=.F.')

      ke = size(gabso, dim=1)

      if (lverbose) then
        do iface = 1, size(Bt%outgoing)
          call ind_1d_to_nd(Bt%da_offsets, Bt%iface(iface), idx)
          associate (d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
            k = local_dims(1) + lk
            i = local_dims(3) + li
            j = local_dims(5) + lj

            print *, myid, 'Building iface', iface, 'idx =>', idx, 'inc', Bt%incoming(iface), 'out', Bt%outgoing(iface)
          end associate
        end do
      end if

      trgt = sb * buildings_temp**4
      do iface = 1, size(Bt%outgoing)
@assertEqual(trgt, Bt%outgoing(iface), trgt*rtol, 'outgoing fluxes should be stefan boltzmann emission because buildings_albedo=0:')
      end do

      do iface = 1, size(Bt%incoming)

        call ind_1d_to_nd(Bt%da_offsets, Bt%iface(iface), idx)
        associate (d => idx(1), lk => idx(2), li => idx(3), lj => idx(4))
          k = local_dims(1) + lk
          i = local_dims(3) + li
          j = local_dims(5) + lj

          if (k .eq. ke .and. i .eq. 3 .and. j .eq. 3) then ! center pyramid box
            @assertEqual(Bt%outgoing(iface), Bt%incoming(iface), rtol, 'all in fluxes of center building should be same as outgoing bc it is encircled by others')
          end if
        end associate
      end do

    end associate
    call destroy_buildings(buildings_solar, ierr); call CHKERR(ierr)
    call destroy_buildings(buildings_thermal, ierr); call CHKERR(ierr)
  end subroutine
end module
