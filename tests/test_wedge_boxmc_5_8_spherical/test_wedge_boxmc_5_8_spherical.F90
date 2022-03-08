module test_wedge_boxmc_5_8_spherical
  use m_boxmc, only: t_boxmc, t_boxmc_wedge_5_8
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_helper_functions, only: itoa, triangle_area_by_vertices
  use m_boxmc_geometry, only: setup_default_unit_wedge_geometry, setup_default_wedge_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(8), T(5), S_target(8), T_target(5)
  real(ireals) :: S_tol(8), T_tol(5)

  type(t_boxmc_wedge_5_8) :: bmc_wedge_5_8

  integer(mpiint) :: myid, mpierr, numnodes, comm

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
  !real(ireal_dp),parameter :: atol=1e-4, rtol=1e-2
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_wedge_5_8%init(comm)

    if (myid .eq. 0) print *, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 45

    dx = 100
    dy = dx
    dz = 50

    S_target = zero
    T_target = zero

  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    myid = this%getProcessRank()
    call PetscFinalize(ierr)
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_direct_src1(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 1
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp) :: Atop, Abot
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], [dx, 0._ireal_dp], [dx / 2, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    Abot = triangle_area_by_vertices(sverts(1:3), sverts(4:6), sverts(7:9))
    Atop = triangle_area_by_vertices(sverts(10:12), sverts(13:15), sverts(16:18))

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    phi = 0; theta = 0
    T_target = zero
    T_target(5) = real(exp(-(bg(1) + bg(2)) * dz) * Abot / Atop, ireals)
    T_target(2:4) = 0.06665
    S_target = [0.0026, 0.0055, 0.0016, 0.0056, 0.0016, 0.0056, 0.0016, 0.0175]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src1_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.009, 0.0046, 0.0040, 0.0046, 0.0040, 0.0046, 0.0040, 0.0071]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src1_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_direct_src5(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 5
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp) :: Atop, Abot
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], [dx, 0._ireal_dp], [dx / 2, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    Abot = triangle_area_by_vertices(sverts(1:3), sverts(4:6), sverts(7:9))
    Atop = triangle_area_by_vertices(sverts(10:12), sverts(13:15), sverts(16:18))

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    phi = 0; theta = 180
    T_target = zero
    T_target(1) = real(exp(-(bg(1) + bg(2)) * dz) * Abot / Atop, ireals)
    T_target(2:4) = 0.06665
    S_target = [0.0175, 0.0016, 0.0055, 0.0016, 0.0055, 0.0016, 0.0055, 0.0026]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src1_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.0071, 0.0040, 0.0046, 0.0040, 0.0046, 0.0040, 0.0046, 0.0090]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src1_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_direct_src2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 2
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374._ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    ! straight outwards from face 2
    phi = 0; theta = 90

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    T_target = zero
    T_target([3, 4]) = 0.4548_ireals

    S_target = [0.0131, 0.0012, 0.0018, 0.0024, 0.0108, 0.0024, 0.0108, 0.0012]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src2_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.015, 0.0029, 0.0065, 0.0063, 0.0027, 0.0063, 0.0027, 0.0015]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src2_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_direct_src3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 3
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374._ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    ! straight outwards from face 2
    phi = 120; theta = 90

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    T_target = zero
    T_target([2, 4]) = 0.4548_ireals

    S_target = [0.0131, 0.0024, 0.0108, 0.0012, 0.0018, 0.0024, 0.0108, 0.0012]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src3_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.015, 0.0063, 0.0027, 0.0029, 0.0065, 0.0063, 0.0027, 0.0015]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src3_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_direct_src4(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 4
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374._ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    ! straight outwards from face 2
    phi = 240; theta = 90

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    T_target = zero
    T_target([2, 3]) = 0.4548_ireals

    S_target = [0.0131, 0.0024, 0.0108, 0.0024, 0.0108, 0.0012, 0.0018, 0.0012]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src4_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.015, 0.0063, 0.0027, 0.0063, 0.0027, 0.0029, 0.0065, 0.0015]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .true., phi, theta, sverts, &
                                 S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_direct_src4_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src1(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 1
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], [dx, 0._ireal_dp], [dx / 2, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    S_target = [0.0000, 0.2455, 0.0015, 0.2455, 0.0015, 0.2455, 0.0015, 0.2220]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src1_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src1_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    S_target = [0.0041, 0.2437, 0.0033, 0.2437, 0.0033, 0.2437, 0.0033, 0.2179]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src1_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 2
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.0019, 0.0000, 0.0000, 0.2362, 0.0712, 0.2362, 0.0712, 0.3474]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src2_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src2_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.0075, 0.0043, 0.0052, 0.2301, 0.0716, 0.2302, 0.0716, 0.3435]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src2_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 3
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.4775, 0.0000, 0.0000, 0.0000, 0.2450, 0.0000, 0.2450, 0.0000]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src3_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src3_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.4702, 0.0047, 0.0043, 0.0029, 0.2388, 0.0029, 0.2388, 0.0042]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src3_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src4(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 4
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.0019, 0.2362, 0.0712, 0.0000, 0.0000, 0.2362, 0.0712, 0.3474]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src4_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src4_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.0075, 0.2301, 0.0716, 0.0043, 0.0052, 0.2302, 0.0716, 0.3435]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src4_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src5(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 5
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.4775, 0.0000, 0.2450, 0.0000, 0.0000, 0.0000, 0.2450, 0.0000]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src5_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src5_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.4702, 0.0029, 0.2388, 0.0047, 0.0043, 0.0029, 0.2388, 0.0042]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src5_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src6(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 6
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.0019, 0.2362, 0.0712, 0.2362, 0.0712, 0.0000, 0.0000, 0.3474]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src6_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src6_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.0075, 0.2301, 0.0716, 0.2302, 0.0716, 0.0043, 0.0052, 0.3435]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src6_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src7(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 7
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([-dx / 2, 0._ireal_dp], [dx / 2, 0._ireal_dp], [0._ireal_dp, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    S_target = [0.4775, 0.0000, 0.2450, 0.0000, 0.2450, 0.0000, 0.0000, 0.0000]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src7_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src7_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    S_target = [0.4702, 0.0029, 0.2388, 0.0029, 0.2388, 0.0047, 0.0043, 0.0042]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src7_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_spherical_diffuse_src8(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 8
    real(ireal_dp), allocatable :: sverts(:)
    real(ireal_dp), parameter :: R = 6374e1_ireal_dp

    call setup_default_wedge_geometry([0._ireal_dp, 0._ireal_dp], [dx, 0._ireal_dp], [dx / 2, sqrt(dy**2 - (dx / 2)**2)], &
                                      dz, sverts, sphere_radius=R / dx)

    T_target = zero

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    S_target = [0.2812, 0.0000, 0.268, 0.0000, 0.268, 0.0000, 0.268, 0.0000]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src8_1')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src8_2')

    ! with strict forward scattering the same as without ksca
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    S_target = [0.2759, 0.0014, 0.2259, 0.0014, 0.2259, 0.0014, 0.2259, 0.0038]

    call bmc_wedge_5_8%get_coeff(comm, bg, src, .false., phi, theta, sverts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_spherical_diffuse_src8_3')
  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

    real(ireal_dp), parameter :: sigma = 3 ! normal test range for coefficients

    character(len=*), optional :: msg
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals) :: test_atol

    if (myid .eq. 0) then
      print *, ''

      if (present(msg)) then
        write (local_msgS, *) trim(msg), ':: Diffuse boxmc coefficient not as '
        write (local_msgT, *) trim(msg), ':: Direct  boxmc coefficient not as '
        print *, msg
      else
        write (local_msgS, *) 'Diffuse boxmc coefficient not as '
        write (local_msgT, *) 'Direct  boxmc coefficient not as '
      end if

      print *, '---------------------'
      write (*, FMT='( " diffuse ::  :: ",10(es12.5) )') S
      write (*, FMT='( " target  ::  :: ",10(es12.5) )') S_target
      write (*, FMT='( " diff    ::  :: ",10(es12.5) )') S_target - S
      print *, ''
      write (*, FMT='( " direct  ::  :: ", 8(es12.5) )') T
      write (*, FMT='( " target  ::  :: ", 8(es12.5) )') T_target
      write (*, FMT='( " diff    ::  :: ", 8(es12.5) )') T_target - T
      print *, '---------------------'
      print *, ''

      test_atol = real(atol, ireals) * real(sigma, ireals)

      @assertEqual(S_target, S, test_atol, local_msgS)
      @assertLessThanOrEqual(zero, S)
      @assertGreaterThanOrEqual(one, S)

      @assertEqual(T_target, T, test_atol, local_msgT)
      @assertLessThanOrEqual(zero, T)
      @assertGreaterThanOrEqual(one, T)
    end if
  end subroutine

end module
