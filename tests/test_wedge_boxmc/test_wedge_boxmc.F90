module test_wedge_boxmc
  use m_boxmc, only: t_boxmc, t_boxmc_wedge_5_5
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_wedge_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(5), T(5), S_target(5), T_target(5)
  real(ireals) :: S_tol(5), T_tol(5)
  real(ireal_dp) :: vertices(18)

  type(t_boxmc_wedge_5_5) :: bmc_wedge_5_5

  integer(mpiint) :: myid, mpierr, numnodes, comm

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
  !real(ireal_dp),parameter :: atol=1e-4, rtol=1e-3
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_wedge_5_5%init(comm)

    if (myid .eq. 0) print *, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 45

    dx = 100
    dy = dx
    dz = 50

    call setup_default_unit_wedge_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero

    ! computed target with stddev_atol=5e-6, stddev_rtol=1e-4 in optprop_parameters
    ! inp_atol=1e-6_ireal_dp, inp_rtol=1e-4_ireal_dp) !
    !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,vertices,S,T,S_tol,T_tol,inp_atol=1e-6_ireal_dp, inp_rtol=1e-4_ireal_dp) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    logical :: lpetsc_is_initialized
    integer(mpiint) :: ierr
    call PetscInitialized(lpetsc_is_initialized, ierr)
    if (lpetsc_is_initialized) call PetscFinalize(ierr)
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1, 2])
  subroutine test_wedgemc_direct_negative_azimuth_src2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 2
    real(ireal_dp) :: tau

    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]
    S_target = zero

    tau = (bg(1) + bg(2)) * sqrt(dy**2 - (dx / 2)**2)

    phi = 60; theta = 90
    T_target = zero
    T_target([4]) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_negative_azimuth_src2_60')

    phi = -60; theta = 90
    T_target = zero
    T_target([3]) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_negative_azimuth_src2_-60')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_direct_src4(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 4
    real(ireal_dp) :: tau

    ! direct to diffuse tests

    ! down along face 4
    phi = 240; theta = 0

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    tau = (bg(1) + bg(2)) * dz
    T_target = zero
    T_target(5) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    S_target = [0.00017051, 0.00039874, 0.00066926, 0.0208261, 0.00203671]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src4_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.00060366, 0.00064998, 0.00100554, 0.0208439, 0.00097272]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src4_1')

    ! outwards from face 4
    phi = 240; theta = 90

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    tau = (bg(1) + bg(2)) * sqrt(dy**2 - (dx / 2)**2)
    T_target = zero
    T_target([2, 3]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)

    S_target = [0.00624004, 0.0124557, 0.0124458, 0.00244814, 0.00624669]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src4_4')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.0076486, 0.00809211, 0.00811609, 0.00836656, 0.00764847]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src4_3')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_direct_src3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 3
    real(ireal_dp) :: tau

    ! down along face 3
    phi = 120; theta = 0

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    tau = (bg(1) + bg(2)) * dz
    T_target = zero
    T_target(5) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    S_target = [0.00017051, 0.00039874, 0.0208261, 0.00066926, 0.00203671]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src3_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.00060366, 0.00064998, 0.0208439, 0.00100554, 0.00097272]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src3_1')

    ! outwards from face 3
    phi = 120; theta = 90

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]
    tau = (bg(1) + bg(2)) * sqrt(dy**2 - (dx / 2)**2)
    T_target = zero
    T_target([2, 4]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)

    S_target = [0.00624004, 0.0124557, 0.00244814, 0.0124458, 0.00624669]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src3_4')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.0076486, 0.00809211, 0.00836656, 0.00811609, 0.00764847]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src3_3')

  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_direct_src2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 2
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    ! outwards from face 2
    phi = 0; theta = 90

    tau = (bg(1) + bg(2)) * sqrt(dy**2 - (dx / 2)**2)
    T_target = zero
    T_target([3, 4]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)

    S_target = [0.00623344, 0.00244899, 0.01244585, 0.0124478, 0.00623946]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src2_4')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.00764423, 0.00835319, 0.00810088, 0.00809931, 0.00764853]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src2_3')

    ! down along face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz
    T_target = zero
    T_target(5) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    S_target = zero; S_target(2) = 0.0241891

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src2_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = zero; S_target(2) = 0.0241891

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src2_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_direct_src1(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 1

    ! direct to diffuse tests, straight down
    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp / 2]

    phi = 0; theta = 0
    T_target = zero; T_target(5) = real(exp(-(bg(1) + bg(2)) * dz), ireals)
    S_target = [0.00262582, 0.0074815, 0.0074815, 0.0074815, 0.0213178]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src1_2')

    bg = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]
    S_target = [0.0090376, 0.009520906, 0.009520906, 0.009520906, 0.00875678]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_src1_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_direct_lambert_beer(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct tests
    bg = [1e-3, 0., 0.]
    phi = 0; theta = 0
    S_target = zero

    !Should be rotationally symmetric for sza=0
    do iphi = 0, 360, 10
      phi = real(iphi, ireal_dp)
      T_target = zero

      print *, 'downward'
      ! Integral from top face, towards the bottom face
      do src = 1, 1
        T_target(5) = real(exp(-(bg(1) + bg(2)) * dz), ireals)
        call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer')
      end do

      ! Integral along each of the faces, towards the bottom face
      do src = 2, 4
        print *, 'downward along sides', src
        tau = bg(1) * dz
        T_target(5) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)
        call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer')
      end do
    end do

    print *, 'upward'
    ! and the same for upward propagation
    theta = 180
    do iphi = 0, 360, 10
      phi = real(iphi, ireal_dp)
      T_target = zero

      do src = 5, 5
        T_target(1) = real(exp(-(bg(1) + bg(2)) * dz), ireals)
        call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_upward')
      end do

      print *, 'upward Integral along each of the faces, towards the bottom face'
      ! Integral along each of the faces, towards the bottom face
      do src = 2, 4
        tau = bg(1) * dz
        T_target(1) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)
        call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_upward')
      end do
    end do

    ! One check that if we start from a side face with 90 degree zenith, we should have equally much on the two opposite faces
    T_target = zero
    phi = 0; theta = 90
    src = 2

    tau = bg(1) * sqrt(dy**2 - (dx / 2)**2)
    t_target([3, 4]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_sidewards')

    phi = 120
    src = 3
    T_target = zero
    T_target([2, 4]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_sidewards')

    phi = 240
    src = 4
    T_target = zero
    T_target([2, 3]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_sidewards')

    ! Or start the photons at the top and they should still go to the side faces
    T_target = zero
    T_target([3, 4]) = real(4.85805e-01 + 4.85883e-01, ireals) / 2
    phi = 0; theta = 90
    src = 1
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 101')
    @assertEqual(T(3), T(4), 3 * real(atol, ireals), 'stream should be same 101')

    T_target = zero
    T_target([2, 4]) = real(4.85805e-01 + 4.85883e-01, ireals) / 2
    phi = 120; theta = 90
    src = 1
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 102')
    @assertEqual(T(2), T(4), 3 * real(atol, ireals), 'stream should be same 102')

    T_target = zero
    T_target([2, 3]) = real(4.85805e-01 + 4.85883e-01, ireals) / 2
    phi = 240; theta = 90
    src = 1
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 103')
    @assertEqual(T(2), T(3), 3 * real(atol, ireals), 'stream should be same 103')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_diffuse_src1(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 1

    ! ----------------------------------
    bg = [1e-3, 1e-2, .0]
    S_target = [0.090919, 0.254331, 0.254331, 0.254331, 0.111285]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src1_2')

    ! ----------------------------------
    bg = [1e-3, 0., .0]
    S_target = [0., 0.27299, 0.27299, 0.27299, 0.146173]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src1_1')

  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_diffuse_src2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 2

    ! ----------------------------------
    bg = [1e-3, 1e-2, .0]
    S_target = [0.226478, 0.0919335, 0.211099, 0.211049, 0.226539]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src2_3')

    ! ----------------------------------
    bg = [1e-3, 1e-2, 1.0] ! in case of pure forward scattering, it should be the same as just absorption
    S_target = [2.41137e-01, 0.00000e+00, 2.42172e-01, 2.42070e-01, 2.41502e-01]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src2_2')

    ! ----------------------------------
    bg = [1e-3, 0., .0]
    S_target = [2.41137e-01, 0.00000e+00, 2.42172e-01, 2.42070e-01, 2.41502e-01]

    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src2_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_diffuse_src3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 3

    ! ----------------------------------
    bg = [1e-3, 1e-2, .0]
    S_target = [0.226478, 0.211099, 0.0919335, 0.211049, 0.226539]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src3_2')

    ! ----------------------------------
    bg = [1e-3, 0., .0]
    S_target = [2.41482e-01, 2.41252e-01, 0.00000e+00, 2.42897e-01, 2.41190e-01]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src3_1')
  end subroutine

  @test(npes=[1, 2])
  subroutine test_boxmc_select_cases_diffuse_src4(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: src = 4

    ! ----------------------------------
    bg = [1e-3, 1e-2, .0]
    S_target = [0.226478, 0.21119225, 0.21119225, 0.0919335, 0.226539]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src4_2')

    ! ----------------------------------
    bg = [1e-3, 0., .0]
    S_target = [0.241484, 0.242103, 0.242103, 0., 0.241141]
    call bmc_wedge_5_5%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_wedgemc_diffuse_src4_1')
  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

    real(ireal_dp), parameter :: sigma = 3 ! normal test range for coefficients

    character(len=*), optional :: msg
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals), parameter :: test_atol = real(atol, ireals) * real(sigma, ireals)

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

      @assertEqual(S_target, S, test_atol, local_msgS)
      @assertLessThanOrEqual(zero, S)
      @assertGreaterThanOrEqual(one, S)

      @assertEqual(T_target, T, test_atol, local_msgT)
      @assertLessThanOrEqual(zero, T)
      @assertGreaterThanOrEqual(one, T)
    end if
  end subroutine

end module
