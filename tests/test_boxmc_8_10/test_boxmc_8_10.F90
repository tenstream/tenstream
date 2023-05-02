module test_boxmc_8_10
  use m_boxmc, only: t_boxmc, t_boxmc_8_10
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(10), T(8), S_target(10), T_target(8)
  real(ireals) :: S_tol(10), T_tol(8)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_8_10) :: bmc_8_10

  integer(mpiint) :: myid, mpierr, numnodes, comm
  character(len=120) :: msg

  real(ireal_dp), parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
  !real(ireal_dp),parameter :: atol=1e-5, rtol=1e-3
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_8_10%init(comm)

    if (myid .eq. 0) print *, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 0

    dx = 100
    dy = dx
    dz = 50

    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1])
  subroutine test_boxmc_symmetry_in_phi(this)   ! Check that we have symmetry for total transmission for e.g. phi 0==90 or 10==80 etc.
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi, itheta
    real(ireal_dp) :: tau, Tsum(2)

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! from top to bot face
    tau = (bg(1) + bg(2)) * dz

    do iphi = 0, 45, 20
      do itheta = 0, 85, 40
        phi = real(iphi, ireal_dp)
        theta = real(itheta, ireal_dp)
        Tsum = 0._ireal_dp

        do src = 1, 8
          call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
          Tsum(1) = Tsum(1) + sum(T)
        end do

        phi = 90._ireal_dp - real(iphi, ireal_dp)
        do src = 1, 8
          call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
          Tsum(2) = Tsum(2) + sum(T)
        end do

        write (msg, *) 'test_boxmc_symmetry_in_phi : ', iphi, phi, theta
        @assertEqual(Tsum(1), Tsum(2), atol * sigma, msg)
        print *, msg, Tsum
      end do
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! from top to bot face
    theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero

    do iphi = 0, 90, 10
      phi = real(iphi, ireal_dp)
      do src = 1, 4
        T_target = zero
        T_target(src) = real(exp(-tau), ireals)
        call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        write (msg, *) 'test_boxmc_select_cases_direct_srctopface top_to_bot', src, ':: phi', phi
        call check(S_target, T_target, S, T, msg='')
      end do
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface_45(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! outwards from face 2
    phi = 0; theta = 45 * 1._ireal_dp

    tau = (bg(1) + bg(2)) * dz * sqrt(2 * 1._ireal_dp)

    S_target = zero

    do src = 1, 2
      T_target = zero
      T_target(2 + src) = real(exp(-tau), ireals)
      print *, 'phi,theta test', phi, theta
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      print *, 'phi,theta test2', phi, theta
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface_45', src
      call check(S_target, T_target, S, T, msg=msg)
    end do

    phi = 90
    do src = 1, 3, 2
      T_target = zero
      T_target(src + 1) = real(exp(-tau), ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface_45', src
      call check(S_target, T_target, S, T, msg=msg)
    end do

  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srcsidefaces(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! along each of the side faces
    phi = 0; theta = 90

    tau = (bg(1) + bg(2)) * dx

    S_target = zero
    T_target = zero

    do src = 5, 6, 1
      phi = 0

      T_target = zero

      T_target(src + 2) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')

      phi = 90

      T_target = zero

      T_target(src) = real(exp(-tau), ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')
    end do

    do src = 7, 8
      phi = 0

      T_target = zero

      T_target(src) = real(exp(-tau), ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')

      phi = 90

      T_target = zero

      T_target(src - 2) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')
    end do

    do iphi = 0, 360, 60
      phi = real(iphi, ireal_dp)
      theta = 0

      tau = (bg(1) + bg(2)) * dz / 2

      src = 5
      T_target = zero
      T_target([1, 3]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 6
      T_target = zero
      T_target([1, 3]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2 * exp(-tau), ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 7
      T_target = zero
      T_target([1, 2]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2, ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 8
      T_target = zero
      T_target([1, 2]) = real((sinh(tau) - cosh(tau) + 1) / tau / 2 * exp(-tau), ireals)
      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = [0.0, 0.39, 0.1404, 0.1404, 0.0, 0.0, 0.1404, 0.1404, 0.0, 0.0]

    T_target = zero

    src = 2   !top face
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srcbottomface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    T_target = zero

    S_target = [0.39, 0.0, 0.0, 0.0, 0.1404, 0.1404, 0.0, 0.0, 0.1404, 0.1404]

    src = 1
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbottomface')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srcsideface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau
    real(ireals), parameter :: top = 0.56164175, side1 = 0.1047592, side2 = 0.14247725

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz
    T_target = zero

    src = 3
    S_target = [zero, top, side1, zero, zero, zero, side2, side2, zero, zero]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 3')

    src = 4
    S_target = [zero, top, zero, side1, zero, zero, side2, side2, zero, zero]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 4')

    src = 5
    S_target = [top, zero, zero, zero, side1, zero, zero, zero, side2, side2]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 5')

    src = 6
    S_target = [top, zero, zero, zero, zero, side1, zero, zero, side2, side2]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 6')

    src = 7
    S_target = [zero, top, side2, side2, zero, zero, side1, zero, zero, zero]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 7')

    src = 8
    S_target = [zero, top, side2, side2, zero, zero, zero, side1, zero, zero]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 8')

    src = 9
    S_target = [top, zero, zero, zero, side2, side2, zero, zero, side1, zero]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 9')

    src = 10
    S_target = [top, zero, zero, zero, side2, side2, zero, zero, zero, side1]
    call bmc_8_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 10')

  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

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
