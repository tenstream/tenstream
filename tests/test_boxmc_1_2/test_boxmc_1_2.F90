module test_boxmc_1_2
  use m_boxmc, only: t_boxmc, t_boxmc_1_2
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_helper_functions, only: toStr, cstr, colored_str_by_range, deg2rad
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry
  use m_eddington, only: eddington_coeff_zdun, eddington_coeff_bm, eddington_coeff_ec

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals), target :: S(2), T(1), S_target(2), T_target(1)
  real(ireals) :: S_tol(2), T_tol(1)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_1_2) :: bmc

  integer(mpiint) :: myid, mpierr, numnodes, comm
  character(len=default_str_len) :: msg

  real(ireal_dp), parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if (myid .eq. 0) print *, numnodes, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 0

    dx = 100
    dy = dx
    dz = 100

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
  subroutine test_boxmc_select_cases_diffuse_src_z(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), parameter :: c1 = 1 ! opposite side, same stream
    !real(ireals) :: tau, w0, g, mu0, a11, a12, a13, a23, a33

    bg = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp]

    theta = 0
    phi = 0
    T_target = zero

    ! should send diffuse radiation from bot face towards top face
    do src = 1, 2
      S_target = 0
      S_target(src) = c1

      call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_diffuse', src
      call check(S_target, T_target, S, T, msg=msg)
    end do

    bg = [1._ireal_dp / dz, 1._ireal_dp / dz, .8_ireal_dp]

    !tau = real(sum(bg(1:2))*dz, ireals)
    !w0  = real(bg(2) / sum(bg(1:2)), ireals)
    !g   = real(bg(3), ireals)
    !mu0 = real(cos(deg2rad(theta)), ireals)

    !call eddington_coeff_bm(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'bm  ', a11, a12, a13, a23, a33, ':', a33+a13+a23
    !call eddington_coeff_ec(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'ec  ', a11, a12, a13, a23, a33, ':', a33+a13+a23
    !call eddington_coeff_zdun(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'zdun', a11, a12, a13, a23, a33, ':', a33+a13+a23

    src = 1

    S_target = [real(ireals) :: 0.18391, 0.03514]
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src, 'norm', sum(S)
    call check(S_target, T_target, S, T, msg=msg)

    src = 2
    S_target = S_target(2:1:-1)
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src, 'norm', sum(S)
    call check(S_target, T_target, S, T, msg=msg)
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    !real(ireals) :: tau, w0, g, mu0, a11, a12, a13, a23, a33
    integer(iintegers) :: iphi

    src = 1
    bg = [0.1_ireal_dp / dz, 1._ireal_dp / dz, .8_ireal_dp]

    !tau = real(sum(bg(1:2))*dz, ireals)
    !w0  = real(bg(2) / sum(bg(1:2)), ireals)
    !g   = real(bg(3), ireals)
    !mu0 = real(cos(deg2rad(theta)), ireals)

    !call eddington_coeff_bm(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'bm', a11, a12, a13, a23, a33, ':', a33+a13+a23
    !call eddington_coeff_ec(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'ec', a11, a12, a13, a23, a33, ':', a33+a13+a23
    !call eddington_coeff_zdun(tau, w0, g, mu0, a11, a12, a13, a23, a33)
    !print *,'zdun', a11, a12, a13, a23, a33, ':', a33+a13+a23

    do iphi = 0, 5
      phi = iphi * 75

      theta = 60
      T_target = real(exp(-sum(bg(1:2)) * dz / cos(deg2rad(theta))), ireals)
      S_target = [real(ireals) :: 0.15668, 0.524605]
      call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta, 'norm', sum(T) + sum(S)
      call check(S_target, T_target, S, T, msg=msg)

      theta = 180 - 60
      S_target = S_target(2:1:-1)
      call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta, 'norm', sum(T) + sum(S)
      call check(S_target, T_target, S, T, msg=msg)
    end do

    do iphi = 0, 5
      phi = iphi * 75

      theta = 0
      T_target = real(exp(-sum(bg(1:2)) * dz / cos(deg2rad(theta))), ireals)
      S_target = [real(ireals) :: 0.04755, 0.50553]
      call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta, 'norm', sum(T) + sum(S)
      call check(S_target, T_target, S, T, msg=msg)

      theta = 180
      S_target = S_target(2:1:-1)
      call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta, 'norm', sum(T) + sum(S)
      call check(S_target, T_target, S, T, msg=msg)
    end do
  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), target, dimension(:) :: S_target, T_target, S, T
    character(len=*), optional :: msg

    real(ireals) :: S_diff(size(S))
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals), parameter :: test_atol = real(atol, ireals) * real(sigma, ireals)

    real(ireals), parameter :: color_limits(5) = [0., 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: colors(4) = [character(len=10) :: 'black', 'blue', 'green', 'red']

    real(ireals), parameter :: diff_color_limits(8) = [-1., -.1, -1e-3, -1e-5, 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: diff_colors(7) = [character(len=10) :: 'red', 'green', 'blue', 'black', 'blue', 'green', 'red']

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

      S_diff = S_target - S
      print *, '---------------------'
      write (*, FMT='(" diffuse :: "," :: ",A)') colored_str_by_range(S, color_limits, colors)
      !print *,cstr('---', 'blue')
      write (*, FMT='(" target  :: "," :: ",A)') colored_str_by_range(S_target, color_limits, colors)
      !print *,cstr('---', 'blue')
      write (*, FMT='(" diff    :: "," :: ",A)') colored_str_by_range(S_diff, diff_color_limits, diff_colors)
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
