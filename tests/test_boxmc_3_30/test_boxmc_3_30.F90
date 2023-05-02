module test_boxmc_3_30
  use m_boxmc, only: t_boxmc, t_boxmc_3_30
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_helper_functions, only: toStr, cstr, colored_str_by_range, deg2rad
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals), target :: S(30), T(3), S_target(30), T_target(3)
  real(ireals) :: S_tol(30), T_tol(3)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_3_30) :: bmc

  integer(mpiint) :: myid, mpierr, numnodes, comm
  character(len=120) :: msg

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
    real(ireals), parameter :: c1 = .0124 ! opposite side, same stream
    real(ireals), parameter :: c2 = .2170 ! left and adjacent sides center streams
    real(ireals), parameter :: c3 = .27625! left and adjacent sides upward stream

    bg = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp]

    theta = 0
    phi = 0
    T_target = zero

    ! should send diffuse radiation from bot face towards north west directions
    src = 3
    S_target = 0
    S_target(src) = c1
    S_target([12, 21]) = c2
    S_target([16, 23]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from top face towards north west directions
    src = 4
    S_target = 0
    S_target(src) = c1
    S_target([12, 21]) = c2
    S_target([20, 27]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from bot face towards north east directions
    src = 5
    S_target = 0
    S_target(src) = c1
    S_target([11, 21]) = c2
    S_target([15, 25]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from top face towards north east directions
    src = 6
    S_target = 0
    S_target(src) = c1
    S_target([11, 21]) = c2
    S_target([19, 29]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from bot face towards south west directions
    src = 7
    S_target = 0
    S_target(src) = c1
    S_target([12, 22]) = c2
    S_target([14, 24]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from top face towards south west directions
    src = 8
    S_target = 0
    S_target(src) = c1
    S_target([12, 22]) = c2
    S_target([18, 28]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from bot face towards south east directions
    src = 9
    S_target = 0
    S_target(src) = c1
    S_target([11, 22]) = c2
    S_target([13, 26]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation from top face towards south east directions
    src = 10
    S_target = 0
    S_target(src) = c1
    S_target([11, 22]) = c2
    S_target([17, 30]) = c3

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diffuse_mid(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), parameter :: c1 = 0.4407 ! straight
    real(ireals), parameter :: c2 = 0.0699 ! side

    bg = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp]

    theta = 0
    phi = 0
    T_target = zero

    ! should send diffuse radiation into vertical angle stream and a bit to the sides
    src = 1
    S_target = 0
    S_target(1) = c1
    S_target([13, 14, 15, 16, 23, 24, 25, 26]) = c2

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    src = 2
    S_target = 0
    S_target(2) = c1
    S_target([17, 18, 19, 20, 27, 28, 29, 30]) = c2

    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation into horizontal angle stream and a bit to the sides
    src = 11
    S_target = 0
    S_target(11) = c1
    S_target([5, 9, 6, 10, 25, 29, 26, 30]) = c2
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    src = 12
    S_target = 0
    S_target(12) = c1
    S_target([3, 7, 4, 8, 23, 27, 24, 28]) = c2
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    ! should send diffuse radiation into horizontal angle stream and a bit to the sides
    src = 21
    S_target = 0
    S_target(21) = c1
    S_target([3, 5, 4, 6, 15, 19, 16, 20]) = c2
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)

    src = 22
    S_target = 0
    S_target(22) = c1
    S_target([7, 9, 8, 10, 13, 17, 14, 18]) = c2
    call bmc%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_diffuse', src
    call check(S_target, T_target, S, T, msg=msg)
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    bg = [1._ireal_dp / dz, 1e-1_ireal_dp, 1._ireal_dp]

    ! should send diffuse radiation into vertical angle stream
    theta = 0
    phi = 0; src = 1
    T_target = [real(ireals) :: exp(real(-sum(bg(1:2)) * dz, ireals)), 0, 0]
    S_target = 0; S_target(2) = .367879
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta

    call check(S_target, T_target, S, T, msg=msg)

    theta = 180
    phi = 0; src = 1
    T_target = [real(ireals) :: exp(real(-sum(bg(1:2)) * dz, ireals)), 0, 0]
    S_target = 0; S_target(1) = .367879
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    ! photons to the north east
    theta = 45
    phi = 45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([6]) = 0.0208
    S_target([19, 29]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    theta = 180 - 45
    S_target = 0
    S_target([5]) = 0.0208
    S_target([15, 25]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    ! to the south east
    theta = 45
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    phi = 90 + 45; src = 1
    S_target = 0
    S_target([10]) = 0.0208
    S_target([17, 30]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    theta = 180 - 45
    phi = 90 + 45; src = 1
    S_target = 0
    S_target([9]) = 0.0208
    S_target([13, 26]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    ! to the south west
    theta = 45
    phi = 180 + 45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([8]) = 0.0208
    S_target([18, 28]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    theta = 180 - 45
    phi = 180 + 45; src = 1
    S_target = 0
    S_target([7]) = 0.0208
    S_target([14, 24]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    ! to the north west
    theta = 45
    phi = -45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([4]) = 0.0208
    S_target([20, 27]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)

    theta = 180 - 45
    phi = -45; src = 1
    S_target = 0
    S_target([3]) = 0.0208
    S_target([16, 23]) = 0.2318
    call bmc%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src, phi, theta
    call check(S_target, T_target, S, T, msg=msg)
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
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 1, colored_str_by_range(S(1:10:2), color_limits, colors)
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 2, colored_str_by_range(S(2:10:2), color_limits, colors)
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 3, colored_str_by_range(S(12:20:2), color_limits, colors)
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 4, colored_str_by_range(S(11:20:2), color_limits, colors)
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 5, colored_str_by_range(S(22:30:2), color_limits, colors)
      write (*, FMT='(" diffuse :: ",I0," :: ",A)') 6, colored_str_by_range(S(21:30:2), color_limits, colors)
      print *, cstr('---', 'blue')
      write (*, FMT='(" target  :: ",I0," :: ",A)') 1, colored_str_by_range(S_target(1:10:2), color_limits, colors)
      write (*, FMT='(" target  :: ",I0," :: ",A)') 2, colored_str_by_range(S_target(2:10:2), color_limits, colors)
      write (*, FMT='(" target  :: ",I0," :: ",A)') 3, colored_str_by_range(S_target(12:20:2), color_limits, colors)
      write (*, FMT='(" target  :: ",I0," :: ",A)') 4, colored_str_by_range(S_target(11:20:2), color_limits, colors)
      write (*, FMT='(" target  :: ",I0," :: ",A)') 5, colored_str_by_range(S_target(22:30:2), color_limits, colors)
      write (*, FMT='(" target  :: ",I0," :: ",A)') 6, colored_str_by_range(S_target(21:30:2), color_limits, colors)
      print *, cstr('---', 'blue')
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 1, colored_str_by_range(S_diff(1:10:2), diff_color_limits, diff_colors)
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 2, colored_str_by_range(S_diff(2:10:2), diff_color_limits, diff_colors)
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 3, colored_str_by_range(S_diff(12:20:2), diff_color_limits, diff_colors)
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 4, colored_str_by_range(S_diff(11:20:2), diff_color_limits, diff_colors)
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 5, colored_str_by_range(S_diff(22:30:2), diff_color_limits, diff_colors)
      write (*, FMT='(" diff    :: ",I0," :: ",A)') 6, colored_str_by_range(S_diff(21:30:2), diff_color_limits, diff_colors)
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
