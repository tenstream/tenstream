module test_boxmc_8_18
  use m_boxmc, only: t_boxmc, t_boxmc_8_18
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry
  use m_helper_functions, only: deg2rad, itoa, ftoa

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(18), T(8), S_target(18), T_target(8)
  real(ireals) :: S_tol(18), T_tol(8)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_8_18) :: bmc_8_18

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

    call bmc_8_18%init(comm)

    if (myid .eq. 0) print *, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 0

    dx = 100
    dy = dx
    dz = 5

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
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 1e-1_ireal_dp, 1._ireal_dp]

    ! should send diffuse radiation into vertical angle stream
    theta = 20
    phi = 0; src = 1
    T_target = [0.56289, 0.00000, 0.02127, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.41053, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
      & 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface', src
    call check(S_target, T_target, S, T, msg=msg)

    ! send diffuse radiation into lower angle stream, we have four of those, one for each quadrant
    theta = 60
    phi = 45; src = 1
    T_target = [0.28057, 0.03912, 0.03911, 0.00545, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.62579, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
      & 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface src='//itoa(src)//' phi='//ftoa(phi)
    call check(S_target, T_target, S, T, msg=msg)

    phi = 270 + 45; src = 2
    T_target = [0.03912, 0.28057, 0.00545, 0.03911, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.62588, 0.00000, 0.00000, 0.00000, 0.00000, &
      & 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface src='//itoa(src)//' phi='//ftoa(phi)
    call check(S_target, T_target, S, T, msg=msg)

    phi = 180 + 45; src = 4
    T_target = [0.00545, 0.03912, 0.03912, 0.28057, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.62588, 0.00000, 0.00000, &
      & 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface src='//itoa(src)//' phi='//ftoa(phi)
    call check(S_target, T_target, S, T, msg=msg)

    phi = 90 + 45; src = 3
    T_target = [0.03912, 0.00545, 0.28036, 0.03917, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.62587, &
      & 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    write (msg, *) ' test_boxmc_select_cases_direct_srctopface src='//itoa(src)//' phi='//ftoa(phi)
    call check(S_target, T_target, S, T, msg=msg)
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg = [1e-6_ireal_dp, 1e-3_ireal_dp, .99_ireal_dp]

    T_target = zero

    src = 2   ! top face -- forward stream
    S_target = [0.00000, 0.96691, 0.00000, 0.00006, 0.00000, 0.00005, 0.00000, 0.00006, 0.00000, 0.00005, &
      & 0.00817, 0.00827, 0.00000, 0.00000, 0.00819, 0.00824, 0.00000, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface_vertical_angle')

    src = 4   ! top face -- low angle stream upper right quadrant
    S_target = [0.00001, 0.00016, 0.00005, 0.86118, 0.00001, 0.00010, 0.00000, 0.00000, 0.00001, 0.00010, &
      & 0.01920, 0.01920, 0.00000, 0.00002, 0.00000, 0.09960, 0.00000, 0.00002]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface_low_angle_quadrant_1')

    src = 6   ! top face -- low angle stream upper left quadrant
    S_target = [0.00001, 0.00016, 0.00001, 0.00010, 0.00005, 0.86114, 0.00001, 0.00010, 0.00000, 0.00000, &
      & 0.09960, 0.00000, 0.00002, 0.00000, 0.01920, 0.01920, 0.00000, 0.00002]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface_low_angle_quadrant_2')

    src = 8   ! top face -- low angle stream lower left quadrant
    s_target = [0.00001, 0.00016, 0.00000, 0.00000, 0.00001, 0.00010, 0.00006, 0.86123, 0.00001, 0.00009, &
      & 0.01920, 0.01920, 0.00002, 0.00000, 0.09960, 0.00000, 0.00002, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, s, t, s_tol, t_tol, inp_atol=atol, inp_rtol=rtol)
    call check(s_target, t_target, s, t, msg=' test_boxmc_select_cases_diff_srctopface_low_angle_quadrant_3')

    src = 10   ! top face -- low angle stream lower right quadrant
    s_target = [0.00001, 0.00016, 0.00001, 0.00010, 0.00000, 0.00000, 0.00001, 0.00010, 0.00005, 0.86115, &
      & 0.00000, 0.09960, 0.00000, 0.00002, 0.01920, 0.01920, 0.00002, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, s, t, s_tol, t_tol, inp_atol=atol, inp_rtol=rtol)
    call check(s_target, t_target, s, t, msg=' test_boxmc_select_cases_diff_srctopface_low_angle_quadrant_4')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srcbotface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg = [1e-6_ireal_dp, 1e-3_ireal_dp, .99_ireal_dp]

    T_target = zero

    src = 1   ! bot face -- vertical stream
    S_target = [0.96691, 0.00000, 0.00000, 0.00006, 0.00000, 0.00005, 0.00000, 0.00006, 0.00000, 0.00005, &
      & 0.00000, 0.00000, 0.00817, 0.00827, 0.00000, 0.00000, 0.00819, 0.00824]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbotface_vertical_angle')

    src = 3   ! bot face -- low angle stream upper quadrant
    S_target = [0.00021, 0.00002, 0.86106, 0.00004, 0.00010, 0.00001, 0.00001, 0.00000, 0.00010, 0.00001, &
      & 0.00001, 0.00001, 0.01926, 0.01926, 0.00000, 0.00004, 0.00000, 0.09911]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbotface_low_angle_quadrant_1')

    src = 5   ! bot face -- low angle stream left quadrant
    S_target = [0.00017, 0.00001, 0.00014, 0.00000, 0.86051, 0.00004, 0.00013, 0.00001, 0.00001, 0.00000, &
      & 0.00001, 0.00000, 0.09911, 0.00000, 0.00000, 0.00001, 0.01920, 0.01920]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbotface_low_angle_quadrant_2')

    src = 7   ! bot face -- low angle stream bot quadrant
    s_target = [0.00014, 0.00001, 0.00001, 0.00000, 0.00010, 0.00001, 0.86030, 0.00003, 0.00013, 0.00001, &
      & 0.00003, 0.00001, 0.01920, 0.01920, 0.00002, 0.00000, 0.09911, 0.00000]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, s, t, s_tol, t_tol, inp_atol=atol, inp_rtol=rtol)
    call check(s_target, t_target, s, t, msg=' test_boxmc_select_cases_diff_srcbotface_low_angle_quadrant_3')

    src = 9   ! bot face -- low angle stream right quadrant
    s_target = [0.00020, 0.00001, 0.00011, 0.00002, 0.00001, 0.00001, 0.00007, 0.00001, 0.86021, 0.00007, &
      & 0.00000, 0.00001, 0.00001, 0.09911, 0.00004, 0.00000, 0.01920, 0.01920]
    call bmc_8_18%get_coeff(comm, bg, src, .false., phi, theta, vertices, s, t, s_tol, t_tol, inp_atol=atol, inp_rtol=rtol)
    call check(s_target, t_target, s, t, msg=' test_boxmc_select_cases_diff_srcbotface_low_angle_quadrant_4')
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
      write (*, FMT='( " diffuse ::  :: ",18(f9.5) )') S
      write (*, FMT='( " target  ::  :: ",18(f9.5) )') S_target
      write (*, FMT='( " diff    ::  :: ",18(f9.5) )') S_target - S
      print *, ''
      write (*, FMT='( " direct  ::  :: ", 8(f9.5) )') T
      write (*, FMT='( " target  ::  :: ", 8(f9.5) )') T_target
      write (*, FMT='( " diff    ::  :: ", 8(f9.5) )') T_target - T
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
