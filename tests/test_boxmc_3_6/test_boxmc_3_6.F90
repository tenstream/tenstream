module test_boxmc_3_6
  use m_boxmc, only: t_boxmc, t_boxmc_3_6
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(6), T(3), S_target(6), T_target(3)
  real(ireals) :: S_tol(6), T_tol(3)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_3_6) :: bmc_3_6

  integer(mpiint) :: myid, mpierr, numnodes, comm

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
  ! real(ireal_dp),parameter :: atol=1e-4, rtol=1e-3
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_3_6%init(comm)

    if (myid .eq. 0) print *, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 45

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
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! from top to bot face
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_direct_srctopface top_to_bot')

  end subroutine

  !@test(npes =[1,2])
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

    T_target = zero
    T_target(1) = real(exp(-tau) / 2, ireals)
    T_target(3) = real((1 - exp(-tau)) / (2 * tau), ireals)

    src = 1

    call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_direct_srctopface_45')

  end subroutine

  !@test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_srcsidefaces(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! along each of the side faces
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    do iphi = 0, 360, 30
      phi = real(iphi, ireal_dp)
      do src = 2, 3
        call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')
      end do
    end do
  end subroutine

  !@test(npes =[1,2])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero
    S_target = [0., 0.242562, 0.177445, 0.177294, 0.177382, 0.177347]

    T_target = zero

    do src = 2, 2
      call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface')
    end do
  end subroutine

  !@test(npes =[1])
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

    S_target = [0.242562, 0., 0.177445, 0.177294, 0.177382, 0.177347]

    src = 1
    call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbottomface')

    S_target = [0., 0.242562, 0.177445, 0.177294, 0.177382, 0.177347]
    !src = 1
    !call bmc_3_6%get_coeff(comm,bg,src,.False.,phi,theta,1e3_ireals, 1e3_ireals, 1e2_ireals,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    !call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbottomface')
  end subroutine

  !@test(npes =[1,2])
  subroutine test_boxmc_select_cases_diff_srcsideface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz
    T_target = zero

    src = 4
    S_target = [0.296326, 0.296192, 0.0, 0.056321, 0.154418, 0.154580]
    call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface')

    src = 3
    S_target = [0.296326, 0.296192, 0.056321, 0.0, 0.154418, 0.154580]
    call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface')

    src = 6
    S_target = [0.296326, 0.296192, 0.154418, 0.154580, 0.0, 0.056321]
    call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface')

    src = 5
    S_target = [0.296326, 0.296192, 0.154418, 0.154580, 0.056321, 0.0]
    call bmc_3_6%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface')
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
