module test_boxmc_3_6_tau_scaling
  use m_boxmc, only: t_boxmc, t_boxmc_3_6
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry
  use m_helper_functions, only: itoa

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(6), T(3), S_target(6), T_target(3)
  real(ireals) :: S_tol(6), T_tol(3)
  real(ireal_dp) :: vertices(24)

  type(t_boxmc_3_6) :: bmc_3_6

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

    call bmc_3_6%init(comm)

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
  subroutine test_boxmc_tau_scaling3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! from top to bot face
    phi = 0; theta = 0

    ! direct to diffuse tests
    bg = [1e-4_ireal_dp, 1e-1_ireal_dp, 0._ireal_dp]

    tau = (bg(1) + bg(2)) * dz

    S_target = [4.99310e-01, 1.03362e-01, 9.60417e-02, 9.60534e-02, 9.60340e-02, 9.60392e-02]

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, &
                           inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_tau_scaling_3')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_tau_scaling1(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! from top to bot face
    phi = 0; theta = 0

    ! direct to diffuse tests
    bg = [1e-9_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp]

    tau = (bg(1) + bg(2)) * dz

    S_target(1:2) = [(1.21608e-02 + 1.21106e-02) / 2, 1.207e-02]
    S_target(3:6) = (6.10860e-03 + 6.07660e-03 + 6.12260e-03 + 6.14420e-03) / 4

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, &
                           inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_tau_scaling_1')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_tau_scaling2(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! from top to bot face
    phi = 0; theta = 0

    ! direct to diffuse tests
    bg = [1e-30_ireal_dp, 1e-5_ireal_dp / dz, 0._ireal_dp]

    tau = (bg(1) + bg(2)) * dz

    S_target(1:2) = (2.56000e-06 + 2.39000e-06 + 2.47300e-06 + 2.50900e-06) / 4
    S_target(3:6) = (1.26000e-06 + 1.25000e-06 + 1.38000e-06 + 1.32000e-06 &
                     + 1.24400e-06 + 1.23000e-06 + 1.28000e-06 + 1.26000e-06) / 8

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_6%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, &
                           inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_tau_scaling_2')
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

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
      write (*, FMT='( " diffuse ::  :: ",6(es12.5) )') S
      write (*, FMT='( " target  ::  :: ",6(es12.5) )') S_target
      write (*, FMT='( " diff    ::  :: ",6(es12.5) )') S_target - S
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
