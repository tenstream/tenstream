module test_boxmc
  use m_boxmc, only: t_boxmc_8_10
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
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_8_10) :: bmc_8_10

  integer(mpiint) :: myid, mpierr, numnodes, comm

  real(ireal_dp), parameter :: atol = 1e-2, rtol = 1e-1
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
    theta = 45

    dx = 100
    dy = dx
    dz = 50

    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    ! computed target with stddev_atol=5e-6, stddev_rtol=1e-4 in optprop_parameters
    ! inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) !
    !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,vertices,S,T,S_tol,T_tol,inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct(this)
    class(MpiTestMethod), intent(inout) :: this

    ! direct tests
    bg = [1e-3, 1e-3, 0.]
    phi = 0; theta = 0
    S_target = [1.16098e-02, 1.13492e-02, 4.69646e-03, 1.11393e-03, 4.61549e-03,&
      & 1.08210e-03, 4.70091e-03, 1.11544e-03, 4.61404e-03, 1.08313e-03]
    T_target = [9.04842e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00]

    call bmc_8_10%get_coeff(comm, bg, i1, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_direct_lambert_beer(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: src

    ! direct tests
    bg = [1e-2, 0., 0.]
    phi = 0; theta = 0
    S_target = zero

    do src = 1, 4
      T_target = zero
      T_target(src) = real(exp(-(bg(1) + bg(2)) * dz), ireals)

      call bmc_8_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target, T_target, S, T, msg='test_boxmc_direct_lambert_beer')
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diffuse(this)
    class(MpiTestMethod), intent(inout) :: this

    bg = [1e-3, 1e-16, .0]
    S_target = [0.00000e+00, 3.90217e-01, 1.40431e-01, 1.40417e-01, 0.00000e+00,&
      & 0.00000e+00, 1.40416e-01, 1.40429e-01, 0.00000e+00, 0.00000e+00]
    T_target = zero

    call bmc_8_10%get_coeff(comm, bg, i1, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_diffuse')
! --------------------------------------------

    bg = [1e-3, 1e-4, .5]
    S_target = [5.91803e-04, 3.89499e-01, 1.40297e-01, 1.40317e-01, 1.44767e-04,&
      & 1.45072e-04, 1.40312e-01, 1.40306e-01, 1.45083e-04, 1.44973e-04]

    call bmc_8_10%get_coeff(comm, bg, i1, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_diffuse')
! --------------------------------------------

    bg = [1e-3, 1e-4, .999999]
    S_target = [1.85844e-09, 3.90217e-01, 1.40431e-01, 1.40429e-01, 0.00000e+00,&
      & 0.00000e+00, 1.40417e-01, 1.40417e-01, 0.00000e+00, 9.64722e-10]

    call bmc_8_10%get_coeff(comm, bg, i1, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_diffuse')
  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

    real(ireals), parameter :: sigma = 3 ! normal test range for coefficients

    character(len=*), optional :: msg
    character(default_str_len) :: local_msgS, local_msgT

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

      @assertEqual(S_target, S, real(atol, ireals) * sigma, local_msgS)
      @assertLessThanOrEqual(zero, S)
      @assertGreaterThanOrEqual(one, S)

      @assertEqual(T_target, T, real(atol, ireals) * sigma, local_msgT)
      @assertLessThanOrEqual(zero, T)
      @assertGreaterThanOrEqual(one, T)
    end if
  end subroutine

end module
