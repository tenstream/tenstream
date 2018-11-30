module test_boxmc_8_16
  use m_boxmc, only : t_boxmc, t_boxmc_8_16
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_helper_functions, only: deg2rad

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(16),T(8), S_target(16), T_target(8)
  real(ireals) :: S_tol(16),T_tol(8)
  real(ireals), allocatable :: vertices(:)

  type(t_boxmc_8_16) :: bmc_8_16

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireals),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireals),parameter :: atol=1e-3, rtol=1e-2
  !real(ireals),parameter :: atol=1e-5, rtol=1e-3
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_8_16%init(comm)

    if(myid.eq.0) print *,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 5

    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-3_ireals, 1e-1_ireals, one ]
    phi = 0; src = 1

    ! should send diffuse radiation into vertical angle stream
    theta = 20
    T_target = [0.56285, 0.00000, 0.02126, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.41059, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src
    call check(S_target,T_target, S,T, msg=msg)


    ! send diffuse radiation into mid angle stream
    theta = 50
    T_target = [0.40155, 0.00000, 0.05429, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.53641, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation into lower angle stream
    theta = 70
    T_target = [0.16569, 0.00000, 0.06262, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.75717, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src
    call check(S_target,T_target, S,T, msg=msg)

    ! theta 85 should send diffuse radiation into lowest angle stream
    theta = 85
    T_target = [0.00000, 0.00000, 0.00263, 0.00000, 0.00000, 0.00000, 0.00063, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.80653, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.13494, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-6_ireals, 1e-3_ireals, .99_ireals ]

    T_target = zero

    src = 2   ! top face -- forward stream
    S_target = [0.00000, 0.96704, 0.00000, 0.00018, 0.00000, 0.00002, 0.00000, 0.00001, 0.00819, 0.00818, 0.00000, 0.00000, 0.00819, 0.00817, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_vertical_angle')


    src = 4   ! top face -- low angle stream
    S_target = [0.00000, 0.00026, 0.00000, 0.92175, 0.00001, 0.00029, 0.00001, 0.00002, 0.01940, 0.01942, 0.00000, 0.00000, 0.01942, 0.01941, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_mid_angle')


    src = 6   ! top face -- low angle stream
    S_target = [0.00001, 0.00004, 0.00001, 0.00049, 0.00002, 0.84728, 0.00003, 0.00043, 0.03789, 0.03795, 0.00000, 0.00001, 0.03793, 0.03788, 0.00000, 0.00001]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_low_angle')


    src = 8   ! top face -- lowest angle stream
    S_target = [0.00002, 0.00004, 0.00004, 0.00010, 0.00009, 0.00128, 0.00028, 0.59641, 0.10037, 0.10037, 0.00008, 0.00007, 0.10033, 0.10032, 0.00008, 0.00008]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_lowest_angle')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srcbotface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-6_ireals, 1e-3_ireals, .99_ireals ]

    T_target = zero

    src = 1   ! bot face -- forward stream
    S_target = [0.96703, 0.00000, 0.00018, 0.00000, 0.00002, 0.00000, 0.00001, 0.00000, 0.00000, 0.00000, 0.00818, 0.00819, 0.00000, 0.00000, 0.00819, 0.00819]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_vert_angle')

    src = 3   ! bot face -- mid stream
    S_target = [0.00026, 0.00000, 0.92172, 0.00000, 0.00029, 0.00001, 0.00002, 0.00001, 0.00000, 0.00000, 0.01941, 0.01941, 0.00000, 0.00000, 0.01944, 0.01942]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_mid_angle')

    src = 5   ! bot face -- lower stream
    S_target = [0.00004, 0.00001, 0.00049, 0.00001, 0.84715, 0.00002, 0.00044, 0.00003, 0.00000, 0.00001, 0.03803, 0.03792, 0.00001, 0.00001, 0.03787, 0.03796]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_low_angle')

    src = 7   ! bot face -- lowest stream
    S_target = [0.00004, 0.00002, 0.00010, 0.00004, 0.00129, 0.00008, 0.59656, 0.00028, 0.00008, 0.00008, 0.10031, 0.10032, 0.00008, 0.00008, 0.10031, 0.10030]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_lowest_angle')
  end subroutine

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

    character(len=*),optional :: msg
    character(default_str_len) :: local_msgS, local_msgT

    if(myid.eq.0) then
      print*,''

      if( present(msg) ) then
        write(local_msgS,*) trim(msg),':: Diffuse boxmc coefficient not as '
        write(local_msgT,*) trim(msg),':: Direct  boxmc coefficient not as '
        print *,msg
      else
        write(local_msgS,*) 'Diffuse boxmc coefficient not as '
        write(local_msgT,*) 'Direct  boxmc coefficient not as '
      endif

      print*,'---------------------'
      write(*, FMT='( " diffuse ::  :: ",16(f9.5) )' ) S
      write(*, FMT='( " target  ::  :: ",16(f9.5) )' ) S_target
      write(*, FMT='( " diff    ::  :: ",16(f9.5) )' ) S_target-S
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(f9.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(f9.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(f9.5) )' ) T_target-T
      print*,'---------------------'
      print*,''

      @assertEqual(S_target, S, atol*sigma, local_msgS )
      @assertLessThanOrEqual   (zero, S)
      @assertGreaterThanOrEqual(one , S)

      @assertEqual(T_target, T, atol*sigma, local_msgT )
      @assertLessThanOrEqual   (zero, T)
      @assertGreaterThanOrEqual(one , T)
    endif
  end subroutine

end module
