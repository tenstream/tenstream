module test_boxmc_8_16
  use m_boxmc, only : t_boxmc, t_boxmc_8_16
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_helper_functions, only: deg2rad, itoa, ftoa

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

    theta = 60
    phi = 0; src=1
    T_target = [0.30206, 0.00000, 0.06325, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.62473, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface'//itoa(src)//' phi='//ftoa(phi)
    call check(S_target,T_target, S,T, msg=msg)

    phi = 270; src=2
    T_target = [0.06317, 0.30127, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.62561, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface'//itoa(src)//' phi='//ftoa(phi)
    call check(S_target,T_target, S,T, msg=msg)

    phi = 180; src=3
    T_target = [0.06300, 0.00000, 0.30112, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.62593, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface'//itoa(src)//' phi='//ftoa(phi)
    call check(S_target,T_target, S,T, msg=msg)

    phi = 90; src=3
    T_target = [0.00000, 0.00000, 0.30130, 0.06258, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.62618, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface'//itoa(src)//' phi='//ftoa(phi)
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-6_ireals, 1e-3_ireals, .99_ireals ]

    T_target = zero

    src = 2
    S_target = [0.00004, 0.90780, 0.00001, 0.00013, 0.00001, 0.00001, 0.00001, 0.00012, 0.01306, 0.01308, 0.00001, 0.00000, 0.00000, 0.06572, 0.00000, 0.00001]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_quadrant_north')

    src = 4
    S_target = [0.00001, 0.00013, 0.00003, 0.90747, 0.00001, 0.00012, 0.00000, 0.00001, 0.06606, 0.00000, 0.00002, 0.00000, 0.01299, 0.01313, 0.00000, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_quadrant_west')

    src = 6
    S_target = [0.00000, 0.00002, 0.00000, 0.00013, 0.00003, 0.90742, 0.00001, 0.00012, 0.01294, 0.01324, 0.00001, 0.00000, 0.06604, 0.00000, 0.00002, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_quadrant_south')

    src = 8
    S_target = [0.00001, 0.00012, 0.00000, 0.00001, 0.00001, 0.00011, 0.00005, 0.90754, 0.00000, 0.06622, 0.00000, 0.00002, 0.01317, 0.01272, 0.00000, 0.00001]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface_quadrant_east')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srcbotface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-6_ireals, 1e-3_ireals, .99_ireals ]

    T_target = zero

    src = 1
    S_target = [0.90792, 0.00003, 0.00013, 0.00001, 0.00002, 0.00000, 0.00009, 0.00001, 0.00000, 0.00000, 0.01283, 0.01301, 0.00000, 0.00002, 0.00000, 0.06591]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_quadrant_north')

    src = 3
    S_target = [0.00013, 0.00001, 0.90763, 0.00003, 0.00013, 0.00001, 0.00000, 0.00000, 0.00002, 0.00000, 0.06609, 0.00000, 0.00000, 0.00000, 0.01290, 0.01304]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_quadrant_west')

    src = 5
    S_target = [0.00001, 0.00001, 0.00013, 0.00000, 0.90761, 0.00004, 0.00012, 0.00001, 0.00001, 0.00001, 0.01296, 0.01327, 0.00001, 0.00000, 0.06581, 0.00000]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_quadrant_south')

    src = 7
    S_target = [0.00013, 0.00001, 0.00001, 0.00000, 0.00011, 0.00001, 0.90678, 0.00003, 0.00000, 0.00002, 0.00000, 0.06684, 0.00000, 0.00000, 0.01284, 0.01321]
    call bmc_8_16%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbotface_quadrant_east')

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
