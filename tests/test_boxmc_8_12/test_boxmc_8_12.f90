module test_boxmc_8_12
  use m_boxmc, only : t_boxmc, t_boxmc_8_12
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
  real(ireals) :: S(12),T(8), S_target(12), T_target(8)
  real(ireals) :: S_tol(12),T_tol(8)
  real(ireals), allocatable :: vertices(:)

  type(t_boxmc_8_12) :: bmc_8_12

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

    call bmc_8_12%init(comm)

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
    bg  = [1e-3_ireals, 1e-3_ireals, one ]
    phi = 0; src = 1

    ! should send diffuse radiation into vertical angle stream
    theta = 55
    T_target = [0.84242, 0.00000, 0.14042, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00848, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_12%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface_45',src
    call check(S_target,T_target, S,T, msg=msg)


    ! theta 65 should send diffuse radiation into low angle stream
    theta = 65
    T_target = [0.76688, 0.00000, 0.20966, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    S_target = [0.00000, 0.00000, 0.00000, 0.01170, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
    call bmc_8_12%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface_45',src
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
    S_target = [0.00001, 0.94849, 0.00001, 0.00015, 0.01285, 0.01287, 0.00000, 0.00000, 0.01278, 0.01283, 0.00000, 0.00000]
    call bmc_8_12%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')


    src = 4   ! top face -- low angle stream
    S_target = [0.00005, 0.00046, 0.00016, 0.78702, 0.05313, 0.05313, 0.00006, 0.00002, 0.05332, 0.05264, 0.00002, 0.00001]
    call bmc_8_12%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srcbotface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    ! direct to diffuse tests
    bg  = [1e-6_ireals, 1e-3_ireals, .99_ireals ]

    T_target = zero

    src = 1   ! bot face -- forward stream
    S_target = [0.94849, 0.00001, 0.00015, 0.00001, 0.00000, 0.00000, 0.01285, 0.01287, 0.00000, 0.00000, 0.01278, 0.01283]
    call bmc_8_12%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')


    src = 3   ! bot face -- low angle stream
    S_target = [0.00046, 0.00005, 0.78702, 0.00016, 0.00006, 0.00002, 0.05313, 0.05313, 0.00002, 0.00001, 0.05332, 0.05264]
    call bmc_8_12%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')
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
      write(*, FMT='( " diffuse ::  :: ",12(f9.5) )' ) S
      write(*, FMT='( " target  ::  :: ",12(f9.5) )' ) S_target
      write(*, FMT='( " diff    ::  :: ",12(f9.5) )' ) S_target-S
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
