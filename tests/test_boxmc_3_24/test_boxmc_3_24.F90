module test_boxmc_3_24
  use m_boxmc, only : t_boxmc, t_boxmc_3_24
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_helper_functions, only: toStr, cstr, colored_str_by_range, deg2rad
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals), target :: S(24),T(3), S_target(24), T_target(3)
  real(ireals) :: S_tol(24),T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_24) :: bmc

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireal_dp),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-2
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()


    call init_mpi_data_parameters(comm)

    call bmc%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if(myid.eq.0) print *,numnodes,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 100

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
  subroutine test_boxmc_select_cases_diffuse_src_z(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), parameter :: c1=.2 ! opposite side, same stream
    real(ireals), parameter :: c2=.4 ! left and adjacent sides upward stream

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from bot face towards north west directions
    src = 1
    S_target = 0
    S_target(src) = c1
    S_target([12,17]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from top face towards north west directions
    src = 2
    S_target = 0
    S_target(src) = c1
    S_target([16,21]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from bot face towards north east directions
    src = 3
    S_target = 0
    S_target(src) = c1
    S_target([11,19]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from top face towards north east directions
    src = 4
    S_target = 0
    S_target(src) = c1
    S_target([15,23]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from bot face towards south west directions
    src = 5
    S_target = 0
    S_target(src) = c1
    S_target([10,18]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from top face towards south west directions
    src = 6
    S_target = 0
    S_target(src) = c1
    S_target([14,22]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from bot face towards south east directions
    src = 7
    S_target = 0
    S_target(src) = c1
    S_target([ 9,20]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from top face towards south east directions
    src = 8
    S_target = 0
    S_target(src) = c1
    S_target([13,24]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src_x(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), parameter :: c1=.2 ! opposite side, same stream
    real(ireals), parameter :: c2=.4 ! left and adjacent sides upward stream

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    src = 9
    S_target = 0
    S_target(src) = c1
    S_target([7,20]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 10
    S_target = 0
    S_target(src) = c1
    S_target([5,18]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 11
    S_target = 0
    S_target(src) = c1
    S_target([3,19]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 12
    S_target = 0
    S_target(src) = c1
    S_target([1,17]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 13
    S_target = 0
    S_target(src) = c1
    S_target([8,24]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 14
    S_target = 0
    S_target(src) = c1
    S_target([6,22]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 15
    S_target = 0
    S_target(src) = c1
    S_target([ 4,23]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 16
    S_target = 0
    S_target(src) = c1
    S_target([ 2,21]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src_y(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), parameter :: c1=.2 ! opposite side, same stream
    real(ireals), parameter :: c2=.4 ! left and adjacent sides upward stream

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    src = 17
    S_target = 0
    S_target(src) = c1
    S_target([1,12]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation from top face towards north west directions
    src = 18
    S_target = 0
    S_target(src) = c1
    S_target([5,10]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 19
    S_target = 0
    S_target(src) = c1
    S_target([3,11]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 20
    S_target = 0
    S_target(src) = c1
    S_target([7,9]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 21
    S_target = 0
    S_target(src) = c1
    S_target([2,16]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 22
    S_target = 0
    S_target(src) = c1
    S_target([6,14]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 23
    S_target = 0
    S_target(src) = c1
    S_target([ 4,15]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 24
    S_target = 0
    S_target(src) = c1
    S_target([ 8,13]) = c2

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src

    bg  = [1._ireal_dp/dz, 1e-1_ireal_dp, 1._ireal_dp ]

    ! photons to the north east
    theta = 45
    phi = 45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([4]) = 0.0208
    S_target([15,23]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    S_target = 0
    S_target([3]) = 0.0208
    S_target([11,19]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! to the south east
    theta = 45
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    phi = 90+45; src = 1
    S_target = 0
    S_target([8]) = 0.0208
    S_target([13,24]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 90+45; src = 1
    S_target = 0
    S_target([7]) = 0.0208
    S_target([9,20]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! to the south west
    theta = 45
    phi = 180+45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([6]) = 0.0208
    S_target([14,22]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 180+45; src = 1
    S_target = 0
    S_target([5]) = 0.0208
    S_target([10,18]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! to the north west
    theta = 45
    phi = -45; src = 1
    T_target = [real(ireals) :: 0, 0.0434, 0.0434]
    S_target = 0
    S_target([2]) = 0.0208
    S_target([16,21]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = -45; src = 1
    S_target = 0
    S_target([1]) = 0.0208
    S_target([12,17]) = 0.2318
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),target,dimension(:) :: S_target,T_target, S,T
    character(len=*),optional :: msg

    real(ireals) :: S_diff(size(S))
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals), parameter :: test_atol = real(atol, ireals) * real(sigma, ireals)

    real(ireals), parameter :: color_limits(5) = [0., 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: colors(4) = [character(len=10):: 'black', 'blue', 'green', 'red']

    real(ireals), parameter :: diff_color_limits(8) = [-1., -.1, -1e-3, -1e-5, 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: diff_colors(7) = [character(len=10)::'red', 'green', 'blue', 'black', 'blue', 'green', 'red']

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

      S_diff = S_target - S
      print*,'---------------------'
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 1, colored_str_by_range(S( 1: 8:2), color_limits, colors)
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 2, colored_str_by_range(S( 2: 8:2), color_limits, colors)
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 3, colored_str_by_range(S(10:16:2), color_limits, colors)
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 4, colored_str_by_range(S( 9:16:2), color_limits, colors)
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 5, colored_str_by_range(S(18:24:2), color_limits, colors)
      write(*,FMT='(" diffuse :: ",I0," :: ",A)' ) 6, colored_str_by_range(S(17:24:2), color_limits, colors)
      print *,cstr('---', 'blue')
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 1, colored_str_by_range(S_target( 1: 8:2), color_limits, colors)
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 2, colored_str_by_range(S_target( 2: 8:2), color_limits, colors)
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 3, colored_str_by_range(S_target(10:16:2), color_limits, colors)
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 4, colored_str_by_range(S_target( 9:16:2), color_limits, colors)
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 5, colored_str_by_range(S_target(18:24:2), color_limits, colors)
      write(*,FMT='(" target  :: ",I0," :: ",A)' ) 6, colored_str_by_range(S_target(17:24:2), color_limits, colors)
      print *,cstr('---', 'blue')
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 1, colored_str_by_range(S_diff( 1: 8:2), diff_color_limits, diff_colors)
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 2, colored_str_by_range(S_diff( 2: 8:2), diff_color_limits, diff_colors)
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 3, colored_str_by_range(S_diff(10:16:2), diff_color_limits, diff_colors)
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 4, colored_str_by_range(S_diff( 9:16:2), diff_color_limits, diff_colors)
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 5, colored_str_by_range(S_diff(18:24:2), diff_color_limits, diff_colors)
      write(*,FMT='(" diff    :: ",I0," :: ",A)' ) 6, colored_str_by_range(S_diff(17:24:2), diff_color_limits, diff_colors)
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
      print*,'---------------------'
      print*,''

      @assertEqual(S_target, S, test_atol, local_msgS )
      @assertLessThanOrEqual   (zero, S)
      @assertGreaterThanOrEqual(one , S)

      @assertEqual(T_target, T, test_atol, local_msgT )
      @assertLessThanOrEqual   (zero, T)
      @assertGreaterThanOrEqual(one , T)
    endif
  end subroutine

end module
