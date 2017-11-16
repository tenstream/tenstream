module test_boxmc_8_10
  use m_boxmc, only : t_boxmc, t_boxmc_8_10
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(8), S_target(10), T_target(8)
  real(ireals) :: S_tol(10),T_tol(6)

  type(t_boxmc_8_10) :: bmc_8_10

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-3, rtol=1e-2
  ! real(ireals),parameter :: atol=1e-4, rtol=1e-3
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_8_10%init(comm)

    if(myid.eq.0) print *,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 50

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown


  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireals, zero, one/2 ]

    ! from top to bot face
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    S_target = zero


    do src = 1,4
      T_target = zero
      T_target(src) = exp(-tau)
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_direct_srctopface top_to_bot')
    enddo
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_srctopface_45(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireals, zero, one/2 ]

    ! outwards from face 2
    phi = 0; theta = 45*one

    tau = (bg(1)+bg(2)) * dz * sqrt(2*one)

    S_target = zero

    T_target = zero
    T_target(3) = exp(-tau)

    src = 1
    call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_direct_srctopface_45')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_srcsidefaces(this)
    class (MpiTestMethod), intent(inout)  :: this
    integer(iintegers)                    :: src, iphi
    real(ireals)                          :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireals, zero, one/2 ]


    ! along each of the side faces
    phi = 0; theta = 90

    tau = (bg(1)+bg(2)) * dx

    S_target = zero
    T_target = zero

    do src = 5,6,1
      phi = 0

      T_target = zero

      T_target(src+2) =  (sinh(tau)-cosh(tau)+1)/tau
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces')

      phi = 90

      T_target = zero

      T_target(src) =  exp(-tau)
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces')
    enddo

    do src = 7,8
      phi = 0

      T_target = zero

      T_target(src) = exp(-tau)
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces')

      phi = 90

      T_target = zero

      T_target(src-2) = (sinh(tau)-cosh(tau)+1)/tau
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces')
    enddo

    do iphi = 0,360,60
      phi = iphi
      theta= 0

      tau = (bg(1)+bg(2)) * dz/2

      src = 5
      T_target = zero
      T_target([1,3]) = (sinh(tau)-cosh(tau)+1)/tau/2
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 6
      T_target = zero
      T_target([1,3]) = (sinh(tau)-cosh(tau)+1)/tau/2 * exp(-tau)
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 7
      T_target = zero
      T_target([1,2]) = (sinh(tau)-cosh(tau)+1)/tau/2
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')

      src = 8
      T_target = zero
      T_target([1,2]) = (sinh(tau)-cosh(tau)+1)/tau/2 * exp(-tau)
      call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces theta = 0 down along sidefaces')
    enddo
  end subroutine


  @test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireals, zero, zero ]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    S_target = [0.0, 0.24254, 0.17735, 0.17735, 0.0, 0.0, 0.17735, 0.17735, 0.0, 0.0]

    T_target = zero


    src = 2   !top face
    call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diff_srcbottomface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireals, zero, zero ]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    T_target = zero

    S_target = [0.24254, 0.0, 0.0, 0.0, 0.17735, 0.17735, 0.0, 0.0, 0.17735, 0.17735]

    src = 1
    call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbottomface')
  end subroutine

   @test(npes =[1,2])
   subroutine test_boxmc_select_cases_diff_srcsideface(this)
     class (MpiTestMethod), intent(inout) :: this
     integer(iintegers) :: src
     real(ireals) :: tau

     ! direct to diffuse tests
     bg  = [1e-3_ireals, zero, zero ]

     ! outwards from face 2
     phi = 0; theta = 0

     tau = (bg(1)+bg(2)) * dz
     T_target = zero

     src = 3
     S_target = [0.0, 0.592867, 0.05618, 0.0, 0.0, 0.0, 0.15440, 0.15440, 0.0, 0.0]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 3')

     src = 4
     S_target = [0.0, 0.592867, 0.0, 0.05618, 0.0, 0.0, 0.15440, 0.15440, 0.0, 0.0]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 4')

     src = 5
     S_target = [0.592867, 0.0, 0.0, 0.0, 0.05618, 0.0, 0.0, 0.0, 0.15440, 0.15440]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 5')

     src = 6
     S_target = [0.592867, 0.0, 0.0, 0.0, 0.0, 0.05618, 0.0, 0.0, 0.15440, 0.15440]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 6')

     src = 7
     S_target = [0.0, 0.592867, 0.15440, 0.15440, 0.0, 0.0, 0.05618, 0.0, 0.0, 0.0]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 7')

     src = 8
     S_target = [0.0, 0.592867, 0.15440, 0.15440, 0.0, 0.0, 0.0, 0.05618, 0.0, 0.0]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 8')

     src = 9
     S_target = [0.592867, 0.0, 0.0, 0.0, 0.15440, 0.15440, 0.0, 0.0, 0.05618, 0.0]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 9')

     src = 10
     S_target = [0.592867, 0.0, 0.0, 0.0, 0.15440, 0.15440, 0.0, 0.0, 0.0, 0.05618]
     call bmc_8_10%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 10')

   end subroutine


  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

    real(ireals),parameter :: sigma = 3 ! normal test range for coefficients

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
      write(*, FMT='( " diffuse ::  :: ",10(es12.5) )' ) S
      write(*, FMT='( " target  ::  :: ",10(es12.5) )' ) S_target
      write(*, FMT='( " diff    ::  :: ",10(es12.5) )' ) S_target-S
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
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
