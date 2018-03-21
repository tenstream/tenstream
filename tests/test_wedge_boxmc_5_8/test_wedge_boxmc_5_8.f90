module test_wedge_boxmc_5_8
  use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_8
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_helper_functions, only : itoa

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(8),T(5), S_target(8), T_target(5)
  real(ireals) :: S_tol(8),T_tol(5)

  type(t_boxmc_wedge_5_8) :: bmc_wedge_5_8

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-3, rtol=1e-2
  !real(ireals),parameter :: atol=1e-4, rtol=1e-2
contains

  @before
  subroutine setup(this)
      class (MpiTestMethod), intent(inout) :: this
      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      call init_mpi_data_parameters(comm)

      call bmc_wedge_5_8%init(comm)

      if(myid.eq.0) print *,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

      phi   =  0
      theta = 45

      dx = 100
      dy = dx
      dz = 50

      S_target = zero
      T_target = zero

      ! computed target with stddev_atol=5e-6, stddev_rtol=1e-4 in optprop_parameters
      ! inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) !
      !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol,inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      myid     = this%getProcessRank()
      if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1,2])
  subroutine test_wedgemc_direct_negative_azimuth_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      real(ireals) :: tau
      myid     = this%getProcessRank()

      bg  = [1e-3_ireals, zero, zero]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)

      phi = 60; theta = 90
      T_target = zero
      T_target([4]) = (sinh(tau)-cosh(tau)+1)/tau

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_60')

      phi = -60; theta = 90
      T_target = zero
      T_target([3]) = (sinh(tau)-cosh(tau)+1)/tau

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_-60')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_src4(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=4
      real(ireals) :: tau
      myid     = this%getProcessRank()

      ! direct to diffuse tests

      ! down along face 4
      phi = 240; theta = 0

      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau

      S_target = [0.0006, 0.0012, 0.0007, 0.0012, 0.0007, 0.0101, 0.0021, 0.0075]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_1')

      ! up along face 4
      theta=180
      T_target = zero
      T_target(1) = (sinh(tau)-cosh(tau)+1)/tau
      S_target = [0.0075, 0.0007, 0.0012, 0.0007, 0.0012, 0.0021, 0.0101, 0.0006]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_2')


      ! down along face 4
      phi = 240; theta = 0
      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      T_target = zero
      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau
      S_target = [0.0022, 0.0012, 0.0018, 0.0012, 0.0017, 0.0061, 0.0061, 0.0036]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_3')

      ! up along face 4
      theta=180
      T_target = zero
      T_target(1) = (sinh(tau)-cosh(tau)+1)/tau
      S_target = [0.0036, 0.0018, 0.0012, 0.0017, 0.0012, 0.0061, 0.0061, 0.0022]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_4')


      ! outwards from face 4
      phi = 240; theta = 90

      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([2,3]) = (sinh(tau)-cosh(tau)+1)/tau/2

      S_target = [0.0062, 0.0063, 0.0062, 0.0062, 0.0062, 0.0012, 0.0012, 0.0063]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_5')

      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.0076, 0.0041, 0.004, 0.0041, 0.004, 0.0042, 0.0042, 0.0077]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_6')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_src3(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=3
      real(ireals) :: tau
      myid     = this%getProcessRank()

      ! going towards the src face should not give any fluxes to anywhere
      phi = 30; theta = 0

      bg  = [1e-3_ireals, zero, zero ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero

      S_target = zero

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_4')

      ! down along face 3
      phi = 120; theta = 0

      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau

      S_target = [0.0006, 0.0012, 0.0007, 0.0101, 0.0021, 0.0012, 0.0007, 0.0074]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_2')


      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.0022, 0.0011, 0.0018, 0.0061, 0.0061, 0.0011, 0.0017, 0.0035]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_1')

      ! outwards from face 3
      phi = 120; theta = 90

      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([2,4]) = (sinh(tau)-cosh(tau)+1)/tau/2

      S_target = [0.0062, 0.0062, 0.0062, 0.0012, 0.0012, 0.0062, 0.0062, 0.0062]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_4')


      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.0076, 0.0041, 0.004, 0.0042, 0.0042, 0.0041, 0.0041, 0.0076]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_3')


  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      real(ireals) :: tau
      myid     = this%getProcessRank()

      ! direct to diffuse tests
      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]

      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau

      ! down along face 2
      phi = 0; theta = 0

      S_target = [0.0006, 0.0101, 0.0021, 0.0012, 0.0007, 0.0012, 0.0007, 0.0074]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_1')


      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.0022, 0.0061, 0.0061, 0.0011, 0.0018, 0.0011, 0.0017, 0.0035]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_2')

      ! straight outwards from face 2
      phi = 0; theta = 90

      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([3,4]) = (sinh(tau)-cosh(tau)+1)/tau/2

      S_target = [0.0062, 0.0012, 0.0012, 0.0062, 0.0062, 0.0062, 0.0062, 0.0062]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_3')


      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.0076, 0.0042, 0.0042, 0.0041, 0.004, 0.0041, 0.0041, 0.0076]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_4')

      ! outwards from face 2 towards the right face(4)
      phi = 60; theta = 90

      bg  = [1e-3_ireals, zero, one/2 ]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([4]) = (sinh(tau)-cosh(tau)+1)/tau

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_3')


      ! outwards from face 2 towards the right face(3)
      phi = -60; theta = 90

      bg  = [1e-3_ireals, zero, zero ]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([3]) = (sinh(tau)-cosh(tau)+1)/tau

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_4')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_direct_src1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1

      ! direct to diffuse tests, straight down
      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]

      phi = 0; theta = 0
      T_target = zero; T_target(5) = exp(-(bg(1)+bg(2))*dz)
      S_target = [0.0026, 0.0057, 0.0018, 0.0056, 0.0018, 0.0057, 0.0019, 0.0212]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_2')


      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
      S_target = [0.009, 0.0048, 0.0047, 0.0048, 0.0047, 0.0048, 0.0047, 0.0088]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_1')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_direct_lambert_beer(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: src, iphi
      real(ireals) :: tau
      logical, parameter :: lcheckdownward=.True., lcheckupward=.True., lchecksideward=.True.

      ! direct tests
      bg  = [1e-3, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      if(lcheckdownward) then
        !Should be rotationally symmetric for sza=0
        do iphi = 0, 360, 10
          phi = real(iphi, ireals)
          T_target = zero

          print *,'downward phi', phi
          ! Integral from top face, towards the bottom face
          do src=1,1
            T_target(5) = exp(- (bg(1)+bg(2))*dz )
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer phi '//itoa(iphi))
          enddo
        enddo

        do iphi = -89, 89, 1
          ! Integral along each of the faces, towards the bottom face
          do src=2,4
            print *,'downward along sides', src, ': phi' ,phi
            tau = bg(1) * dz
            T_target(5)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
        do iphi = 31, 209, 1
          ! Integral along each of the faces, towards the bottom face
          do src=2,2
            print *,'downward along sides', src, ': phi' ,phi
            tau = bg(1) * dz
            T_target(5)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
        do iphi = 121, 329, 1
          ! Integral along each of the faces, towards the bottom face
          do src=2,2
            print *,'downward along sides', src, ': phi' ,phi
            tau = bg(1) * dz
            T_target(5)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(5)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
      endif ! lcheckdownward

      if(lcheckupward) then
        ! and the same for upward propagation
        theta = 180
        !Should be rotationally symmetric for sza=180
        do iphi = 0, 360, 10
          phi = real(iphi, ireals)
          T_target = zero

          print *,'upward phi', phi
          ! Integral from bot face, towards the top face
          do src=5,5
            T_target(1) = exp(- (bg(1)+bg(2))*dz )
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer upward phi '//itoa(iphi))
          enddo
        enddo

        do iphi = -89, 89, 1
          ! Integral along each of the faces, towards the top face
          do src=2,4
            tau = bg(1) * dz
            T_target(1)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
        do iphi = 31, 209, 1
          ! Integral along each of the faces, towards the top face
          do src=2,2
            print *,'downward along sides', src, ': phi' ,phi
            tau = bg(1) * dz
            T_target(1)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
        do iphi = 121, 329, 1
          ! Integral along each of the faces, towards the top face
          do src=2,2
            print *,'downward along sides', src, ': phi' ,phi
            tau = bg(1) * dz
            T_target(1)=zero
            if(phi.ge.  0._ireals.and.phi.lt. 90._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.270._ireals.and.phi.le.360._ireals.and.src.eq.2) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt. 30._ireals.and.phi.lt.210._ireals.and.src.eq.3) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            if(phi.gt.120._ireals.and.phi.lt.330._ireals.and.src.eq.4) T_target(1)=(sinh(tau)-cosh(tau)+1)/tau
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides '//itoa(src)//' phi '//itoa(iphi))
          enddo
        enddo
      endif

      if(lchecksideward) then
        ! One check that if we start from a side face with 90 degree zenith, we should have equally much on the two opposite faces
        T_target = zero
        phi = 0; theta = 90
        src = 2

        tau = bg(1) * sqrt(dy**2 - (dx/2)**2)
        t_target([3,4]) = (sinh(tau)-cosh(tau)+1) / tau / 2
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')

        phi = 120
        src = 3
        T_target = zero
        T_target([2,4]) = (sinh(tau)-cosh(tau)+1) / tau / 2
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')


        phi = 240
        src = 4
        T_target = zero
        T_target([2,3]) = (sinh(tau)-cosh(tau)+1) / tau / 2
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')


        ! Or start the photons at the top and they should still go to the side faces
        T_target = zero
        T_target([3,4]) = 0.485865
        phi = 0; theta = 90
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 101')
        @assertEqual(T(3), T(4), 3*atol, 'stream should be same 101')

        T_target = zero
        T_target([2,4]) = 0.485865
        phi = 120; theta = 90
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 102')
        @assertEqual(T(2), T(4), 3*atol, 'stream should be same 102')


        T_target = zero
        T_target([2,3]) = 0.485865
        phi = 240; theta = 90
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 103')
        @assertEqual(T(2), T(3), 3*atol, 'stream should be same 103')
      endif ! lchecksideward
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0909, 0.2363, 0.018, 0.2364, 0.0179, 0.2364, 0.018, 0.1114]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.2730, 0.0, 0.2730, 0.0, 0.2730, 0.0, 0.1461]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_1')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src8(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=8
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.1114, 0.0180, 0.2363, 0.0179, 0.2364, 0.0184, 0.2364, 0.0909]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src8_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.1461, 0.0, 0.2730, 0.0, 0.2730, 0.0, 0.2730, 0.0]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src8_1')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.0, 0.0, 0.2421, 0.0, 0.2421, 0.0, 0.4826]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0347, 0.0435, 0.0484, 0.1899, 0.0212, 0.1901, 0.0212, 0.4182]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, 0.0, 0.0, 0.242, 0.0, 0.242, 0.0, 0.4828]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_3')

  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src3(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=3
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.4826, 0.0, 0.0, 0.0, 0.2421, 0.0, 0.2421, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.4182, 0.0484, 0.0435, 0.0212, 0.1899, 0.0212, 0.1901, 0.0347]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.4826, 0.0, 0.0, 0.0, 0.2421, 0.0, 0.2421, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_3')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src4(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=4
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.2421, 0.0, 0.0, 0.0, 0.2421, 0.0, 0.4826]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0347, 0.1899, 0.0212, 0.0435, 0.0484, 0.1901, 0.0212, 0.4182]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, 0.242, 0.0, 0.0, 0.0, 0.242, 0.0, 0.4828]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_3')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src5(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=5
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.4826, 0.0, 0.2421, 0.0, 0.0, 0.0, 0.2421, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.4182, 0.0212, 0.1899, 0.0484, 0.0435, 0.0212, 0.1901, 0.0347]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.4826, 0.0, 0.2421, 0.0, 0.0, 0.0, 0.2421, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_3')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src6(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=6
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.2421, 0.0, 0.2421, 0.0, 0.0, 0.0, 0.4826]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0347, 0.1899, 0.0212, 0.1901, 0.0212, 0.0435, 0.0484, 0.4182]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, 0.242, 0.0, 0.242, 0.0, 0.0, 0.0, 0.4828]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_3')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src7(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=7
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.4826, 0.0, 0.2421, 0.0, 0.2421, 0.0, 0.0, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.4182, 0.0212, 0.1899, 0.0212, 0.1901, 0.0484, 0.0435, 0.0347]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.4826, 0.0, 0.2421, 0.0, 0.2421, 0.0, 0.0, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_3')
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
