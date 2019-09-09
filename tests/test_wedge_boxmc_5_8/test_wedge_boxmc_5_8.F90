module test_wedge_boxmc_5_8
  use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_8
  use m_data_parameters, only :          &
    mpiint, iintegers,                   &
    ireals, ireal_dp, irealLUT, ireal_params, &
    one, zero, i1, default_str_len,      &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_helper_functions, only : itoa, ftoa, triangle_area_by_vertices, rad2deg
  use m_boxmc_geometry, only : setup_default_unit_wedge_geometry, setup_default_wedge_geometry
  use m_optprop_LUT, only : azimuth_from_param_phi
  use m_LUT_param_phi, only: LUT_wedge_dz

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(8),T(5), S_target(8), T_target(5)
  real(ireals) :: S_tol(8),T_tol(5)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_wedge_5_8) :: bmc_wedge_5_8

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-1, theta_eps=1e-3
  !real(ireal_dp),parameter :: atol=1e-4, rtol=1e-2
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

      call setup_default_unit_wedge_geometry(dx, dy, dz, vertices)

      S_target = zero
      T_target = zero

      ! computed target with stddev_atol=5e-6, stddev_rtol=1e-4 in optprop_parameters
      ! inp_atol=1e-6_ireal_dp, inp_rtol=1e-4_ireal_dp) !
      !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,vertices,S,T,S_tol,T_tol,inp_atol=1e-6_ireal_dp, inp_rtol=1e-4_ireal_dp) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(mpiint) :: ierr
      myid     = this%getProcessRank()
      call PetscFinalize(ierr)
      if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1])
  subroutine test_wedgemc_no_absorption(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: isrc, itau
      real(ireal_dp) :: tau, w0, g, aspectzx, phi, theta, Cx, Cy, verts(3*6)

      myid     = this%getProcessRank()

      w0       = 1
      g        = 0
      aspectzx = 1
      !Cx       = 0.500000000
      !Cy       = 0.866
      Cx       = 1.00000000
      Cy       = 0.5
      phi      = 270
      theta    = 20
      verts = [&
         0._ireal_dp , 0._ireal_dp , 0._ireal_dp       , &
         1._ireal_dp  , 0._ireal_dp , 0._ireal_dp       , &
         Cx   , Cy   , 0._ireal_dp       , &
         0._ireal_dp , 0._ireal_dp , aspectzx   , &
         1._ireal_dp  , 0._ireal_dp , aspectzx   , &
         Cx   , Cy   , aspectzx ]


      do itau=-10,2,4
        tau = 10._ireal_dp**(itau)
        do isrc=1,5
          bg  = [tau*(1-w0), tau*w0, g]
          call bmc_wedge_5_8%get_coeff(comm,bg,isrc,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          print *,'tau', tau, 'src', isrc, 'T', sum(T), 'S', sum(S), 'divergence', 1 - (sum(T)+sum(S))
          if(any(T_tol.gt.0) .or. any(S_tol.gt.0)) then
            @assertLessThan(1 - (sum(T)+sum(S)), 1e-6_ireals, 'given w0 is 1, there should be no divergence' )
          endif
        enddo
      enddo
  end subroutine

  @test(npes =[1])
  subroutine test_wedgemc_direct_custom_ex1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: isrc=2
      real(ireal_dp) :: tau, w0, g, aspectzx, phi, theta, Cx, Cy, verts(3*6)
      myid     = this%getProcessRank()

      tau      = 9.60528996E-06
      w0       = 0.00000000
      g        = 0.00000000
      aspectzx = 0.217637643
      Cx       = 0.500000000
      Cy       = 0.948499978
      phi      = rad2deg(azimuth_from_param_phi(.75_ireal_params, real([Cx, Cy], ireal_params)))
      theta    = 56.3100014
      !verts = [&
      !   0.00000000  ,     0.00000000  ,     0.00000000,&
      !   1.00000000  ,     0.00000000  ,     0.00000000,&
      !   0.514306128 ,     0.947683156 ,     0.00000000,&
      !   0.00000000  ,     0.00000000  ,     1.08398283,&
      !   1.00000000  ,     0.00000000  ,     1.08398283,&
      !   0.514306128 ,     0.947683156 ,     1.08398283 ]
      verts = [&
         0._ireal_dp , 0._ireal_dp , 0._ireal_dp       , &
         1._ireal_dp  , 0._ireal_dp , 0._ireal_dp       , &
         Cx   , Cy   , 0._ireal_dp       , &
         0._ireal_dp , 0._ireal_dp , aspectzx   , &
         1._ireal_dp  , 0._ireal_dp , aspectzx   , &
         Cx   , Cy   , aspectzx ]

      bg  = [tau*(1._ireal_dp-w0), tau*w0, g]

      S_target = zero
      T_target = zero
      T_target([4]) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      do isrc=1,5
        call bmc_wedge_5_8%get_coeff(comm,bg,isrc,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        print *,'src', isrc, 'T', T
      enddo
  end subroutine

  @test(npes =[1])
  subroutine test_wedgemc_direct_negative_azimuth_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      real(ireal_dp) :: tau
      myid     = this%getProcessRank()

      bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)

      phi = 60; theta = 90
      T_target = zero
      T_target([4]) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_60')

      phi = -60; theta = 90
      T_target = zero
      T_target([3]) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_-60')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_src4(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=4
      real(ireal_dp) :: tau
      myid     = this%getProcessRank()

      ! direct to diffuse tests

      ! down along face 4
      phi = 240; theta = theta_eps

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      S_target = [0.0006, 0.0012, 0.0007, 0.0012, 0.0007, 0.0101, 0.0021, 0.0075]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_1')

      ! up along face 4
      theta = 180 - theta_eps
      T_target = zero
      T_target(1) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)
      S_target = [0.0075, 0.0007, 0.0012, 0.0007, 0.0012, 0.0021, 0.0101, 0.0006]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_2')


      ! down along face 4
      phi = 240; theta = theta_eps
      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      T_target = zero
      T_target(5) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)
      S_target = [0.0022, 0.0012, 0.0018, 0.0012, 0.0017, 0.0061, 0.0061, 0.0036]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_3')

      ! up along face 4
      theta = 180 - theta_eps
      T_target = zero
      T_target(1) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)
      S_target = [0.0036, 0.0018, 0.0012, 0.0017, 0.0012, 0.0061, 0.0061, 0.0022]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_4')


      ! outwards from face 4
      phi = 240; theta = 90 - theta_eps

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([2,3]) = real((sinh(tau)-cosh(tau)+1)/tau/2, ireals)

      S_target = [0.0062, 0.0063, 0.0062, 0.0062, 0.0062, 0.0012, 0.0012, 0.0063]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_5')

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.0076, 0.0041, 0.004, 0.0041, 0.004, 0.0042, 0.0042, 0.0077]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_6')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_src3(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=3
      real(ireal_dp) :: tau
      myid     = this%getProcessRank()

      ! going towards the src face should not give any fluxes to anywhere
      phi = -30; theta = 0

      bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero

      S_target = zero

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_4')

      ! down along face 3
      phi = 120; theta = theta_eps

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]
      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      S_target = [0.0006, 0.0012, 0.0007, 0.0101, 0.0021, 0.0012, 0.0007, 0.0074]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_2')


      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.0022, 0.0011, 0.0018, 0.0061, 0.0061, 0.0011, 0.0017, 0.0035]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_1')

      ! outwards from face 3
      phi = 120; theta = 90

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([2,4]) = real((sinh(tau)-cosh(tau)+1)/tau/2, ireals)

      S_target = [0.0062, 0.0062, 0.0062, 0.0012, 0.0012, 0.0062, 0.0062, 0.0062]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_4')


      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.0076, 0.0041, 0.004, 0.0042, 0.0042, 0.0041, 0.0041, 0.0076]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_3')


  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      real(ireal_dp) :: tau
      myid     = this%getProcessRank()

      ! direct to diffuse tests
      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]

      tau = (bg(1)+bg(2)) * dz
      T_target = zero
      T_target(5) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      ! down along face 2
      phi = 0; theta = theta_eps

      S_target = [0.0006, 0.0101, 0.0021, 0.0012, 0.0007, 0.0012, 0.0007, 0.0074]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_1')


      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.0022, 0.0061, 0.0061, 0.0011, 0.0018, 0.0011, 0.0017, 0.0035]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_2')

      ! straight outwards from face 2
      phi = 0; theta = 90

      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]
      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([3,4]) = real((sinh(tau)-cosh(tau)+1)/tau/2, ireals)

      S_target = [0.0062, 0.0012, 0.0012, 0.0062, 0.0062, 0.0062, 0.0062, 0.0062]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_3')


      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.0076, 0.0042, 0.0042, 0.0041, 0.004, 0.0041, 0.0041, 0.0076]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_4')

      ! outwards from face 2 towards the right face(4)
      phi = 60; theta = 90 - theta_eps

      bg  = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp/2 ]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([4]) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_3')


      ! outwards from face 2 towards the right face(3)
      phi = -60; theta = 90

      bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp ]
      S_target = zero

      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
      T_target = zero
      T_target([3]) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_4')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_src1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1

      ! direct to diffuse tests, straight down
      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 1._ireal_dp/2 ]

      phi = 0; theta = 0
      T_target = zero; T_target(5) = real(exp(-(bg(1)+bg(2))*dz), ireals)
      S_target = [0.0026, 0.0057, 0.0018, 0.0056, 0.0018, 0.0057, 0.0019, 0.0212]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_2')


      bg  = [1e-3_ireal_dp, 1e-3_ireal_dp, 0._ireal_dp ]
      S_target = [0.009, 0.0048, 0.0047, 0.0048, 0.0047, 0.0048, 0.0047, 0.0088]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_1')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_direct_lambert_beer(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: src, iphi
      real(ireal_dp) :: tau
      logical, parameter :: lcheckdownward=.True., lcheckupward=.True., lchecksideward=.True.
      !logical, parameter :: lcheckdownward=.True., lcheckupward=.False., lchecksideward=.False.
      !logical, parameter :: lcheckdownward=.False., lcheckupward=.True., lchecksideward=.False.
      !logical, parameter :: lcheckdownward=.False., lcheckupward=.False., lchecksideward=.True.

      ! direct tests
      bg  = [1e-3, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      if(lcheckdownward) then
        !Should be rotationally symmetric for sza=0
        do iphi = 0, 360, 10
          phi = real(iphi, ireal_dp)
          T_target = zero

          ! Integral from top face, towards the bottom face
          do src=1,1
            T_target(5) = real(exp(- (bg(1)+bg(2))*dz ), ireals)
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer phi '//itoa(iphi))
          enddo
        enddo

        theta = theta_eps
        tau = bg(1) * dz
        do iphi = 0, 360, 10
          phi = real(iphi, ireal_dp)
          ! Integral along each of the faces, towards the bottom face
          do src=2,4
            T_target(5)=zero
            if(phi.gt.-90._ireal_dp.and.phi.lt. 90._ireal_dp.and.src.eq.2) T_target(5)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt.270._ireal_dp.and.phi.le.430._ireal_dp.and.src.eq.2) T_target(5)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt. 30._ireal_dp.and.phi.lt.210._ireal_dp.and.src.eq.3) T_target(5)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt.150._ireal_dp.and.phi.lt.330._ireal_dp.and.src.eq.4) T_target(5)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides / downward '//itoa(src)//' phi '//ftoa(phi))
          enddo
        enddo
      endif ! lcheckdownward

      if(lcheckupward) then
        ! and the same for upward propagation
        theta = 180
        !Should be rotationally symmetric for sza=180
        do iphi = 0, 360, 10
          phi = real(iphi, ireal_dp)
          T_target = zero

          ! Integral from bot face, towards the top face
          do src=5,5
            T_target(1) = real(exp(- (bg(1)+bg(2))*dz ), ireals)
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer upward phi '//itoa(iphi))
          enddo
        enddo

        theta = 180 - theta_eps
        tau = bg(1) * dz
        do iphi = 0, 360, 10
          phi = real(iphi, ireal_dp)
          ! Integral along each of the faces, towards the top face
          do src=2,4
            T_target(1)=zero
            if(phi.gt.-90._ireal_dp.and.phi.lt. 90._ireal_dp.and.src.eq.2) T_target(1)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt.270._ireal_dp.and.phi.lt.430._ireal_dp.and.src.eq.2) T_target(1)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt. 30._ireal_dp.and.phi.lt.210._ireal_dp.and.src.eq.3) T_target(1)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            if(phi.gt.150._ireal_dp.and.phi.lt.330._ireal_dp.and.src.eq.4) T_target(1)=real((sinh(tau)-cosh(tau)+1)/tau, ireals)
            call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
            call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_along_sides / upward '//itoa(src)//' phi '//ftoa(phi))
          enddo
        enddo
      endif

      if(lchecksideward) then
        ! One check that if we start from a side face with 90 degree zenith, we should have equally much on the two opposite faces
        T_target = zero
        phi = 0; theta = 90
        src = 2

        tau = bg(1) * sqrt(dy**2 - (dx/2)**2)
        t_target([3,4]) = real((sinh(tau)-cosh(tau)+1) / tau / 2, ireals)
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards src='//itoa(src)//' phi='//ftoa(phi))

        phi = 120
        src = 3
        T_target = zero
        T_target([2,4]) = real((sinh(tau)-cosh(tau)+1) / tau / 2, ireals)
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards src='//itoa(src)//' phi='//ftoa(phi))


        phi = 240
        src = 4
        T_target = zero
        T_target([2,3]) = real((sinh(tau)-cosh(tau)+1) / tau / 2, ireals)
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards src='//itoa(src)//' phi='//ftoa(phi))


        ! Or start the photons at the top and they should still go to the side faces
        T_target = zero
        T_target([3,4]) = 0.485865_ireals
        phi = 0; theta = 90 - theta_eps
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 101')
        @assertEqual(T(3), T(4), real(2*3*atol, ireals), 'stream should be same 101')

        T_target = zero
        T_target([2,4]) = 0.485865_ireals
        phi = 120; theta = 90 - theta_eps
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 102')
        @assertEqual(T(2), T(4), real(2*3*atol, ireals), 'stream should be same 102')


        T_target = zero
        T_target([2,3]) = 0.485865_ireals
        phi = 240; theta = 90 - theta_eps
        src = 1
        call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 103')
        @assertEqual(T(2), T(3), real(2*3*atol, ireals), 'stream should be same 103')
      endif ! lchecksideward
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0833, 0.2107, 0.0226, 0.2107, 0.0226, 0.2107, 0.0226, 0.1791]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.2374, 0.0, 0.2374, 0.0, 0.2374, 0.0, 0.2511]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_1')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src8(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=8
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.1792, 0.0226, 0.2107, 0.0226, 0.2107, 0.0226, 0.2107, 0.0833]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src8_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.2516, 0.0, 0.2370, 0.0, 0.2370, 0.0, 0.2370, 0.0]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src8_1')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.0, 0.0, c1, 0.0, c1, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c8, c7, c6, c5, c4, c5, c4, c3]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, 0.0, 0.0, c1, 0.0, c1, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_3')

  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src3(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=3
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [c2, 0.0, 0.0, 0.0, c1, 0.0, c1, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c3, c6, c7, c4, c5, c4, c5, c8]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [c2, 0.0, 0.0, 0.0, c1, 0.0, c1, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_3')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src4(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=4
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, c1, 0.0, 0.0, 0.0, c1, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c8, c5, c4, c7, c6, c5, c4, c3]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, c1, 0.0, 0.0, 0.0, c1, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_3')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src5(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=5
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [c2, 0.0, c1, 0.0, 0.0, 0.0, c1, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c3, c4, c5, c6, c7, c4, c5, c8]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [c2, 0.0, c1, 0.0, 0.0, 0.0, c1, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src5_3')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src6(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=6
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, c1, 0.0, c1, 0.0, 0.0, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c8, c5, c4, c5, c4, c7, c6, c3]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, c1, 0.0, c1, 0.0, 0.0, 0.0, c2]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src6_3')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src7(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=7
      real,parameter :: c1=0.2772, c2=0.4107, c3=0.3651, c4=0.0237, c5=0.2155, c6=0.0437, c7=0.0389, c8=0.0390
      myid     = this%getProcessRank()

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [c2, 0.0, c1, 0.0, c1, 0.0, 0.0, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_1')

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [c3, c4, c5, c4, c5, c6, c7, c8]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_2')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [c2, 0.0, c1, 0.0, c1, 0.0, 0.0, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src7_3')
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src1_rectilinear_C(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1
      real(ireal_dp) :: aspect_zx, dz, Atop, dtau, w0
      real(ireal_dp), dimension(2) :: A, B, C
      myid     = this%getProcessRank()

      aspect_zx = 0.473405033
      dtau = 7.63205662E-06
      w0 = 0.109288104

      A = [0._ireal_dp, 0._ireal_dp]
      B = [1._ireal_dp, 0._ireal_dp]
      C = [.0_ireal_dp, 1.0_ireal_dp]

      Atop = triangle_area_by_vertices(A,B,C)
      dz = LUT_wedge_dz(real(Atop, irealLUT), real(aspect_zx, irealLUT))
      if(allocated(vertices)) deallocate(vertices)
      call setup_default_wedge_geometry(&
        A, B, C, &
        real(dz, ireal_dp), &
        vertices)

      ! ----------------------------------

      bg = [dtau*(1._ireal_dp-w0), dtau*w0, 0._ireal_dp ]
      S_target = [0.0000, 0.2114, 0.0000, 0.2114, 0.0000, 0.2828, 0.0000, 0.2945]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_C01')

      deallocate(vertices)


      C = [1.0_ireal_dp, 1.0_ireal_dp]

      Atop = triangle_area_by_vertices(A,B,C)
      dz = LUT_wedge_dz(real(Atop, irealLUT), real(aspect_zx, irealLUT))
      if(allocated(vertices)) deallocate(vertices)
      call setup_default_wedge_geometry(&
        A, B, C, &
        real(dz, ireal_dp), &
        vertices)

      ! ----------------------------------

      bg = [dtau*(1._ireal_dp-w0), dtau*w0, 0._ireal_dp ]
      S_target = [0.0000, 0.2114, 0.0000, 0.2828, 0.0000, 0.2114, 0.0000, 0.2945]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src1_C01')

      deallocate(vertices)
  end subroutine

  subroutine check(S_target,T_target, S,T, msg)
      real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

      real(ireal_dp),parameter :: sigma = 6 ! normal test range for coefficients

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

        @assertEqual(S_target, S, real(atol*sigma, ireals), local_msgS )
        @assertLessThanOrEqual   (zero, S)
        @assertGreaterThanOrEqual(one , S)

        @assertEqual(T_target, T, real(atol*sigma, ireals), local_msgT )
        @assertLessThanOrEqual   (zero, T)
        @assertGreaterThanOrEqual(one , T)
      endif
  end subroutine

end module
