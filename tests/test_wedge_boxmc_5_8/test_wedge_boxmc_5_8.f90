module test_wedge_boxmc_5_8
  use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_8
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(8),T(5), S_target(8), T_target(5)
  real(ireals) :: S_tol(8),T_tol(5)

  type(t_boxmc_wedge_5_8) :: bmc_wedge_5_8

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-3, rtol=1e-2
  !real(ireals),parameter :: atol=1e-4, rtol=1e-3
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
      if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

!  @test(npes =[1,2])
!  subroutine test_wedgemc_direct_negative_azimuth_src2(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=2
!      real(ireals) :: tau
!
!
!      bg  = [1e-3_ireals, zero, zero]
!      S_target = zero
!
!      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
!
!      phi = 60; theta = 90
!      T_target = zero
!      T_target([4]) = (sinh(tau)-cosh(tau)+1)/tau
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_60')
!
!      phi = -60; theta = 90
!      T_target = zero
!      T_target([3]) = (sinh(tau)-cosh(tau)+1)/tau
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_negative_azimuth_src2_-60')
!  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_select_cases_direct_src4(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=4
!      real(ireals) :: tau
!
!      ! direct to diffuse tests
!
!      ! down along face 4
!      phi = 240; theta = 0
!
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!      tau = (bg(1)+bg(2)) * dz
!      T_target = zero
!      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau
!
!      S_target = [0.00017051, 0.00039874 , 0.00066926 , 0.0208261, 0.00203671]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_2')
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [0.00060366,  0.00064998 ,  0.00100554,  0.0208439,  0.00097272]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_1')
!
!      ! outwards from face 4
!      phi = 240; theta = 90
!
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
!      T_target = zero
!      T_target([2,3]) = (sinh(tau)-cosh(tau)+1)/tau/2
!
!      S_target = [ 0.00624004,  0.0124557 ,  0.0124458 ,  0.00244814,  0.00624669]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_4')
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [ 0.0076486 ,  0.00809211,  0.00811609,  0.00836656,  0.00764847 ]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src4_3')
!  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_select_cases_direct_src3(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=3
!      real(ireals) :: tau
!
!      ! down along face 3
!      phi = 120; theta = 0
!
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!      tau = (bg(1)+bg(2)) * dz
!      T_target = zero
!      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau
!
!      S_target = [0.00017051,  0.00039874,  0.0208261 ,  0.00066926,  0.00203671]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_2')
!
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [0.00060366,  0.00064998,  0.0208439 ,  0.00100554,  0.00097272]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_1')
!
!      ! outwards from face 3
!      phi = 120; theta = 90
!
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
!      T_target = zero
!      T_target([2,4]) = (sinh(tau)-cosh(tau)+1)/tau/2
!
!      S_target = [ 0.00624004,  0.0124557 ,  0.00244814,  0.0124458 ,  0.00624669]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_4')
!
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [ 0.0076486 ,  0.00809211,  0.00836656,  0.00811609,  0.00764847 ]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src3_3')
!
!
!  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_select_cases_direct_src2(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=2
!      real(ireals) :: tau
!
!      ! direct to diffuse tests
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!
!
!      ! outwards from face 2
!      phi = 0; theta = 90
!
!      tau = (bg(1)+bg(2)) * sqrt(dy**2 - (dx/2)**2)
!      T_target = zero
!      T_target([3,4]) = (sinh(tau)-cosh(tau)+1)/tau/2
!
!      S_target = [0.00623344,  0.00244899,  0.01244585,  0.0124478 ,  0.00623946]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_4')
!
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [0.00764423,  0.00835319,  0.00810088,  0.00809931,  0.00764853]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_3')
!
!
!      ! down along face 2
!      phi = 0; theta = 0
!
!      tau = (bg(1)+bg(2)) * dz
!      T_target = zero
!      T_target(5) = (sinh(tau)-cosh(tau)+1)/tau
!
!      S_target = zero; S_target(2) = 0.0241891
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_2')
!
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = zero; S_target(2) = 0.0241891
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src2_1')
!  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_select_cases_direct_src1(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=1
!
!      ! direct to diffuse tests, straight down
!      bg  = [1e-3_ireals, 1e-3_ireals, one/2 ]
!
!      phi = 0; theta = 0
!      T_target = zero; T_target(5) = exp(-(bg(1)+bg(2))*dz)
!      S_target = [ 0.00262582, 0.0074815, 0.0074815, 0.0074815, 0.0213178]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_2')
!
!
!      bg  = [1e-3_ireals, 1e-3_ireals, zero ]
!      S_target = [0.0090376, 0.009520906, 0.009520906, 0.009520906, 0.00875678]
!
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_src1_1')
!  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_direct_lambert_beer(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers) :: src, iphi
!      real(ireals) :: tau
!
!      ! direct tests
!      bg  = [1e-3, 0., 0. ]
!      phi = 0; theta = 0
!      S_target = zero
!
!      !Should be rotationally symmetric for sza=0
!      do iphi = 0, 360, 10
!        phi = iphi
!        T_target = zero
!
!        print *,'downward'
!        ! Integral from top face, towards the bottom face
!        do src=1,1
!          T_target(5) = exp(- (bg(1)+bg(2))*dz )
!          call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!          call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer')
!        enddo
!
!        ! Integral along each of the faces, towards the bottom face
!        do src=2,4
!          print *,'downward along sides', src
!          tau = bg(1) * dz
!          T_target(5) = (sinh(tau)-cosh(tau)+1)/tau
!          call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!          call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer')
!        enddo
!      enddo
!
!      print *,'upward'
!      ! and the same for upward propagation
!      theta = 180
!      do iphi = 0, 360, 10
!        phi = iphi
!        T_target = zero
!
!        do src=5,5
!          T_target(1) = exp(- (bg(1)+bg(2))*dz )
!          call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!          call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_upward')
!        enddo
!
!        print *,'upward Integral along each of the faces, towards the bottom face'
!        ! Integral along each of the faces, towards the bottom face
!        do src=2,4
!          tau = bg(1) * dz
!          T_target(1) = (sinh(tau)-cosh(tau)+1)/tau
!          call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!          call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_upward')
!        enddo
!      enddo
!
!
!      ! One check that if we start from a side face with 90 degree zenith, we should have equally much on the two opposite faces
!      T_target = zero
!      phi = 0; theta = 90
!      src = 2
!
!      tau = bg(1) * sqrt(dy**2 - (dx/2)**2)
!      t_target([3,4]) = (sinh(tau)-cosh(tau)+1) / tau / 2
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')
!
!      phi = 120
!      src = 3
!      T_target = zero
!      T_target([2,4]) = (sinh(tau)-cosh(tau)+1) / tau / 2
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')
!
!
!      phi = 240
!      src = 4
!      T_target = zero
!      T_target([2,3]) = (sinh(tau)-cosh(tau)+1) / tau / 2
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_sidewards')
!
!
!      ! Or start the photons at the top and they should still go to the side faces
!      T_target = zero
!      T_target([3,4]) = (4.85805E-01+4.85883E-01)/2
!      phi = 0; theta = 90
!      src = 1
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 101')
!      @assertEqual(T(3), T(4), 3*atol, 'stream should be same 101')
!
!      T_target = zero
!      T_target([2,4]) = (4.85805E-01+4.85883E-01)/2
!      phi = 120; theta = 90
!      src = 1
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 102')
!      @assertEqual(T(2), T(4), 3*atol, 'stream should be same 102')
!
!
!      T_target = zero
!      T_target([2,3]) = (4.85805E-01+4.85883E-01)/2
!      phi = 240; theta = 90
!      src = 1
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_direct_lambert_beer_top_plate_towards sidefaces 103')
!      @assertEqual(T(2), T(3), 3*atol, 'stream should be same 103')
!  end subroutine
!
!  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src1(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=1

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

  !@test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src2(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=2

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.0347, 0.0435, 0.0484, 0.1899, 0.0212, 0.1901, 0.0212, 0.4182]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_3')

      ! ----------------------------------
      bg = [1e-3, 1e-2, 1.0 ] ! in case of pure forward scattering, it should be the same as just absorption
      S_target = [0.0, 0.0, 0.0, 0.242, 0.0, 0.242, 0.0, 0.4828]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.0, 0.0, 0.0, 0.2421, 0.0, 0.2421, 0.0, 0.4826]

      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src2_1')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_select_cases_diffuse_src3(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers),parameter :: src=3

      ! ----------------------------------
      bg = [1e-3, 1e-2, .0 ]
      S_target = [0.4182, 0.0484, 0.0435, 0.0212, 0.1899, 0.0212, 0.1901, 0.4182]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_2')

      ! ----------------------------------
      bg = [1e-3, 0., .0 ]
      S_target = [0.4826, 0.0, 0.0, 0.0, 0.2421, 0.0, 0.2421, 0.0]
      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src3_1')
  end subroutine
!
!  @test(npes =[1,2])
!  subroutine test_boxmc_select_cases_diffuse_src4(this)
!      class (MpiTestMethod), intent(inout) :: this
!      integer(iintegers),parameter :: src=4
!
!      ! ----------------------------------
!      bg = [1e-3, 1e-2, .0 ]
!      S_target = [ 0.226478, 0.21119225, 0.21119225, 0.0919335, 0.226539 ]
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_2')
!
!      ! ----------------------------------
!      bg = [1e-3, 0., .0 ]
!      S_target = [0.241484, 0.242103, 0.242103, 0., 0.241141]
!      call bmc_wedge_5_8%get_coeff(comm,bg,src,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
!      call check(S_target,T_target, S,T, msg='test_wedgemc_diffuse_src4_1')
!  end subroutine

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
