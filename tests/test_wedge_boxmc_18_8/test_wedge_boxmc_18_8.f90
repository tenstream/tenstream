module test_wedge_boxmc_18_8
  use m_boxmc, only : t_boxmc,t_boxmc_wedge_18_8
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_helper_functions, only : itoa, triangle_area_by_vertices
  use m_boxmc_geometry, only : setup_default_unit_wedge_geometry, setup_default_wedge_geometry

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(8),T(18), S_target(8), T_target(18)
  real(ireals) :: S_tol(8),T_tol(18)
  real(ireals), allocatable :: vertices(:)

  type(t_boxmc_wedge_18_8) :: bmc_wedge_18_8

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

      call bmc_wedge_18_8%init(comm)

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
      ! inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) !
      !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,vertices,S,T,S_tol,T_tol,inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(mpiint) :: ierr
      myid     = this%getProcessRank()
      call PetscFinalize(ierr)
      if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1,2])
  subroutine test_boxmc_direct_lambert_beer_top_to_bottom(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: src

      ! direct tests
      bg  = [1e-3, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      do src = 1, 3
        T_target = zero
        T_target(15+src) = exp(-(bg(1)+bg(2))*dz)
        call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
          phi,theta,vertices,S,T,S_tol,T_tol, &
          inp_atol=atol, inp_rtol=rtol)

        call check(S_target,T_target, S,T, &
          msg='test_wedgemc_direct_lambert_beer_top_to_bottom '//itoa(src))
      enddo


      theta = 180
      do src = 16, 18
        T_target = zero
        T_target(-15+src) = exp(-(bg(1)+bg(2))*dz)
        call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
          phi,theta,vertices,S,T,S_tol,T_tol, &
          inp_atol=atol, inp_rtol=rtol)

        call check(S_target,T_target, S,T, &
          msg='test_wedgemc_direct_lambert_beer_bottom_to_top '//itoa(src))
      enddo
  end subroutine


  @test(npes =[1,2])
  subroutine test_boxmc_direct_lambert_beer_side_to_bottom(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: src
      real(ireals) :: tau

      ! direct tests
      bg  = [1e-3, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      tau = (bg(1)+bg(2))*dz / 2

      src = 4
      T_target = zero; T_target(16) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 5
      T_target = zero; T_target(17) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))



      src = 6
      T_target = zero; T_target(16) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 7
      T_target = zero; T_target(17) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))


      phi = 120
      src = 8
      T_target = zero; T_target(16) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 9
      T_target = zero; T_target(18) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))



      src = 10
      T_target = zero; T_target(16) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 11
      T_target = zero; T_target(18) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      phi = 240
      src = 12
      T_target = zero; T_target(18) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 13
      T_target = zero; T_target(17) = (sinh(tau)-cosh(tau)+1) / tau * exp(-tau)
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))



      src = 14
      T_target = zero; T_target(18) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))

      src = 15
      T_target = zero; T_target(17) = (sinh(tau)-cosh(tau)+1) / tau
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_side_to_bottom '//itoa(src))
  end subroutine


  @test(npes =[1,2])
  subroutine test_boxmc_direct_lambert_beer_top_to_side(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: src

      ! direct tests
      bg  = [1e-3, 0., 0. ]
      phi = 0; theta = 60
      S_target = zero

      src = 1
      T_target = [0.00000, 0.00000, 0.00000, &
                  0.00000, 0.00000, 0.00000, 0.00000, &
                  0.36898, 0.27204, 0.00000, 0.32114, &
                  0.00000, 0.00000, 0.00000, 0.00000, &
                  0.00000, 0.00000, 0.00000]
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_top_to_side '//itoa(src))



      src = 2
      T_target = [0.00000, 0.00000, 0.00000, &
                  0.00000, 0.00000, 0.00000, 0.00000, &
                  0.00000, 0.00000, 0.00000, 0.00000, &
                  0.27207, 0.36867, 0.32140, 0.00000, &
                  0.00000, 0.00000, 0.00000]
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_top_to_side '//itoa(src))


      src = 3
      T_target = [0.00000, 0.00000, 0.00000, &
                  0.00000, 0.00000, 0.00000, 0.00000, &
                  0.00000, 0.45943, 0.00000, 0.02952, &
                  0.45960, 0.00000, 0.02958, 0.00000, &
                  0.00000, 0.00000, 0.00000]
      call bmc_wedge_18_8%get_coeff(comm,bg,src,.True.,&
        phi,theta,vertices,S,T,S_tol,T_tol, &
        inp_atol=atol, inp_rtol=rtol)

      call check(S_target,T_target, S,T, &
        msg='test_wedgemc_direct_lambert_beer_top_to_side '//itoa(src))
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
        write(*, FMT='( " diffuse ::  :: ",  8(f10.5) )' ) S
        write(*, FMT='( " target  ::  :: ",  8(f10.5) )' ) S_target
        write(*, FMT='( " diff    ::  :: ",  8(f10.5) )' ) S_target-S
        print*,''
        write(*, FMT='( " direct  ::  :: ", 18(f10.5) )' ) T
        write(*, FMT='( " target  ::  :: ", 18(f10.5) )' ) T_target
        write(*, FMT='( " diff    ::  :: ", 18(f10.5) )' ) T_target-T
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
