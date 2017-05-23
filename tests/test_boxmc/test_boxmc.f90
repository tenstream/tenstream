module test_boxmc
  use m_boxmc, only : t_boxmc,t_boxmc_8_10,t_boxmc_1_2,t_boxmc_3_10
  use m_data_parameters, only :     &
    mpiint, ireals, iintegers,      &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol

  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(8), S_target(10), T_target(8)
  real(ireals) :: S_tol(10),T_tol(8)

  type(t_boxmc_8_10) :: bmc_8_10

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-2, rtol=1e-1
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
      theta = 45

      dx = 100
      dy = dx
      dz = 50

      ! computed target with stddev_atol=5e-6, stddev_rtol=1e-4 in optprop_parameters
      ! inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) !
      !    call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol,inp_atol=1e-6_ireals, inp_rtol=1e-4_ireals) ! inp_atol=atol, inp_rtol=rtol)

  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1,8])
  subroutine test_boxmc_select_cases_direct(this)
      class (MpiTestMethod), intent(inout) :: this

      ! direct tests
      bg  = [1e-3, 1e-3, 0. ]
      phi = 0; theta = 0
      S_target = [1.16098E-02,1.13492E-02,4.69646E-03,1.11393E-03,4.61549E-03,1.08210E-03,4.70091E-03,1.11544E-03,4.61404E-03,1.08313E-03]
      T_target = [9.04842E-01,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00]

      call bmc_8_10%get_coeff(comm,bg,1,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct')
  end subroutine

  @test(npes =[1,2])
  subroutine test_boxmc_direct_lambert_beer(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: src

      ! direct tests
      bg  = [1e-2, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      do src=1,4
        T_target = zero
        T_target(src) = exp(- (bg(1)+bg(2))*dz )

        call bmc_8_10%get_coeff(comm,bg,src,.True.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_boxmc_direct_lambert_beer')
      enddo
  end subroutine

  @test(npes =[1,8,16])
  subroutine test_boxmc_select_cases_diffuse(this)
      class (MpiTestMethod), intent(inout) :: this

      bg = [1e-3, 1e-16, .0 ]
      S_target = [0.00000E+00,3.90217E-01,1.40431E-01,1.40417E-01,0.00000E+00,0.00000E+00,1.40416E-01,1.40429E-01,0.00000E+00,0.00000E+00]
      T_target = zero

      call bmc_8_10%get_coeff(comm,bg,1,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_diffuse')
! --------------------------------------------

      bg = [1e-3, 1e-4, .5 ]
      S_target = [5.91803E-04,3.89499E-01,1.40297E-01,1.40317E-01,1.44767E-04,1.45072E-04,1.40312E-01,1.40306E-01,1.45083E-04,1.44973E-04]

      call bmc_8_10%get_coeff(comm,bg,1,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_diffuse')
! --------------------------------------------

      bg = [1e-3, 1e-4, .999999 ]
      S_target = [1.85844E-09,3.90217E-01,1.40431E-01,1.40429E-01,0.00000E+00,0.00000E+00,1.40417E-01,1.40417E-01,0.00000E+00,9.64722E-10]

      call bmc_8_10%get_coeff(comm,bg,1,.False.,phi,theta,dx,dy,dz,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_diffuse')
  end subroutine

  subroutine check(S_target,T_target, S,T, msg)
      real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

      real(ireals),parameter :: sigma = 3 ! normal test range for coefficients

      integer(iintegers) :: i
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
