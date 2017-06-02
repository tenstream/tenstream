module test_LUT_8_10
  use m_boxmc, only : t_boxmc,t_boxmc_8_10,t_boxmc_1_2,t_boxmc_3_10
  use m_data_parameters, only : mpiint, ireals, iintegers, &
    one, zero, init_mpi_data_parameters, i1, default_str_len
  use m_optprop_LUT, only : t_optprop_LUT_8_10
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse

#include "petsc/finclude/petsc.h"
  use petsc


  use pfunit_mod
  implicit none

  real(ireals) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(8), S_target(10), T_target(8)
  real(ireals) :: S_tol(10),T_tol(8)

  real(ireals) :: BMC_diff2diff(100), BMC_dir2diff(8*10), BMC_dir2dir(8*8)
  real(ireals) :: LUT_diff2diff(100), LUT_dir2diff(8*10), LUT_dir2dir(8*8)

  type(t_boxmc_8_10) :: bmc_8_10
  type(t_optprop_LUT_8_10) :: OPP

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=5e-2, rtol=1e-1

  integer(mpiint) :: ierr

  @testParameter(constructor = newTest)
  type, extends(MpiTestParameter) :: peCase
    real(ireals) :: kabs,ksca,g,phi,theta
  contains
    procedure :: toString
  end type peCase

  @TestCase(constructor = newTest)
  type, extends(MPITestCase) :: parameterized_Test
    real(ireals) :: kabs,ksca,g,phi,theta
  contains
    procedure :: setUp
    procedure :: teardown
  end type parameterized_Test

contains

  ! Constructor for parameter test
  function newTest(testParameter) result (tst)
      type (parameterized_Test) :: tst
      type (peCase), intent(in) :: testParameter

      tst%kabs = testParameter%kabs
      tst%ksca = testParameter%ksca
      tst%g    = testParameter%g
      tst%phi  = testParameter%phi
      tst%theta= testParameter%theta
  end function newTest

  function newPeCase(kabs,ksca,g,phi,theta)
      type (peCase) :: newPeCase
      real(ireals) :: kabs,ksca,g,phi,theta
      newPeCase%kabs  = kabs
      newPeCase%ksca  = ksca
      newPeCase%g     = g
      newPeCase%phi   = phi
      newPeCase%theta = theta
  end function newPeCase

  ! Define the parameters over which to be cycled...
  function getParameters() result(params)
      type(peCase), allocatable :: params(:)

      integer(iintegers) :: ikabs,iksca,ig,iphi,itheta
      real(ireals)       ::  kabs, ksca, g, phi, theta

      integer(iintegers) :: itest,iloop

      do iloop=1,2
        if(iloop.eq.2) allocate(params(itest))
        itest=0

        do ikabs=0,8,5
          do iksca=4,8,2
            do ig=0,5,5
              do iphi=0,0,50
                do itheta=0,5,1

                  itest = itest+1

                  kabs = 10.**(-ikabs)
                  ksca = 10.**(-iksca)
                  g    = ig/10.       
                  phi  = 1.*iphi     
                  theta= 1.*itheta   
                  if(iloop.eq.2) params(itest) = newPeCase(kabs,ksca,g,phi,theta)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
  end function getParameters

  ! Override the parent's version of pECase
  function toString(this) result(string)
      class(pECase), intent(in) :: this
      character(:), allocatable :: string
      allocate(character(default_str_len) :: string)
      write(string,FMT='( 3E8.2, 2I0 )') &
          this%kabs,this%ksca,this%g,int(this%phi),int(this%theta) !,':ranks',this%getNumProcessesRequested()
  end function toString

  @before
  subroutine setup(this)
      class (parameterized_test), intent(inout) :: this

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      !if(myid.eq.0) print *,'Testing LUT coefficients against BOXMC with tolerances atol/rtol :: ',atol,rtol,' :: on ',numnodes,'ranks'

      PETSC_COMM_WORLD = comm
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

      call init_mpi_data_parameters(comm)
      call read_commandline_options()

      call bmc_8_10%init(comm)

      phi   =  0
      theta =  0

      dx = 100
      dy = dx
      dz = 50
  end subroutine setup

  @after
  subroutine teardown(this)
      class (parameterized_test), intent(inout) :: this
      call OPP%destroy()
      call PetscFinalize(ierr)
  end subroutine teardown


  @test( npes=[8,16], testParameters={getParameters()} )
  subroutine test_LUT_direct_coeff(this)
      class (parameterized_test), intent(inout) :: this

      integer(iintegers) :: src
      real(ireals) :: taux, tauz, w0

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      associate( &
            kabs => this%kabs, &
            ksca => this%ksca, &
            g    => this%g,    &
            phi  => this%phi,  &
            theta=> this%theta )

        if(myid.eq.0) print *,'Echo Test for :: ',kabs,ksca,g,phi,theta
        taux = (kabs+ksca) * dx
        tauz = (kabs+ksca) * dz
        w0   = ksca / (kabs+ksca)

        call OPP%init([phi], [theta], comm)

        call OPP%LUT_get_dir2dir (taux, tauz, w0, g , phi, theta, LUT_dir2dir)
        call OPP%LUT_get_dir2diff(taux, tauz, w0, g , phi, theta, LUT_dir2diff)

        do src=1,8

          call bmc_8_10%get_coeff(comm,[kabs,ksca,g],src,.True.,phi,theta,dx,dy,dz,S_target,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)

          ! Rearrange coeffs from dst_ordering to src ordering:
          BMC_dir2diff(src : 8*10 : 8) = S_target
          BMC_dir2dir (src : 8*8  : 8) = T_target
        enddo

        call check(BMC_dir2diff,BMC_dir2dir,LUT_dir2diff,LUT_dir2dir, msg='test_LUT_direct_coeffs')

      end associate
  endsubroutine



  @test( npes=[1,2] )
  subroutine test_LUT_direct_lambert_beer(this)
      !  class (MpiTestMethod), intent(inout) :: this
      class (parameterized_test), intent(inout) :: this

      integer(iintegers) :: src
      real(ireals) :: taux, tauz, w0

      ! direct tests
      bg  = [1e-2, 0., 0. ]
      phi = 0; theta = 0
      S_target = zero

      taux = (bg(1)+bg(2)) * dx
      tauz = (bg(1)+bg(2)) * dz
      w0   = bg(2) / (bg(1)+bg(2))

      call OPP%init([phi], [theta], this%getMpiCommunicator())

      call OPP%LUT_get_dir2dir (taux, tauz, w0, bg(3), phi, theta, LUT_dir2dir)
      call OPP%LUT_get_dir2diff(taux, tauz, w0, bg(3), phi, theta, LUT_dir2diff)

      do src=1,4
        T_target = zero
        T_target(src) = exp(- (bg(1)+bg(2))*dz )

        S = LUT_dir2diff((src-1)*10+1:src*10)
        T = LUT_dir2dir ((src-1)*8+1:src*8)

        T(5:8) = zero ! hard to know that with lambert beer -- use raytracer as test instead

        call check(S_target,T_target, S,T, msg='test_LUT_direct_lambert_beer')
      enddo
  end subroutine


  subroutine check(S_target,T_target, S,T, msg)
      real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

      real(ireals),parameter :: sigma = 6 ! normal test range for coefficients

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

        !print*,'---------------------'
        !write(*, FMT='( " diffuse ::  :: ",10(es12.5) )' ) S
        !write(*, FMT='( " target  ::  :: ",10(es12.5) )' ) S_target
        !write(*, FMT='( " diff    ::  :: ",10(es12.5) )' ) S_target-S
        print*,'RMSE ::: ',rmse(S,S_target)
        !print*,''
        !write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
        !write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
        !write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
        print*,'RMSE ::: ',rmse(T,T_target)
        !print*,'---------------------'
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
