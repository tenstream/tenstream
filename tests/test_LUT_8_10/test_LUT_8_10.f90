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
    use m_optprop_parameters, only: Ntheta, Nphi, &
      preset_tau21, preset_w015, preset_g3

      type(peCase), allocatable :: params(:)

      integer(iintegers) :: itau, iw0, ig, itheta, iphi
      real(ireals)       ::  kabs, ksca, g, phi, theta

      integer(iintegers) :: itest,iloop
      integer(iintegers) :: Ntau, Nw0, Ng

      associate ( preset_tau => preset_tau21, &
          preset_w0  => preset_w015, &
          preset_g   => preset_g3 )

        Ntau = size(preset_tau)
        Nw0  = size(preset_w0)
        Ng   = size(preset_g)

        do iloop=1,2
          if(iloop.eq.2) allocate(params(itest*2))
          itest=0

          do itau = 1, Ntau-1, Ntau/2
            do iw0 = 1, Nw0-1, Nw0/2
              do ig = 1, Ng-1, Ng/2
                do iphi = 1, 85, 40
                  do itheta = 1, 85, 40

                    itest = itest+1

                    ! Do Lookup Tests directly on supports of LUT
                    kabs = preset_tau(itau) * (one-preset_w0(iw0))
                    ksca = preset_tau(itau) * preset_w0(iw0)
                    g    = preset_g(ig)
                    phi  = real(iphi, ireals)
                    theta= real(itheta, ireals)
                    if(iloop.eq.2) params(2*(itest-1)+1) = newPeCase(kabs,ksca,g,phi,theta)

                    ! Again, do interlaced tests between support points of LUT
                    kabs = (preset_tau(itau)+preset_tau(itau+1))/2. * (one-(preset_w0(iw0)+preset_w0(iw0+1))/2.)
                    ksca = (preset_tau(itau)+preset_tau(itau+1))/2. * (preset_w0(iw0)+preset_w0(iw0+1))/2.
                    g    = (preset_g(ig)+preset_g(ig+1))/2.
                    if(iloop.eq.2) params(2*(itest-1)+2) = newPeCase(kabs,ksca,g,phi,theta)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      end associate
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
      !call OPP%destroy()
      call PetscFinalize(ierr)
  end subroutine teardown


  @test( npes=[2,1], testParameters={getParameters()} )
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

        call OPP%init(comm)

        call OPP%LUT_get_dir2dir ([tauz, w0, g , tauz/taux, phi, theta], LUT_dir2dir)
        call OPP%LUT_get_dir2diff([tauz, w0, g , tauz/taux, phi, theta], LUT_dir2diff)
        print*,taux, tauz
        do src=1,8

          call bmc_8_10%get_coeff(comm,[kabs,ksca,g],src,.True.,phi,theta,dx,dy,dz,S_target,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)

          ! Rearrange coeffs from dst_ordering to src ordering:
          BMC_dir2diff(src : 8*10 : 8) = S_target
          BMC_dir2dir (src : 8*8  : 8) = T_target
        enddo

        call check(BMC_dir2diff,BMC_dir2dir,LUT_dir2diff,LUT_dir2dir, msg='test_LUT_direct_coeffs')

      end associate
  endsubroutine

  @test( npes=[2,1], testParameters={getParameters()} )
  subroutine test_LUT_diff_coeff(this)
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

        call OPP%init(comm)

        call OPP%LUT_get_diff2diff([tauz, w0, g, tauz/taux], LUT_diff2diff)
        do src=1,10

          call bmc_8_10%get_coeff(comm,[kabs,ksca,g],src,.False.,phi,theta,dx,dy,dz,S_target,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)

          ! Rearrange coeffs from dst_ordering to src ordering:
          BMC_diff2diff(src : 10*10 : 10) = S_target
        enddo
        BMC_dir2dir = zero
        LUT_dir2dir = zero

        call check(BMC_diff2diff,BMC_dir2dir,LUT_diff2diff,LUT_dir2dir, msg='test_LUT_diff_coeffs')

      end associate
  endsubroutine



  @test( npes=[1] )
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

      call OPP%init(this%getMpiCommunicator())

      call OPP%LUT_get_dir2dir ([tauz, w0, bg(3), tauz/taux, phi, theta], LUT_dir2dir)
      call OPP%LUT_get_dir2diff([tauz, w0, bg(3), tauz/taux, phi, theta], LUT_dir2diff)

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

      character(len=*),optional :: msg
      character(default_str_len) :: local_msgS, local_msgT
      logical, parameter :: ldetail=.True.

      if(myid.eq.0) then
        if(ldetail) then
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
          write(*, FMT='( " diffuse :: ",10(f12.5) )' ) S
          print*,''
          write(*, FMT='( " target  :: ",10(f12.5) )' ) S_target
          print*,''
          write(*, FMT='( " diff    :: ",10(f12.5) )' ) S_target-S
          print*,'RMSE ::: ',rmse(S,S_target)*[one, 100._ireals],'%'
          print*,''
          write(*, FMT='( " direct  :: ", 8(f12.5) )' ) T
          print*,''
          write(*, FMT='( " target  :: ", 8(f12.5) )' ) T_target
          print*,''
          write(*, FMT='( " diff    :: ", 8(f12.5) )' ) T_target-T
          print*,'RMSE ::: ',rmse(T,T_target)*[one, 100._ireals],'%'
          print*,'---------------------'
          print*,''
        else
          print*,'RMSE ::: ',rmse(S,S_target)*[one, 100._ireals],'% ', &
                 'direct:::',rmse(T,T_target)*[one, 100._ireals],'%'
          print *,''
        endif

        @assertEqual(S_target, S, atol*sigma, local_msgS )
        @assertLessThanOrEqual   (zero, S)
        @assertGreaterThanOrEqual(one , S)

        @assertEqual(T_target, T, atol*sigma, local_msgT )
        @assertLessThanOrEqual   (zero, T)
        @assertGreaterThanOrEqual(one , T)
      endif
  end subroutine

end module
