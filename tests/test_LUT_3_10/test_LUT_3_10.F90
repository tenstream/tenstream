module test_LUT_3_10
  use m_boxmc, only: t_boxmc_3_10
  use m_data_parameters, only: mpiint, iintegers, &
                               ireals, irealLUT, ireal_dp, &
                               init_mpi_data_parameters, i1, default_str_len
  use m_optprop_LUT, only: t_optprop_LUT_3_10
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

#include "petsc/finclude/petsc.h"
  use petsc

  use pfunit_mod
  implicit none

  integer(iintegers), parameter :: Ndir = 3, Ndiff = 10

  real(irealLUT) :: bg(3), phi, theta, dx, dy, dz
  real(irealLUT) :: S(Ndiff), T(Ndir)
  real(ireals) :: S_target(Ndiff), T_target(Ndir)
  real(ireals) :: S_tol(Ndiff), T_tol(Ndir)
  real(ireals), allocatable :: vertices(:)

  real(ireals) :: BMC_diff2diff(Ndiff * Ndiff), BMC_dir2diff(Ndir * Ndiff), BMC_dir2dir(Ndir * Ndir)
  real(irealLUT) :: LUT_diff2diff(Ndiff * Ndiff), LUT_dir2diff(Ndir * Ndiff), LUT_dir2dir(Ndir * Ndir)

  type(t_boxmc_3_10) :: bmc_3_10
  type(t_optprop_LUT_3_10) :: OPP

  integer(mpiint) :: myid, mpierr, numnodes, comm

  real(irealLUT), parameter :: atol = 5e-2, rtol = 1e-1, zero = 0, one = 1

  integer(mpiint) :: ierr

  @testParameter(constructor=newTest)
  type, extends(MpiTestParameter) :: peCase
    real(irealLUT) :: kabs, ksca, g, phi, theta
  contains
    procedure :: toString
  end type peCase

  @TestCase(constructor=newTest)
  type, extends(MPITestCase) :: parameterized_Test
    real(irealLUT) :: kabs, ksca, g, phi, theta
  contains
    procedure :: setUp
    procedure :: teardown
  end type parameterized_Test

contains

  ! Constructor for parameter test
  function newTest(testParameter) result(tst)
    type(parameterized_Test) :: tst
    type(peCase), intent(in) :: testParameter

    tst%kabs = testParameter%kabs
    tst%ksca = testParameter%ksca
    tst%g = testParameter%g
    tst%phi = testParameter%phi
    tst%theta = testParameter%theta
  end function newTest

  function newPeCase(kabs, ksca, g, phi, theta)
    type(peCase) :: newPeCase
    real(irealLUT) :: kabs, ksca, g, phi, theta
    newPeCase%kabs = kabs
    newPeCase%ksca = ksca
    newPeCase%g = g
    newPeCase%phi = phi
    newPeCase%theta = theta
  end function newPeCase

  ! Define the parameters over which to be cycled...
  function getParameters() result(params)
    use m_optprop_parameters, only: preset_tau31, preset_w010, preset_g3

    type(peCase), allocatable :: params(:)

    integer(iintegers) :: itau, iw0, ig, itheta, iphi
    real(irealLUT) :: kabs, ksca, g, phi, theta

    integer(iintegers) :: itest, iloop
    integer(iintegers) :: Ntau, Nw0, Ng

    associate (preset_tau => preset_tau31, &
               preset_w0 => preset_w010, &
               preset_g => preset_g3)

      Ntau = size(preset_tau)
      Nw0 = size(preset_w0)
      Ng = size(preset_g)

      do iloop = 1, 2
        if (iloop .eq. 2) allocate (params(itest * 2))
        itest = 0

        do itau = 1, Ntau - 1, Ntau / 2
          do iw0 = 1, Nw0 - 1, Nw0 / 2
            do ig = 1, Ng - 1, Ng / 2
              do iphi = 1, 85, 40
                do itheta = 1, 85, 40

                  itest = itest + 1

                  ! Do Lookup Tests directly on supports of LUT
                  kabs = preset_tau(itau) * (one - preset_w0(iw0))
                  ksca = preset_tau(itau) * preset_w0(iw0)
                  g = preset_g(ig)
                  phi = real(iphi, irealLUT)
                  theta = real(itheta, irealLUT)
                  if (iloop .eq. 2) params(2 * (itest - 1) + 1) = newPeCase(kabs, ksca, g, phi, theta)

                  ! Again, do interlaced tests between support points of LUT
                  kabs = (preset_tau(itau) + preset_tau(itau + 1)) / 2.*(one - (preset_w0(iw0) + preset_w0(iw0 + 1)) / 2.)
                  ksca = (preset_tau(itau) + preset_tau(itau + 1)) / 2.*(preset_w0(iw0) + preset_w0(iw0 + 1)) / 2.
                  g = (preset_g(ig) + preset_g(ig + 1)) / 2.
                  if (iloop .eq. 2) params(2 * (itest - 1) + 2) = newPeCase(kabs, ksca, g, phi, theta)
                end do
              end do
            end do
          end do
        end do
      end do
    end associate
  end function getParameters

  ! Override the parent's version of pECase
  function toString(this) result(string)
    class(pECase), intent(in) :: this
    character(:), allocatable :: string
    allocate (character(default_str_len) :: string)
    write (string, FMT='( 3E9.2, 2I0 )') &
      this%kabs, this%ksca, this%g, int(this%phi), int(this%theta) !,':ranks',this%getNumProcessesRequested()
  end function toString

  @before
  subroutine setup(this)
    class(parameterized_test), intent(inout) :: this

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    !if(myid.eq.0) print *,'Testing LUT coefficients against BOXMC with tolerances atol/rtol :: ',atol,rtol,' :: on ',numnodes,'ranks'

    PETSC_COMM_WORLD = comm
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)

    call bmc_3_10%init(comm)

    phi = 0
    theta = 0

    dx = 100
    dy = dx
    dz = 50

    call setup_default_unit_cube_geometry(real(dx, ireals), real(dy, ireals), real(dz, ireals), vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class(parameterized_test), intent(inout) :: this
    !call OPP%destroy()
    call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes=[2, 1], testParameters={getParameters() })
  subroutine test_LUT_direct_coeff(this)
    class(parameterized_test), intent(inout) :: this

    integer(iintegers) :: src
    real(irealLUT) :: taux, tauz, w0

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    associate ( &
      kabs => this%kabs, &
      ksca => this%ksca, &
      g => this%g, &
      phi => this%phi, &
      theta => this%theta)

      if (myid .eq. 0) print *, 'Echo Test for :: ', kabs, ksca, g, phi, theta
      taux = (kabs + ksca) * dx
      tauz = (kabs + ksca) * dz
      w0 = ksca / (kabs + ksca)

      call OPP%init(comm)

      call OPP%get_dir2dir([tauz, w0, tauz / taux, g, phi, theta], LUT_dir2dir)
      call OPP%get_dir2diff([tauz, w0, tauz / taux, g, phi, theta], LUT_dir2diff)
      print *, taux, tauz
      do src = 1, Ndir

        call bmc_3_10%get_coeff(comm, &
                                real([kabs, ksca, g], ireal_dp), &
                                src, .true., &
                                real(phi, ireal_dp), &
                                real(theta, ireal_dp), &
                                real(vertices, ireal_dp), &
                                S_target, T_target, S_tol, T_tol, &
                                inp_atol=real(atol, ireal_dp), inp_rtol=real(rtol, ireal_dp))

        ! Rearrange coeffs from dst_ordering to src ordering:
        BMC_dir2diff(src:Ndir * Ndiff:Ndir) = S_target
        BMC_dir2dir(src:Ndir * Ndir:Ndir) = T_target
      end do

      call check(BMC_dir2diff, BMC_dir2dir, LUT_dir2diff, LUT_dir2dir, msg='test_LUT_direct_coeffs')

    end associate
  end subroutine

  @test(npes=[2, 1], testParameters={getParameters() })
  subroutine test_LUT_diff_coeff(this)
    class(parameterized_test), intent(inout) :: this

    integer(iintegers) :: src
    real(irealLUT) :: taux, tauz, w0

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    associate ( &
      kabs => this%kabs, &
      ksca => this%ksca, &
      g => this%g, &
      phi => this%phi, &
      theta => this%theta)

      if (myid .eq. 0) print *, 'Echo Test for :: ', kabs, ksca, g, phi, theta
      taux = (kabs + ksca) * dx
      tauz = (kabs + ksca) * dz
      w0 = ksca / (kabs + ksca)

      call OPP%init(comm)

      call OPP%get_diff2diff([tauz, w0, tauz / taux, g], LUT_diff2diff)
      do src = 1, Ndiff

        call bmc_3_10%get_coeff(comm, &
                                real([kabs, ksca, g], ireal_dp), &
                                src, .false., &
                                real(phi, ireal_dp), &
                                real(theta, ireal_dp), &
                                real(vertices, ireal_dp), &
                                S_target, T_target, S_tol, T_tol, &
                                inp_atol=real(atol, ireal_dp), inp_rtol=real(rtol, ireal_dp))

        ! Rearrange coeffs from dst_ordering to src ordering:
        BMC_diff2diff(src:Ndiff * Ndiff:Ndiff) = S_target
      end do
      BMC_dir2dir = zero
      LUT_dir2dir = zero

      call check(BMC_diff2diff, BMC_dir2dir, LUT_diff2diff, LUT_dir2dir, msg='test_LUT_diff_coeffs')

    end associate
  end subroutine

  @test(npes=[1])
  subroutine test_LUT_direct_lambert_beer(this)
    class(parameterized_test), intent(inout) :: this

    integer(iintegers) :: src
    real(irealLUT) :: taux, tauz, w0

    ! direct tests
    bg = [1e-2, 0., 0.]
    phi = 0; theta = 0
    S_target = zero

    taux = (bg(1) + bg(2)) * dx
    tauz = (bg(1) + bg(2)) * dz
    w0 = bg(2) / (bg(1) + bg(2))

    call OPP%init(this%getMpiCommunicator())

    call OPP%get_dir2dir([tauz, w0, tauz / taux, bg(3), phi, theta], LUT_dir2dir)
    call OPP%get_dir2diff([tauz, w0, tauz / taux, bg(3), phi, theta], LUT_dir2diff)

    do src = 1, 1
      T_target = zero
      T_target(src) = exp(-(bg(1) + bg(2)) * dz)

      S = LUT_dir2diff((src - 1) * Ndiff + 1:src * Ndiff)
      T = LUT_dir2dir((src - 1) * Ndir + 1:src * Ndir)

      T(2:Ndir) = zero ! hard to know that with lambert beer -- use raytracer as test instead

      call check(S_target, T_target, S, T, msg='test_LUT_direct_lambert_beer')
    end do
  end subroutine

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target
    real(irealLUT), intent(in), dimension(:) :: S, T

    real(irealLUT), parameter :: sigma = 6 ! normal test range for coefficients

    character(len=*), optional :: msg
    character(default_str_len) :: local_msgS, local_msgT
    logical, parameter :: ldetail = .true.

    if (myid .eq. 0) then
      if (ldetail) then
        print *, ''

        if (present(msg)) then
          write (local_msgS, *) trim(msg), ':: Diffuse boxmc coefficient not as '
          write (local_msgT, *) trim(msg), ':: Direct  boxmc coefficient not as '
          print *, msg
        else
          write (local_msgS, *) 'Diffuse boxmc coefficient not as '
          write (local_msgT, *) 'Direct  boxmc coefficient not as '
        end if

        print *, '---------------------'
        write (*, FMT='( " diffuse :: ",10(f12.5) )') S
        print *, ''
        write (*, FMT='( " target  :: ",10(f12.5) )') S_target
        print *, ''
        write (*, FMT='( " diff    :: ",10(f12.5) )') S_target - S
        print *, 'RMSE ::: ', rmse(S, real(S_target, irealLUT))*[one, 100._ireallut], '%'
        print *, ''
        write (*, FMT='( " direct  :: ", 8(f12.5) )') T
        print *, ''
        write (*, FMT='( " target  :: ", 8(f12.5) )') T_target
        print *, ''
        write (*, FMT='( " diff    :: ", 8(f12.5) )') T_target - T
        print *, 'RMSE ::: ', rmse(T, real(T_target, irealLUT))*[one, 100._ireallut], '%'
        print *, '---------------------'
        print *, ''
      else
        print *, 'RMSE ::: ', rmse(S, real(S_target, irealLUT))*[one, 100._ireallut], '% ', &
          'direct:::', rmse(T, real(T_target, irealLUT))*[one, 100._ireallut], '%'
        print *, ''
      end if

      @assertEqual(real(S_target, irealLUT), S, atol * sigma, local_msgS)
      @assertLessThanOrEqual(zero, S)
      @assertGreaterThanOrEqual(one, S)

      @assertEqual(real(T_target, irealLUT), T, atol * sigma, local_msgT)
      @assertLessThanOrEqual(zero, T)
      @assertGreaterThanOrEqual(one, T)
    end if
  end subroutine

end module
