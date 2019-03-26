module test_wedge_param_phi
  use m_boxmc, only : t_boxmc_wedge_18_8
  use m_data_parameters, only : mpiint, ireals, irealLUT, iintegers, &
    init_mpi_data_parameters, default_str_len, &
    i1, i2, i3, i4, i5, pi_irealLUT
  use m_LUT_param_phi, only: param_phi_from_azimuth, azimuth_from_param_phi, &
    phi_crit, theta_crit, &
    iterative_phi_theta_from_param_phi_and_param_theta, &
    param_phi_param_theta_from_phi_and_theta_withcoords, &
    param_phi_param_theta_from_phi_and_theta_withnormals


  use m_optprop, only : t_optprop_wedge_18_8
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse, CHKERR, get_arg, itoa, &
    ind_nd_to_1d, ind_1d_to_nd, rad2deg, deg2rad
  use m_search, only: find_real_location
  use m_boxmc_geometry, only : setup_default_wedge_geometry

#include "petsc/finclude/petsc.h"
  use petsc

  use pfunit_mod
  implicit none

  integer(mpiint) :: myid,mpierr,numnodes,comm, ierr
  real(irealLUT), parameter :: eps=10*epsilon(eps)

contains

  @before
  subroutine setup(this)
      class (MpiTestMethod), intent(inout) :: this

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      PETSC_COMM_WORLD = comm
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

      call init_mpi_data_parameters(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes=[1])
  subroutine test_azimuth_from_param_phi(this)
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: C_wedge(2)

    C_wedge = [.5, .8660254]

    @assertEqual( 60._irealLUT, rad2deg(azimuth_from_param_phi(-2._irealLUT, C_wedge)))
    @assertEqual(-60._irealLUT, rad2deg(azimuth_from_param_phi( 2._irealLUT, C_wedge)))
    @assertEqual( 30._irealLUT, rad2deg(azimuth_from_param_phi(-1._irealLUT, C_wedge)))
    @assertEqual(-30._irealLUT, rad2deg(azimuth_from_param_phi( 1._irealLUT, C_wedge)))

    @assertEqual(  0._irealLUT, rad2deg(azimuth_from_param_phi( 0._irealLUT, C_wedge)))
    @assertEqual(+45._irealLUT, rad2deg(azimuth_from_param_phi(-1.5_irealLUT, C_wedge)))
    @assertEqual(-45._irealLUT, rad2deg(azimuth_from_param_phi( 1.5_irealLUT, C_wedge)))

    @assertEqual(+15._irealLUT, rad2deg(azimuth_from_param_phi(-0.5_irealLUT, C_wedge)))
    @assertEqual(-15._irealLUT, rad2deg(azimuth_from_param_phi( 0.5_irealLUT, C_wedge)))
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_from_azimuth(this)
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: C_wedge(2)

    C_wedge = [.5, .8660254]

    @assertEqual( -2._irealLUT, param_phi_from_azimuth(deg2rad( 60._irealLUT), C_wedge), eps)
    @assertEqual( -1._irealLUT, param_phi_from_azimuth(deg2rad( 30._irealLUT), C_wedge), eps)
    @assertEqual(  0._irealLUT, param_phi_from_azimuth(deg2rad(  0._irealLUT), C_wedge), eps)
    @assertEqual(  1._irealLUT, param_phi_from_azimuth(deg2rad(-30._irealLUT), C_wedge), eps)
    @assertEqual(  2._irealLUT, param_phi_from_azimuth(deg2rad(-60._irealLUT), C_wedge), eps)

    @assertEqual(+0.5_irealLUT, param_phi_from_azimuth(deg2rad(-15._irealLUT), C_wedge), eps)
    @assertEqual(-0.5_irealLUT, param_phi_from_azimuth(deg2rad(+15._irealLUT), C_wedge), eps)


    C_wedge = [0.511973441, 0.943134367]
    @assertTrue( param_phi_from_azimuth(deg2rad(-32.360667629177840_irealLUT), C_wedge) .gt. 1._irealLUT)
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_from_azimuth_and_back(this)
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: param_phi, phi
    integer(iintegers) :: iphi, iCx, iCy
    real(irealLUT) :: C_wedge(2)
    real(irealLUT), parameter :: eps = sqrt(epsilon(eps))


    do iCx = -100, 100, 5
      do iCy = -20, 20, 5
        C_wedge = [.5+iCx*.01, .8660254+iCy*.01]
        do iphi = -600, 600
          phi = real(iphi, irealLUT)/10
          param_phi = param_phi_from_azimuth(deg2rad(phi), C_wedge)
          @assertEqual(phi, rad2deg(azimuth_from_param_phi( param_phi, C_wedge)), eps)
        enddo
      enddo
    enddo
  end subroutine

  @test(npes=[1])
  subroutine test_phi_crit(this)
    use m_helper_functions, only: norm, angle_between_two_vec
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: side_normal(3), theta, phic

    side_normal = [1._irealLUT, 0._irealLUT, 0._irealLUT]
    theta = pi_irealLUT/2
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(0._irealLUT), phic)

    side_normal = [.8660254_irealLUT, -.5_irealLUT, 0._irealLUT]
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(30._irealLUT), phic)

    side_normal = [.8660254_irealLUT, -.5_irealLUT, 0.1_irealLUT]
    side_normal = side_normal / norm(side_normal)
    theta = deg2rad(10._irealLUT)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertGreaterThan(rad2deg(phic), 60._irealLUT)

    side_normal = [-1._irealLUT, 0._irealLUT, 0.0_irealLUT]
    side_normal = side_normal / norm(side_normal)
    theta = deg2rad(90._irealLUT)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(0._irealLUT), phic)


    side_normal = [-.8660254_irealLUT, -.5_irealLUT, 0._irealLUT]
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(-30._irealLUT), phic)

    side_normal = [-.8660254_irealLUT, -.5_irealLUT, 0.1_irealLUT]
    side_normal = side_normal / norm(side_normal)
    theta = deg2rad(0._irealLUT)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertFalse(0_mpiint.eq.ierr)
    @assertEqual(deg2rad(-30._irealLUT), phic)
  end subroutine

  @test(npes=[1])
  subroutine test_theta_crit(this)
    use m_helper_functions, only: norm, angle_between_two_vec
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: side_normal(3), thetac, phi

    side_normal = [0._irealLUT, 1._irealLUT, 0._irealLUT]
    phi = deg2rad(0._irealLUT)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(0._irealLUT, rad2deg(thetac))

    side_normal = [0._irealLUT, 1._irealLUT, 0.1_irealLUT]
    side_normal = side_normal / norm(side_normal)
    phi = deg2rad(0._irealLUT)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(5.47000599_irealLUT, rad2deg(thetac), eps)

    side_normal = [0._irealLUT, 1._irealLUT, 0.1_irealLUT]
    side_normal = side_normal / norm(side_normal)
    phi = deg2rad(40._irealLUT)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(6.92789030_irealLUT, rad2deg(thetac), eps)

    side_normal = [0._irealLUT, 1._irealLUT, 0.1_irealLUT]
    side_normal = side_normal / norm(side_normal)
    phi = deg2rad(-40._irealLUT)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(6.92789030_irealLUT, rad2deg(thetac), eps)
  end subroutine

  @test(npes=[1])
  subroutine test_iterative_phi_theta_from_param_phi_and_param_theta(this)
    use m_helper_functions, only: norm, angle_between_two_vec
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: param_phi, param_theta, phi, theta
    real(ireals), dimension(2) :: A, B, C
    real(ireals) :: dz, sphere_radius
    real(ireals), allocatable :: vertices(:)
    real(irealLUT), parameter :: eps=sqrt(epsilon(eps))

    A = [0._ireals, 0._ireals]
    B = [1._ireals, 0._ireals]
    C = [0.5_ireals, 0.8660254_ireals]
    dz = 1
    sphere_radius = 100

    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    param_theta = .01_irealLUT
    param_phi = -1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertGreaterThan(rad2deg(phi), 30._irealLUT)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertLessThan(rad2deg(phi), -30._irealLUT)


    param_theta = 1._irealLUT
    param_phi = -1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(30._irealLUT, rad2deg(phi), eps)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(-30._irealLUT, rad2deg(phi), eps)


    C = [0.511973441_ireals, 0.943134367_ireals]
    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(-27.3593781_irealLUT, rad2deg(phi), eps)
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_and_param_theta_back_and_forth(this)
    use m_helper_functions, only: norm, angle_between_two_vec
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: param_phi, param_theta, phi, theta, pphi, ptheta
    real(ireals), dimension(2) :: A, B, C
    real(ireals) :: dz, sphere_radius
    real(ireals), allocatable :: vertices(:)
    integer(iintegers), parameter :: Nphi=200, Ntheta=200
    integer(iintegers) :: iphi, itheta
    real(irealLUT), parameter :: eps=1e-4_irealLUT

    A = [0._ireals, 0._ireals]
    B = [1._ireals, 0._ireals]
    C = [0.5_ireals, 0.8660254_ireals]

    dz = 1
    sphere_radius = 100

    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi = -1._irealLUT + 2._irealLUT / &
          real(Nphi-1, irealLUT) * real(iphi-1, irealLUT)
        param_theta = -1._irealLUT + 2._irealLUT / &
          real(Ntheta-1, irealLUT) * real(itheta-1, irealLUT)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
          param_phi, param_theta, phi, theta, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, irealLUT), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, eps)
        @assertEqual(param_theta, ptheta, eps)
      enddo
    enddo

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi = -2._irealLUT + 1._irealLUT / &
          real(Nphi-1, irealLUT) * real(iphi-1, irealLUT)
        param_theta = 1e-4_irealLUT + (1._irealLUT-1e-4_irealLUT) / &
          real(Ntheta-1, irealLUT) * real(itheta-1, irealLUT)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
          param_phi, param_theta, phi, theta, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, irealLUT), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, eps)
        @assertEqual(param_theta, ptheta, eps)
      enddo
    enddo


    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi = 1._irealLUT + 1._irealLUT / &
          real(Nphi-1, irealLUT) * real(iphi-1, irealLUT)
        param_theta = 1e-4_irealLUT + (1._irealLUT-1e-4_irealLUT) / &
          real(Ntheta-1, irealLUT) * real(itheta-1, irealLUT)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, irealLUT), &
          param_phi, param_theta, phi, theta, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, irealLUT), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, eps)
        @assertEqual(param_theta, ptheta, eps)
      enddo
    enddo
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_and_param_theta_with_coords_vs_with_normals(this)
    use m_helper_functions, only: norm, angle_between_two_vec
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: phi, theta
    real(irealLUT) :: pphi1, ptheta1
    real(irealLUT) :: pphi2, ptheta2
    real(irealLUT), dimension(2) :: A, B, C
    real(ireals) :: dz
    real(ireals), allocatable :: vertices(:)
    integer(iintegers), parameter :: Nphi=200, Ntheta=200
    integer(iintegers) :: iphi, itheta
    real(irealLUT), parameter :: eps=1e-4_irealLUT

    A = [0._irealLUT, 0._irealLUT]
    B = [1._irealLUT, 0._irealLUT]
    C = [0.5_irealLUT, 0.8660254_irealLUT]

    dz = 1

    call setup_default_wedge_geometry(real(A, ireals), real(B, ireals), real(C, ireals), dz, vertices)

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        phi   = deg2rad(-60._irealLUT + 120._irealLUT / &
          real(Nphi-1, irealLUT) * real(iphi-1, irealLUT))
        theta = deg2rad(0._irealLUT + 90._irealLUT / &
          real(Ntheta-1, irealLUT) * real(itheta-1, irealLUT))

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, irealLUT), &
          phi, theta, pphi1, ptheta1, ierr)
        @assertEqual(0_mpiint, ierr)

        call param_phi_param_theta_from_phi_and_theta_withnormals(&
          [0._irealLUT, 1._irealLUT, 0._irealLUT], &
          [ C(2), -C(1), 0._irealLUT], &
          [-C(2), -C(1), 0._irealLUT], &
          C(1), C(2), &
          phi, theta, pphi2, ptheta2, ierr)

        !print *,'phi, theta', phi, theta
        !print *,'pphi1', pphi1, 'ptheta1', ptheta1
        !print *,'pphi2', pphi2, 'ptheta2', ptheta2

        @assertEqual(0_mpiint, ierr)
        @assertEqual(pphi1, pphi2, eps)
        @assertEqual(ptheta1, ptheta2, eps)
      enddo
    enddo
  end subroutine
end module
