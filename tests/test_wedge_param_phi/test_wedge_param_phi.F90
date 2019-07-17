module test_wedge_param_phi
  use m_boxmc, only : t_boxmc_wedge_18_8
  use m_data_parameters, only : mpiint, iintegers, &
    ireals, irealLUT, ireal_dp, ireal_params, &
    init_mpi_data_parameters, default_str_len, &
    i1, i2, i3, i4, i5, pi_ireal_params
  use m_LUT_param_phi, only: param_phi_from_azimuth, azimuth_from_param_phi, &
    phi_crit, theta_crit, &
    iterative_phi_theta_from_param_phi_and_param_theta, &
    param_phi_param_theta_from_phi_and_theta_withcoords, &
    param_phi_param_theta_from_phi_and_theta_withnormals


  use m_optprop, only : t_optprop_wedge_18_8
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse, CHKERR, get_arg, itoa, &
    ind_nd_to_1d, ind_1d_to_nd, rad2deg, deg2rad, &
    angle_between_two_vec, linspace
  use m_search, only: find_real_location
  use m_boxmc_geometry, only : setup_default_wedge_geometry

#include "petsc/finclude/petsc.h"
  use petsc

  use pfunit_mod
  implicit none

  integer(mpiint) :: myid,mpierr,numnodes,comm, ierr
  real(ireal_params), parameter :: eps=1e-4

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
    real(ireal_params) :: C_wedge(2)

    C_wedge = [.5, .8660254]

    @assertEqual( 60._ireal_params, rad2deg(azimuth_from_param_phi(-2.0_ireal_params, C_wedge)), eps)
    @assertEqual(-60._ireal_params, rad2deg(azimuth_from_param_phi( 2.0_ireal_params, C_wedge)), eps)
    @assertEqual( 30._ireal_params, rad2deg(azimuth_from_param_phi(-1.0_ireal_params, C_wedge)), eps)
    @assertEqual(-30._ireal_params, rad2deg(azimuth_from_param_phi( 1.0_ireal_params, C_wedge)), eps)

    @assertEqual(  0._ireal_params, rad2deg(azimuth_from_param_phi( 0.0_ireal_params, C_wedge)), eps)
    @assertEqual(+45._ireal_params, rad2deg(azimuth_from_param_phi(-1.5_ireal_params, C_wedge)), eps)
    @assertEqual(-45._ireal_params, rad2deg(azimuth_from_param_phi( 1.5_ireal_params, C_wedge)), eps)

    @assertEqual(+15._ireal_params, rad2deg(azimuth_from_param_phi(-0.5_ireal_params, C_wedge)), eps)
    @assertEqual(-15._ireal_params, rad2deg(azimuth_from_param_phi( 0.5_ireal_params, C_wedge)), eps)
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_from_azimuth(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: C_wedge(2)

    C_wedge = [.5, .8660254]

    @assertEqual( -2._ireal_params, param_phi_from_azimuth(deg2rad( 60._ireal_params), C_wedge), eps)
    @assertEqual( -1._ireal_params, param_phi_from_azimuth(deg2rad( 30._ireal_params), C_wedge), eps)
    @assertEqual(  0._ireal_params, param_phi_from_azimuth(deg2rad(  0._ireal_params), C_wedge), eps)
    @assertEqual(  1._ireal_params, param_phi_from_azimuth(deg2rad(-30._ireal_params), C_wedge), eps)
    @assertEqual(  2._ireal_params, param_phi_from_azimuth(deg2rad(-60._ireal_params), C_wedge), eps)

    @assertEqual(+0.5_ireal_params, param_phi_from_azimuth(deg2rad(-15._ireal_params), C_wedge), eps)
    @assertEqual(-0.5_ireal_params, param_phi_from_azimuth(deg2rad(+15._ireal_params), C_wedge), eps)


    C_wedge = [0.511973441, 0.943134367]
    @assertTrue( param_phi_from_azimuth(deg2rad(-32.360667629177840_ireal_params), C_wedge) .gt. 1._ireal_params)
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_from_azimuth_and_back(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: param_phi, phi
    integer(iintegers) :: iphi, iCx, iCy
    real(ireal_params) :: C_wedge(2)
    real(ireal_params), parameter :: eps = sqrt(epsilon(eps))


    do iCx = -100, 100, 5
      do iCy = -20, 20, 5
        C_wedge = [.5+real(iCx)*.01, .8660254+real(iCy)*.01]
        do iphi = -600, 600
          phi = real(iphi, ireal_params)/10
          param_phi = param_phi_from_azimuth(deg2rad(phi), C_wedge)
          @assertEqual(phi, rad2deg(azimuth_from_param_phi( param_phi, C_wedge)), eps)
        enddo
      enddo
    enddo
  end subroutine

  @test(npes=[1])
  subroutine test_phi_crit(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: side_normal(3), theta, phic

    side_normal = [1._ireal_params, 0._ireal_params, 0._ireal_params]
    theta = pi_ireal_params/2
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(0._ireal_params), phic)

    side_normal = [.8660254_ireal_params, -.5_ireal_params, 0._ireal_params]
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(30._ireal_params), phic, eps)

    side_normal = [.8660254_ireal_params, -.5_ireal_params, 0.1_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    theta = deg2rad(10._ireal_params)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertGreaterThan(rad2deg(phic), 60._ireal_params, eps)

    side_normal = [-1._ireal_params, 0._ireal_params, 0.0_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    theta = deg2rad(90._ireal_params)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(deg2rad(0._ireal_params), phic, eps)


    side_normal = [-.8660254_ireal_params, -.5_ireal_params, 0._ireal_params]
    call phi_crit(side_normal, theta, phic, ierr)
    @assertEqual(0_mpiint, ierr)
    @assertEqual(-30._ireal_params, rad2deg(phic), eps)

    side_normal = [-.8660254_ireal_params, -.5_ireal_params, 0.1_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    theta = deg2rad(0._ireal_params)
    call phi_crit(side_normal, theta, phic, ierr)
    @assertFalse(0_mpiint.eq.ierr)
  end subroutine

  @test(npes=[1])
  subroutine test_theta_crit(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: side_normal(3), thetac, phi, theta_target

    side_normal = [0._ireal_params, 1._ireal_params, 0._ireal_params]
    phi = deg2rad(0._ireal_params)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(0._ireal_params, rad2deg(thetac))

    side_normal = [0._ireal_params, 1._ireal_params, 0.1_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    phi = deg2rad(0._ireal_params)
    thetac = theta_crit(side_normal, phi)
    theta_target = angle_between_two_vec(side_normal, [0._ireal_params, 0._ireal_params, -1._ireal_params]) - pi_ireal_params/2
    @assertEqual(rad2deg(theta_target), rad2deg(thetac), eps)

    side_normal = [0._ireal_params, 1._ireal_params, 0.1_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    phi = deg2rad(40._ireal_params)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(7.437376_ireal_params, rad2deg(thetac), eps)

    side_normal = [0._ireal_params, 1._ireal_params, 0.1_ireal_params]
    side_normal = side_normal / norm2(side_normal)
    phi = deg2rad(-40._ireal_params)
    thetac = theta_crit(side_normal, phi)
    @assertEqual(7.437376_ireal_params, rad2deg(thetac), eps)
  end subroutine

  @test(npes=[1])
  subroutine test_iterative_phi_theta_from_param_phi_and_param_theta(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: param_phi, param_theta, phi, theta
    real(ireal_dp), dimension(2) :: A, B, C
    real(ireal_dp) :: dz, sphere_radius
    real(ireal_dp), allocatable :: vertices(:)

    A = [0._ireal_dp, 0._ireal_dp]
    B = [1._ireal_dp, 0._ireal_dp]
    C = [0.5_ireal_dp, 0.8660254_ireal_dp]
    dz = 1
    sphere_radius = 100

    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    param_theta = .01_ireal_params
    param_phi = -1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertGreaterThan(rad2deg(phi), 30._ireal_params)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertLessThan(rad2deg(phi), -30._ireal_params)


    param_theta = 1._ireal_params
    param_phi = -1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(30._ireal_params, rad2deg(phi), eps)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(-30._ireal_params, rad2deg(phi), eps)


    C = [0.511973441_ireal_dp, 0.943134367_ireal_dp]
    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    param_phi = 1
    call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
      param_phi, param_theta, phi, theta, ierr)

    @assertEqual(0_mpiint, ierr)
    @assertEqual(-27.35938568_ireal_params, rad2deg(phi), eps)
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_and_param_theta_back_and_forth(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: param_phi, param_theta, phi, theta, pphi, ptheta
    real(ireal_dp), dimension(2) :: A, B, C
    real(ireal_dp) :: dz, sphere_radius
    real(ireal_dp), allocatable :: vertices(:)
    integer(iintegers), parameter :: Nphi=200, Ntheta=200
    integer(iintegers) :: iphi, itheta
    real(ireal_params), parameter :: eps=1e-3_ireal_params, epsres=1e-4_ireal_params

    A = [0._ireal_dp, 0._ireal_dp]
    B = [1._ireal_dp, 0._ireal_dp]
    C = [0.5_ireal_dp, 0.8660254_ireal_dp]

    dz = 1
    sphere_radius = 6370

    call setup_default_wedge_geometry(A, B, C, dz, vertices, sphere_radius)

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi   = linspace(iphi  , [-0._ireal_params+eps, +1._ireal_params-eps], Nphi)
        param_theta = linspace(itheta, [-0._ireal_params+eps, +1._ireal_params-eps], Ntheta)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
          param_phi, param_theta, phi, theta, ierr)
        @assertEqual(0_mpiint, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, ireal_params), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, epsres)
        @assertEqual(param_theta, ptheta, epsres)
      enddo
    enddo

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi   = linspace(iphi  , [-2._ireal_params+eps, -1._ireal_params-eps], Nphi)
        param_theta = linspace(itheta, [ 0._ireal_params+eps, +1._ireal_params-eps], Ntheta)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
          param_phi, param_theta, phi, theta, ierr)
        @assertEqual(0_mpiint, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, ireal_params), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, epsres)
        @assertEqual(param_theta, ptheta, epsres)
      enddo
    enddo


    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        param_phi   = linspace(iphi  , [ 1._ireal_params+eps, +2._ireal_params-eps], Nphi)
        param_theta = linspace(itheta, [ 0._ireal_params+eps, +1._ireal_params-eps], Ntheta)
        call iterative_phi_theta_from_param_phi_and_param_theta(real(vertices, ireal_params), &
          param_phi, param_theta, phi, theta, ierr)
        @assertEqual(0_mpiint, ierr)

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, ireal_params), &
          phi, theta, pphi, ptheta, ierr)
        !print *,'param_phi/theta', param_phi, param_theta, '=>', rad2deg(phi), rad2deg(theta), '=>', pphi, ptheta

        @assertEqual(0_mpiint, ierr)
        @assertEqual(param_phi, pphi, epsres)
        @assertEqual(param_theta, ptheta, epsres)
      enddo
    enddo
  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_and_param_theta_with_coords_vs_with_normals(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: phi, theta
    real(ireal_params) :: pphi1, ptheta1
    real(ireal_params) :: pphi2, ptheta2
    real(ireal_params), dimension(2) :: A, B, C
    real(ireal_dp) :: dz
    real(ireal_dp), allocatable :: vertices(:)
    integer(iintegers), parameter :: Nphi=200, Ntheta=200
    integer(iintegers) :: iphi, itheta
    real(ireal_params), parameter :: eps=1e-4_ireal_params

    A = [0._ireal_params, 0._ireal_params]
    B = [1._ireal_params, 0._ireal_params]
    C = [0.5_ireal_params, 0.8660254_ireal_params]

    dz = 1

    call setup_default_wedge_geometry(real(A, ireal_dp), real(B, ireal_dp), real(C, ireal_dp), dz, vertices)

    do iphi = 1, Nphi
      do itheta = 1, Ntheta
        phi   = deg2rad(-60._ireal_params + 120._ireal_params / &
          real(Nphi-1, ireal_params) * real(iphi-1, ireal_params))
        theta = deg2rad(0._ireal_params + 90._ireal_params / &
          real(Ntheta-1, ireal_params) * real(itheta-1, ireal_params))

        call param_phi_param_theta_from_phi_and_theta_withcoords(real(vertices, ireal_params), &
          phi, theta, pphi1, ptheta1, ierr)
        @assertEqual(0_mpiint, ierr)

        call param_phi_param_theta_from_phi_and_theta_withnormals(&
          [0._ireal_params, 1._ireal_params, 0._ireal_params], &
          [ C(2), -C(1), 0._ireal_params], &
          [-C(2), -C(1), 0._ireal_params], &
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

  @test(npes=[1])
  subroutine test_param_phi_and_param_theta_with_normals(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params), dimension(3) :: local_normal_base, local_normal_left, local_normal_right
    real(ireal_params) :: Cx, Cy, azimuth, zenith, param_phi, param_theta
    integer(mpiint) :: ierr

    local_normal_base  = [ 7.45058060E-09_ireal_params, 0.998655736_ireal_params, 5.18334582E-02_ireal_params]
    local_normal_left  = [ 0.826860547_ireal_params   ,-0.510045707_ireal_params, 0.236970723_ireal_params]
    local_normal_right = [-0.853828311_ireal_params   ,-0.427395433_ireal_params,-0.297170728_ireal_params]

    Cx = 0.552032232_ireal_params
    Cy = 0.894926786_ireal_params

    azimuth = -1.01413071_ireal_params
    zenith  =  0.298544109_ireal_params

    call param_phi_param_theta_from_phi_and_theta_withnormals(&
      local_normal_base, &
      local_normal_left, &
      local_normal_right,&
      Cx, Cy, &
      azimuth, zenith, &
      param_phi, param_theta, ierr)

    !print *,'test result', param_phi, param_theta, ierr

  end subroutine

  @test(npes=[1])
  subroutine test_param_phi_from_azimuth_and_back_regular_mesh(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_params) :: param_phi, phi
    integer(iintegers) :: iCx, iCy
    real(ireal_params) :: C_wedge(2)
    real(ireal_params), parameter :: eps = sqrt(epsilon(eps))

    do iCx = 1, 3
      do iCy = 1, 3
        C_wedge(1) = linspace(iCx, [0._ireal_dp, 1._ireal_dp], 3_iintegers)
        C_wedge(2) = linspace(iCy, [0.5_ireal_dp, 1.1_ireal_dp], 3_iintegers)
          phi = -1e-3
          param_phi = param_phi_from_azimuth(deg2rad(phi), C_wedge)
          !print *,'C_wedge', C_wedge, 'phi', phi, 'param_phi', param_phi, &
          ! '-->', rad2deg(azimuth_from_param_phi( param_phi, C_wedge))
          @assertEqual(phi, rad2deg(azimuth_from_param_phi( param_phi, C_wedge)), eps)
      enddo
    enddo
  end subroutine
end module
