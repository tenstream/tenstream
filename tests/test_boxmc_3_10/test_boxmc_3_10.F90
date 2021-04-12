module test_boxmc_3_10
  use m_boxmc, only : t_boxmc, t_boxmc_3_10
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, irealLUT, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_geometric_coeffs, only : dir2dir3_geometric_coeff_corr
  use m_helper_functions, only : spherical_2_cartesian, cstr, toStr

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(3), S_target(10), T_target(3)
  real(ireals) :: S_tol(10),T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_10) :: bmc_3_10

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireal_dp),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-1
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()


    call init_mpi_data_parameters(comm)

    call bmc_3_10%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if(myid.eq.0) print *,numnodes,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 50
    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  !@test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp/2 ]

    ! from top to bot face
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_direct_srctopface top_to_bot')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srctopface_45(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp/2 ]

    ! outwards from face 2
    phi = 0; theta = 45*1._ireal_dp

    tau = (bg(1)+bg(2)) * dz * sqrt(2*1._ireal_dp)

    S_target = zero

    T_target = zero
    T_target(1) = real(exp(-tau)/2, ireals)
    T_target(3) = real((1-exp(-tau))/(2*tau), ireals)

    src = 1

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_direct_srctopface_45')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srcsidefaces(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp/2 ]


    ! along each of the side faces
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real((sinh(tau)-cosh(tau)+1)/tau, ireals)
    !von der Seite nach unten
    ! ???????? HERE

    do iphi=0,360,30
      phi = real(iphi, ireal_dp)
      do src = 2,3 ! Seitenflächen
        call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target,T_target, S,T, msg='test_boxmc_select_cases_direct_srcsidefaces')
      enddo
    enddo
  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp ]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    S_target = [0.0, 0.390156, 0.1404375, 0.1404375, 0.0, 0.0, 0.1404375, 0.1404375, 0.0, 0.0]

    T_target = zero


    src = 2   !top face
    call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srctopface')
  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srcbottomface(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp ]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1)+bg(2)) * dz

    T_target = zero

    S_target = [0.390156, 0.0, 0.0, 0.0, 0.1404375, 0.1404375, 0.0, 0.0, 0.1404375, 0.1404375]

    src = 1
    call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcbottomface')
  end subroutine


  !@test(npes =[1])
  subroutine test_boxmc_select_cases_diff_srcsideface(this)
     class (MpiTestMethod), intent(inout) :: this
     integer(iintegers) :: src
     real(ireal_dp) :: tau
     real,parameter :: top=0.56173, a=0.104806, b=0.1424402

     ! direct to diffuse tests
     bg  = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp ]

     ! outwards from face 2
     phi = 0; theta = 0

     tau = (bg(1)+bg(2)) * dz
     T_target = zero

     src = 3
     S_target = [0.0, top, a, 0.0, 0.0, 0.0, b, b, 0.0, 0.0]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 3')

     src = 4
     S_target = [0.0, top, 0.0, a, 0.0, 0.0, b, b, 0.0, 0.0]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 4')

     src = 5
     S_target = [top, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, b, b]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 5')

     src = 6
     S_target = [top, 0.0, 0.0, 0.0, 0.0, a, 0.0, 0.0, b, b]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 6')

     src = 7
     S_target = [0.0, top, b, b, 0.0, 0.0, a, 0.0, 0.0, 0.0]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 7')

     src = 8
     S_target = [0.0, top, b, b, 0.0, 0.0, 0.0, a, 0.0, 0.0]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 8')

     src = 9
     S_target = [top, 0.0, 0.0, 0.0, b, b, 0.0, 0.0, a, 0.0]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 9')

     src = 10
     S_target = [top, 0.0, 0.0, 0.0, b, b, 0.0, 0.0, 0.0, a]
     call bmc_3_10%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
     call check(S_target,T_target, S,T, msg=' test_boxmc_select_cases_diff_srcsideface src = 10')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir45_east_west_distortion(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:)
    real(ireal_dp), parameter :: dx=1, dy=dx, dz=dx

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    !right side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([6,12,18,24]) = verts([6,12,18,24]) + dz

    phi = 90; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0.5, 0.5, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

    phi = 90; theta = 45
    src = 2
    T_target = [ real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src2_case1')

    phi = -90; theta = 45
    src = 2
    T_target = [ real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src2_case2')


    !right side down
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([6,12,18,24]) = verts([6,12,18,24]) - dz

    phi = -90; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0.5, 0.5, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_down_src1_case1')

    phi = 90; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_down_src1_case2')

    phi = 90; theta = 45
    src = 2
    T_target = [ real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_down_src2_case1')

    phi = -90; theta = 45
    src = 2
    T_target = [ real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_down_src2_case2')
  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir45_north_south_distortion(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:)
    real(ireal_dp), parameter :: dx=1, dy=dx, dz=dx

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    !back side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([9,12,21,24]) = verts([9,12,21,24]) + dz

    phi = 0; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0.5, 0, 0.5]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src1')

    phi = 0; theta = 45
    src = 3
    T_target = [ real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src2_case1')

    phi = 180; theta = 45
    src = 3
    T_target = [ real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src2_case2')


    !back side down
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([9,12,21,24]) = verts([9,12,21,24]) - dz

    phi = 180; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0.5, 0, 0.5]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src1_case1')

    phi = 0; theta = 45
    src = 1
    T_target = [ real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src1_case2')

    phi = 0; theta = 45
    src = 3
    T_target = [ real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src2_case1')

    phi = 180; theta = 45
    src = 3
    T_target = [ real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src2_case2')
  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_dir2dir3_geometric_coeff_corr(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy= dx, dz=1 * dx
    real(irealLUT) :: v(9), v_mc(9)
    real(ireals) :: sundir(3)
    integer(iintegers) :: itheta, iphi

    !bg  = [ 0e-0_ireal_dp/dz, 0e-0_ireal_dp, 1._ireal_dp/2 ]
    !         absorption            scattering        asymmetry parameter
    bg = [ -log(5e-1_ireal_dp)/dz, 0e-0_ireal_dp/dz, 1._ireal_dp/2 ]
    S_target = zero
    iphi=30
    itheta=5
    do iphi=0,360,30
      do itheta=10,50,20
        phi = real(iphi, ireals)
        theta = real(itheta, ireals)

        call setup_default_unit_cube_geometry(dx, dy, dz, verts)
        verts_dtd = verts
        !verts_dtd([3,9,15,21]) = verts_dtd([3,9,15,21]) + dz
        !verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + 2*dz
        verts_dtd([9,12,21,24]) = verts_dtd([9,12,21,24]) + dz / 2
        !verts_dtd([18,6]) = verts_dtd([18,6]) + dz/4

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v(src:3**2:3) = real(T, irealLUT)
        enddo

        sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]

        print *, cstr('regular not corrected', 'red')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        call dir2dir3_geometric_coeff_corr(verts_dtd, sundir, bg(1), v)

        print *, cstr('regular corrected', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
!        @assertEqual(v_mc([2,3,4]), v([2,3,4]), max(maxval(v_mc([2,3,4]))*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
      enddo
    enddo

  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_dir2dir3_geometric_coeffs(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), dimension(24) :: verts_dtd
    real( ireal_dp), parameter :: dx=1, dy= dx, dz=1 * dx
    real(irealLUT) :: v(9), v_mc(9)
    real(ireals) :: sundir(3)
    integer(iintegers) :: itheta, iphi

    bg = [2.7245631510506538E-006_ireal_dp, 0e-0_ireal_dp/dz, 1._ireal_dp/2 ]
    S_target = zero
    iphi=1
    itheta=40
    !do iphi=0,360,30
    !  do itheta=10,50,20
        phi = real(iphi, ireals)
        theta = real(itheta, ireals)

        verts_dtd(1:3)   = [   0.0000000000000000_ireals,        0.0000000000000000_ireals,        0.0000000000000000_ireals]
        verts_dtd(4:6)   = [   100.00000000000000_ireals,        0.0000000000000000_ireals,        0.0000000000000000_ireals]
        verts_dtd(7:9)   = [   0.0000000000000000_ireals,        100.00000000000000_ireals,        2.2658245368802454E-002_ireals]
        verts_dtd(10:12) = [   100.00000000000000_ireals,        100.00000000000000_ireals,        2.2658245368802454E-002_ireals]
        verts_dtd(13:15) = [   0.0000000000000000_ireals,        0.0000000000000000_ireals,        147.47492293885080_ireals]
        verts_dtd(16:18) = [   100.00000000000000_ireals,        0.0000000000000000_ireals,        147.47492293885080_ireals]
        verts_dtd(19:21) = [   0.0000000000000000_ireals,        100.00000000000000_ireals,        147.49758118421960_ireals]
        verts_dtd(22:24) = [   100.00000000000000_ireals,        100.00000000000000_ireals,        147.49758118421960_ireals]


        sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]

        call dir2dir3_geometric_coeff_corr(verts_dtd, sundir, bg(1), v)

        print *, cstr('regular corrected', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
!        @assertEqual(v_mc([2,3,4]), v([2,3,4]), max(maxval(v_mc([2,3,4]))*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
    !  enddo
    !enddo
    end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

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
