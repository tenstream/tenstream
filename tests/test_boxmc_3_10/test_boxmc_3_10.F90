module test_boxmc_3_10
  use m_boxmc, only : t_boxmc, t_boxmc_3_10
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, irealLUT, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_optprop, only : dir2dir3_coeff_corr_zy, dir2dir3_coeff_corr_zx, &
    dir2dir3_coeff_corr_xx, dir2dir3_coeff_corr_xy, &
    dir2dir3_coeff_corr_yy, dir2dir3_coeff_corr_src_x
  use m_helper_functions, only : spherical_2_cartesian, cstr

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

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-2
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

    do iphi=0,360,30
      phi = real(iphi, ireal_dp)
      do src = 2,3
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
  subroutine test_boxmc_distorted_cube_east_west_distortion_gomtc_corr(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real(ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: V(3)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 270; theta = 20
    src = 1

    !right side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    V = real(T, irealLUT)

    verts_dtd = verts
    verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    V = real(T, irealLUT)
    print *, 'regular not corrected', V(2), V(3), V(1)
    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    call dir2dir3_coeff_corr_zx(V, verts, verts_dtd, sundir)
    print *, 'regular corrected', V(2), V(3), V(1)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, 'distorted', T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_north_south_distortion_gomtc_corr(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real(ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: V(3)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 180; theta = 20
    src = 1

    !right side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    V = real(T, irealLUT)

    verts_dtd = verts
    verts_dtd([9,12,21,24]) = verts_dtd([9,12,21,24]) + dz

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    V = real(T, irealLUT)
    print *, 'regular not corrected', V(2), V(3), V(1)
    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    call dir2dir3_coeff_corr_zy(V, verts, verts_dtd, sundir)
    print *, 'regular corrected', V(2), V(3), V(1)

    !verts([9,12,21,24]) = verts([9,12,21,24]) + dz
    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, 'distorted', T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir2dir2_coeff_corr_xx(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: v(9)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 300; theta = 60
    src = 2

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts
    !verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz / 4
    verts_dtd([18,6]) = verts_dtd([18,6]) + dz/4

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    print *, 'vertices'
    print *, 'A', verts_dtd(1), verts_dtd(2), verts_dtd(3)
    print *, 'B', verts_dtd(4), verts_dtd(5), verts_dtd(6)
    print *, 'C', verts_dtd(7), verts_dtd(8), verts_dtd(9)
    print *, 'D', verts_dtd(10), verts_dtd(11), verts_dtd(12)
    print *, 'E', verts_dtd(13), verts_dtd(14), verts_dtd(15)
    print *, 'F', verts_dtd(16), verts_dtd(17), verts_dtd(18)
    print *, 'G', verts_dtd(19), verts_dtd(20), verts_dtd(21)
    print *, 'H', verts_dtd(22), verts_dtd(23), verts_dtd(24)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    v = 0
    v(2) = real(T(1), irealLUT)
    v(5) = real(T(2), irealLUT)
    v(8) = real(T(3), irealLUT)

    print *, 'regular not corrected', v(5), v(8), v(2)

    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    print *, 'sundir', sundir

    call dir2dir3_coeff_corr_xx(v, verts, verts_dtd, sundir)
    print *, cstr('regular corrected', 'red'), v(5), v(8), v(2)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, 'distorted', T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir2dir2_coeff_corr_xy(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: v(9)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 70; theta = 70
    src = 2

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts
    verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz / 4
    !verts_dtd([18,6]) = verts_dtd([18,6]) + dz/4

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    print *, 'vertices'
    print *, 'A', verts_dtd(1), verts_dtd(2), verts_dtd(3)
    print *, 'B', verts_dtd(4), verts_dtd(5), verts_dtd(6)
    print *, 'C', verts_dtd(7), verts_dtd(8), verts_dtd(9)
    print *, 'D', verts_dtd(10), verts_dtd(11), verts_dtd(12)
    print *, 'E', verts_dtd(13), verts_dtd(14), verts_dtd(15)
    print *, 'F', verts_dtd(16), verts_dtd(17), verts_dtd(18)
    print *, 'G', verts_dtd(19), verts_dtd(20), verts_dtd(21)
    print *, 'H', verts_dtd(22), verts_dtd(23), verts_dtd(24)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    v = 0
    v(2) = real(T(1), irealLUT)
    v(5) = real(T(2), irealLUT)
    v(8) = real(T(3), irealLUT)

    print *, 'regular not corrected', v(5), v(8), v(2)

    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    print *, 'sundir', sundir

    call dir2dir3_coeff_corr_xx(v, verts, verts_dtd, sundir)
    print *, cstr('regular xx corrected', 'yellow'), v(5), v(8), v(2)
    call dir2dir3_coeff_corr_xy(v, verts, verts_dtd, sundir)
    print *, cstr('regular xy corrected', 'red'), v(5), v(8), v(2)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, cstr('distorted', 'green'), T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir2dir2_coeff_corr_yy(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: v(9)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 0; theta = 30
    src = 3

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts
    !verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz / 4
    verts_dtd([18,6]) = verts_dtd([18,6]) + dz/4

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    print *, 'vertices'
    print *, 'A', verts_dtd(1), verts_dtd(2), verts_dtd(3)
    print *, 'B', verts_dtd(4), verts_dtd(5), verts_dtd(6)
    print *, 'C', verts_dtd(7), verts_dtd(8), verts_dtd(9)
    print *, 'D', verts_dtd(10), verts_dtd(11), verts_dtd(12)
    print *, 'E', verts_dtd(13), verts_dtd(14), verts_dtd(15)
    print *, 'F', verts_dtd(16), verts_dtd(17), verts_dtd(18)
    print *, 'G', verts_dtd(19), verts_dtd(20), verts_dtd(21)
    print *, 'H', verts_dtd(22), verts_dtd(23), verts_dtd(24)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    v = 0
    v(3) = real(T(1), irealLUT)
    v(6) = real(T(2), irealLUT)
    v(9) = real(T(3), irealLUT)

    print *, 'regular not corrected', v(6), v(9), v(3)

    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    print *, 'sundir', sundir

    call dir2dir3_coeff_corr_yy(v, verts, verts_dtd, sundir)
    print *, 'regular corrected', v(6), v(9), v(3)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, 'distorted', T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir2dir3_coeff_corr_zx(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: v(9)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 300; theta = 20
    src = 1

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts
    !verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz / 4
    verts_dtd([18,6]) = verts_dtd([18,6]) + dz

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    print *, 'vertices'
    print *, 'A', verts_dtd(1), verts_dtd(2), verts_dtd(3)
    print *, 'B', verts_dtd(4), verts_dtd(5), verts_dtd(6)
    print *, 'C', verts_dtd(7), verts_dtd(8), verts_dtd(9)
    print *, 'D', verts_dtd(10), verts_dtd(11), verts_dtd(12)
    print *, 'E', verts_dtd(13), verts_dtd(14), verts_dtd(15)
    print *, 'F', verts_dtd(16), verts_dtd(17), verts_dtd(18)
    print *, 'G', verts_dtd(19), verts_dtd(20), verts_dtd(21)
    print *, 'H', verts_dtd(22), verts_dtd(23), verts_dtd(24)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    v = 0
    v(1) = real(T(1), irealLUT)
    v(4) = real(T(2), irealLUT)
    v(7) = real(T(3), irealLUT)

    print *, 'regular not corrected', v(4), v(7), v(1)

    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    print *, 'sundir', sundir

    call dir2dir3_coeff_corr_zy(v, verts, verts_dtd, sundir)
    print *, 'regular corrected', v(4), v(7), v(1)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, 'distorted', T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_distorted_cube_dir2dir3_coeff_corr_y_trgt_x_src(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_dtd(:)
    real( ireal_dp), parameter :: dx=1, dy=dx, dz=dx
    real(irealLUT) :: v(9)
    real(ireals) :: sundir(3)

    bg  = [0e-0_ireal_dp/dz, 0._ireal_dp, 1._ireal_dp/2 ]
    S_target = zero

    phi = 90; theta = 50
    src = 2

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts
    verts_dtd([3,9,15,21]) = verts_dtd([3,9,15,21]) + dz
!    verts_dtd([6,12,18,24]) = verts_dtd([6,12,18,24]) + dz / 2
    !verts_dtd([18,6]) = verts_dtd([18,6]) + dz/4

    print *, 'vertices'
    print *, 'A', verts(1), verts(2), verts(3)
    print *, 'B', verts(4), verts(5), verts(6)
    print *, 'C', verts(7), verts(8), verts(9)
    print *, 'D', verts(10), verts(11), verts(12)
    print *, 'E', verts(13), verts(14), verts(15)
    print *, 'F', verts(16), verts(17), verts(18)
    print *, 'G', verts(19), verts(20), verts(21)
    print *, 'H', verts(22), verts(23), verts(24)

    print *, 'vertices'
    print *, 'A', verts_dtd(1), verts_dtd(2), verts_dtd(3)
    print *, 'B', verts_dtd(4), verts_dtd(5), verts_dtd(6)
    print *, 'C', verts_dtd(7), verts_dtd(8), verts_dtd(9)
    print *, 'D', verts_dtd(10), verts_dtd(11), verts_dtd(12)
    print *, 'E', verts_dtd(13), verts_dtd(14), verts_dtd(15)
    print *, 'F', verts_dtd(16), verts_dtd(17), verts_dtd(18)
    print *, 'G', verts_dtd(19), verts_dtd(20), verts_dtd(21)
    print *, 'H', verts_dtd(22), verts_dtd(23), verts_dtd(24)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    v = 0
    v(2) = real(T(1), irealLUT)
    v(5) = real(T(2), irealLUT)
    v(8) = real(T(3), irealLUT)

    print *, cstr('regular not corrected', 'red'), v(5), v(8), v(2)

    sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals)) * [-one, -one, one]
    print *, 'sundir', sundir

    call dir2dir3_coeff_corr_src_x(v, verts_dtd, sundir)
    print *, cstr('regular corrected', 'blue'), v(5), v(8), v(2)

    call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,verts_dtd,S,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    print *, cstr('montecarlo distorted', 'green'), T_target(2), T_target(3), T_target(1)

    !call check(S_target,T_target, S,T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

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
