module test_boxmc_3_10
  use m_boxmc, only: t_boxmc, t_boxmc_3_10
  use m_data_parameters, only: &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_helper_functions, only: toStr, colored_str_by_range
  use m_optprop_parameters, only: stddev_atol
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi, theta, dx, dy, dz
  real(ireals) :: S(10), T(3), S_target(10), T_target(3)
  real(ireals) :: S_tol(10), T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_10) :: bmc_3_10

  integer(mpiint) :: myid, mpierr, numnodes, comm
  character(len=120) :: msg

  real(ireal_dp), parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp), parameter :: atol = 1e-3, rtol = 1e-2
contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    call bmc_3_10%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if (myid .eq. 0) print *, numnodes, 'Testing Box-MonteCarlo model with tolerances atol/rtol :: ', atol, rtol

    phi = 0
    theta = 0

    dx = 100
    dy = dx
    dz = 50
    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! from top to bot face
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real(exp(-tau), ireals)

    src = 1
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_direct_srctopface top_to_bot')

  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srctopface_45(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! outwards from face 2
    phi = 0; theta = 45 * 1._ireal_dp

    tau = (bg(1) + bg(2)) * dz * sqrt(2 * 1._ireal_dp)

    S_target = zero

    T_target = zero
    T_target(1) = real(exp(-tau) / 2, ireals)
    T_target(3) = real((1 - exp(-tau)) / (2 * tau), ireals)

    src = 1

    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_direct_srctopface_45')

  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_direct_srcsidefaces(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, iphi
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 1._ireal_dp / 2]

    ! along each of the side faces
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = zero

    T_target = zero
    T_target(1) = real((sinh(tau) - cosh(tau) + 1) / tau, ireals)

    do iphi = 0, 360, 30
      phi = real(iphi, ireal_dp)
      do src = 2, 3
        call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
        call check(S_target, T_target, S, T, msg='test_boxmc_select_cases_direct_srcsidefaces')
      end do
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srctopface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    S_target = [0.0, 0.390156, 0.1404375, 0.1404375, 0.0, 0.0, 0.1404375, 0.1404375, 0.0, 0.0]

    T_target = zero

    src = 2   !top face
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srctopface')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srcbottomface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz

    T_target = zero

    S_target = [0.390156, 0.0, 0.0, 0.0, 0.1404375, 0.1404375, 0.0, 0.0, 0.1404375, 0.1404375]

    src = 1
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcbottomface')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_select_cases_diff_srcsideface(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp) :: tau
    real, parameter :: top = 0.56173, a = 0.104806, b = 0.1424402

    ! direct to diffuse tests
    bg = [1e-3_ireal_dp, 0._ireal_dp, 0._ireal_dp]

    ! outwards from face 2
    phi = 0; theta = 0

    tau = (bg(1) + bg(2)) * dz
    T_target = zero

    src = 3
    S_target = [0.0, top, a, 0.0, 0.0, 0.0, b, b, 0.0, 0.0]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 3')

    src = 4
    S_target = [0.0, top, 0.0, a, 0.0, 0.0, b, b, 0.0, 0.0]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 4')

    src = 5
    S_target = [top, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, b, b]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 5')

    src = 6
    S_target = [top, 0.0, 0.0, 0.0, 0.0, a, 0.0, 0.0, b, b]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 6')

    src = 7
    S_target = [0.0, top, b, b, 0.0, 0.0, a, 0.0, 0.0, 0.0]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 7')

    src = 8
    S_target = [0.0, top, b, b, 0.0, 0.0, 0.0, a, 0.0, 0.0]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 8')

    src = 9
    S_target = [top, 0.0, 0.0, 0.0, b, b, 0.0, 0.0, a, 0.0]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 9')

    src = 10
    S_target = [top, 0.0, 0.0, 0.0, b, b, 0.0, 0.0, 0.0, a]
    call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_select_cases_diff_srcsideface src = 10')

  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_distorted_cube_dir45_east_west_distortion(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:)
    real(ireal_dp), parameter :: dx = 1, dy = dx, dz = dx

    bg = [0e-0_ireal_dp / dz, 0._ireal_dp, 1._ireal_dp / 2]
    S_target = zero

    !right side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([6, 12, 18, 24]) = verts([6, 12, 18, 24]) + dz

    phi = 90; theta = 45
    src = 1
    T_target = [real(ireals) :: 0.5, 0.5, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_up_src1')

    phi = 90; theta = 45
    src = 2
    T_target = [real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_up_src2_case1')

    phi = -90; theta = 45
    src = 2
    T_target = [real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_up_src2_case2')

    !right side down
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([6, 12, 18, 24]) = verts([6, 12, 18, 24]) - dz

    phi = -90; theta = 45
    src = 1
    T_target = [real(ireals) :: 0.5, 0.5, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_down_src1_case1')

    phi = 90; theta = 45
    src = 1
    T_target = [real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_down_src1_case2')

    phi = 90; theta = 45
    src = 2
    T_target = [real(ireals) :: 0, 1, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_down_src2_case1')

    phi = -90; theta = 45
    src = 2
    T_target = [real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_cube_dir45_down_src2_case2')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_distorted_cube_dir45_north_south_distortion(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:)
    real(ireal_dp), parameter :: dx = 1, dy = dx, dz = dx

    bg = [0e-0_ireal_dp / dz, 0._ireal_dp, 1._ireal_dp / 2]
    S_target = zero

    !back side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([9, 12, 21, 24]) = verts([9, 12, 21, 24]) + dz

    phi = 0; theta = 45
    src = 1
    T_target = [real(ireals) :: 0.5, 0, 0.5]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src1')

    phi = 0; theta = 45
    src = 3
    T_target = [real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src2_case1')

    phi = 180; theta = 45
    src = 3
    T_target = [real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_up_src2_case2')

    !back side down
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts([9, 12, 21, 24]) = verts([9, 12, 21, 24]) - dz

    phi = 180; theta = 45
    src = 1
    T_target = [real(ireals) :: 0.5, 0, 0.5]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src1_case1')

    phi = 0; theta = 45
    src = 1
    T_target = [real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src1_case2')

    phi = 0; theta = 45
    src = 3
    T_target = [real(ireals) :: 0, 0, 1]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src2_case1')

    phi = 180; theta = 45
    src = 3
    T_target = [real(ireals) :: 1, 0, 0]
    call bmc_3_10%get_coeff(comm, bg, src, .true., phi, theta, verts, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
    call check(S_target, T_target, S, T, msg=' test_boxmc_distorted_ns_cube_dir45_down_src2_case2')
  end subroutine

  @test(npes=[1])
  subroutine test_boxmc_distorted_cube_diff2diff(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireal_dp), allocatable :: verts(:), verts_symm(:)
    real(ireal_dp), parameter :: dx = 1, dy = dx, dz = dx
    real(ireals) :: Ss_trgt(size(S_target), size(S_target))

    bg = [1e-1_ireal_dp / dz, 0._ireal_dp, 1._ireal_dp / 2]
    T_target = zero

    !back side up
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    call setup_default_unit_cube_geometry(dx, dy, dz, verts_symm)

    verts_symm([3, 6, 15, 18]) = verts_symm([3, 6, 15, 18]) - dz
    verts([9, 12, 21, 24]) = verts([9, 12, 21, 24]) + dz
    do src = 1, 8
      print *, 'verts', src, ':', verts(src:size(verts):8)
    end do

    Ss_trgt(1, :) = [real(ireals) :: &
                     2.19971e-01, 0.00000e+00, 2.04549e-02, 2.04649e-02, 1.43955e-01, &
                     1.44041e-01, 9.85783e-02, 0.00000e+00, 2.61983e-01, 3.47848e-02 &
                     ]
    Ss_trgt(2, :) = [real(ireals) :: &
                     0.00000e+00, 2.19969e-01, 1.44047e-01, 1.43993e-01, 2.04579e-02, &
                     2.04490e-02, 3.47453e-02, 2.61988e-01, 0.00000e+00, 9.85785e-02 &
                     ]
    Ss_trgt(3, :) = [real(ireals) :: &
                     5.78161e-02, 4.07155e-01, 1.63150e-01, 0.00000e+00, 0.00000e+00, &
                     0.00000e+00, 2.01940e-01, 1.10228e-01, 0.00000e+00, 0.00000e+00 &
                     ]
    Ss_trgt(4, :) = [real(ireals) :: &
                     5.78408e-02, 4.07360e-01, 0.00000e+00, 1.63117e-01, 0.00000e+00, &
                     0.00000e+00, 2.01804e-01, 1.10175e-01, 0.00000e+00, 0.00000e+00 &
                     ]
    Ss_trgt(5, :) = [real(ireals) :: &
                     4.07360e-01, 5.78282e-02, 0.00000e+00, 0.00000e+00, 1.63196e-01, &
                     0.00000e+00, 0.00000e+00, 0.00000e+00, 1.10103e-01, 2.01799e-01 &
                     ]
    Ss_trgt(6, :) = [real(ireals) :: &
                     4.07411e-01, 5.78942e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
                     1.62993e-01, 0.00000e+00, 0.00000e+00, 1.10239e-01, 2.01763e-01 &
                     ]
    Ss_trgt(7, :) = [real(ireals) :: &
                     2.78895e-01, 9.83542e-02, 2.01884e-01, 2.01884e-01, 0.00000e+00, &
                     0.00000e+00, 1.50649e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00 &
                     ]
    Ss_trgt(8, :) = [real(ireals) :: &
                     0.00000e+00, 7.41315e-01, 1.10193e-01, 1.10301e-01, 0.00000e+00, &
                     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00 &
                     ]
    Ss_trgt(9, :) = [real(ireals) :: &
                     7.41330e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.10215e-01, &
                     1.10246e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00 &
                     ]
    Ss_trgt(10, :) = [real(ireals) :: &
                      9.83269e-02, 2.78907e-01, 0.00000e+00, 0.00000e+00, 2.01946e-01, &
                      2.01699e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.50771e-01 &
                      ]

    do src = 1, size(S_target)
      call run_tests(Ss_trgt(src, :), verts, verts_symm)
    end do

    call setup_default_unit_cube_geometry(dx, dy, dz, verts_symm)
    verts_symm([3, 6, 15, 18]) = verts_symm([3, 6, 15, 18]) - dz
    verts([9, 12, 21, 24]) = verts([9, 12, 21, 24]) + dz

  contains
    subroutine run_tests(S_trgt, vs, vs_symm)
      real(ireal_dp), intent(in) :: vs(:), vs_symm(:)
      real(ireals), intent(in) :: S_trgt(:)

      call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vs, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_trgt, T_target, S, T, msg='test_boxmc_distorted_diff2diff_src'//toStr(src))
      call bmc_3_10%get_coeff(comm, bg, src, .false., phi, theta, vs_symm, S, T, S_tol, T_tol, inp_atol=atol, inp_rtol=rtol)
      call check(S_trgt, T_target, S, T, msg='test_boxmc_distorted_diff2diff_src'//toStr(src)//'_backside_symm')
    end subroutine
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(S_target, T_target, S, T, msg)
    real(ireals), intent(in), dimension(:) :: S_target, T_target, S, T

    character(len=*), optional :: msg
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals), parameter :: test_atol = real(atol, ireals) * real(sigma, ireals)
    real(ireals), parameter :: color_limits(5) = [0., 1., 10., 50., 100.]
    character(len=*), parameter :: colors(4) = [character(len=10) :: 'black', 'blue', 'green', 'red']

    if (myid .eq. 0) then
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
      write (*, FMT='( " diffuse  :: ",10(es12.5) )') S
      write (*, FMT='( " target   :: ",10(es12.5) )') S_target
      write (*, FMT='( " abs diff :: ",10(es12.5) )') S_target - S
      write (*, FMT='( " rel diff[%] ",A )') &
        & colored_str_by_range(100 * abs(S_target - S) / max(sqrt(tiny(S_target)), S_target), color_limits, colors)
      print *, ''
      write (*, FMT='( " direct   :: ", 8(es12.5) )') T
      write (*, FMT='( " target   :: ", 8(es12.5) )') T_target
      write (*, FMT='( " abs diff :: ", 8(es12.5) )') T_target - T
      write (*, FMT='( " rel diff[%] ",A )') &
        & colored_str_by_range(100 * abs(T_target - T) / max(sqrt(tiny(T_target)), T_target), color_limits, colors)
      print *, '---------------------'
      print *, ''

      @assertEqual(S_target, S, test_atol, local_msgS)
      @assertLessThanOrEqual(zero, S)
      @assertGreaterThanOrEqual(one, S)

      @assertEqual(T_target, T, test_atol, local_msgT)
      @assertLessThanOrEqual(zero, T)
      @assertGreaterThanOrEqual(one, T)
    end if
  end subroutine

end module
