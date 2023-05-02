module test_pprts_symmetry

  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, irealLUT, zero, one, pi, mpiint

#include "petsc/finclude/petsc.h"
  use petsc

  use m_boxmc_geometry, only: setup_default_unit_cube_geometry
  use m_pprts_base, only: t_solver, t_solver_3_10, t_solver_8_10, destroy_pprts
  use m_pprts, only: init_pprts, set_optical_properties, &
                     solve_pprts, set_angles, pprts_get_result_toZero
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: &
    & CHKERR, &
    & colored_str_by_range, &
    & cstr, &
    & spherical_2_cartesian, &
    & toStr

  use m_optprop, only: &
    & t_optprop, &
    & t_optprop_cube, &
    & t_optprop_3_10, &
    & t_optprop_3_16, &
    & t_optprop_3_24, &
    & t_optprop_8_10, &
    & t_optprop_8_16

  use pfunit_mod

  implicit none

  class(t_optprop_cube), allocatable :: OPP
  type(t_solver_3_10) :: solver_3_10
  type(t_solver_8_10) :: solver_8_10

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    if (allocated(OPP)) then
      call OPP%destroy(ierr); call CHKERR(ierr)
      deallocate (OPP)
    end if

    if (solver_8_10%linitialized) then
      call destroy_pprts(solver_8_10, lfinalizepetsc=.true.)
    end if
    if (solver_3_10%linitialized) then
      call destroy_pprts(solver_3_10, lfinalizepetsc=.true.)
    end if
  end subroutine teardown

  @test(npes=[1])
  subroutine test_pprts_symmetry_check_against_bmc(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr
    real(irealLUT), parameter :: eps = 1e-2

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    allocate (t_optprop_3_10 :: OPP)
    call test_east_west_symmetry(OPP)
    call test_north_south_symmetry(OPP)
    deallocate (OPP)

    allocate (t_optprop_8_10 :: OPP)
    call test_east_west_symmetry(OPP)
    call test_north_south_symmetry(OPP)
    deallocate (OPP)

    allocate (t_optprop_3_16 :: OPP)
    call test_east_west_symmetry(OPP)
    call test_north_south_symmetry(OPP)
    deallocate (OPP)

    allocate (t_optprop_8_16 :: OPP)
    call test_east_west_symmetry(OPP)
    call test_north_south_symmetry(OPP)
    deallocate (OPP)

    allocate (t_optprop_3_24 :: OPP)
    call test_east_west_symmetry(OPP)
    call test_north_south_symmetry(OPP)
    deallocate (OPP)

  contains
    subroutine test_east_west_symmetry(OPP)
      class(t_optprop_cube) :: OPP
      real(ireals), parameter :: tauz = 1, w0 = 1, g = 1
      real(ireals) :: vertices(24)
      real(ireals), parameter :: dx = 1, dy = 1, dz = 2
      real(irealLUT) :: angles(2)
      logical, parameter :: ldir = .false.

      real(irealLUT), allocatable :: dir2dir(:), dir2diff(:), tmp(:)
      real(irealLUT), allocatable :: dir2dir_rot(:), dir2diff_rot(:)
      integer(iintegers) :: i, j, d0, d1
      integer(iintegers) :: N

      real(irealLUT), parameter :: color_limits(5) = [0., 1e-5, 1e-3, .1, 1.]
      character(len=*), parameter :: colors(4) = [character(len=10) :: 'black', 'blue', 'green', 'red']

      call init_mpi_data_parameters(comm)

      call read_commandline_options(comm)

      call OPP%init(comm)

      call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

      allocate (dir2dir(OPP%LUT%dir_streams**2))
      allocate (dir2diff(OPP%LUT%diff_streams * OPP%LUT%dir_streams))
      allocate (dir2dir_rot(size(dir2dir)))
      allocate (dir2diff_rot(size(dir2diff)))

      N = OPP%LUT%diff_streams / 3
      allocate (tmp(OPP%LUT%diff_streams))

      !------------------------------
      angles = [90.001, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ') '// &
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do
      print *, ''

      angles = [-90.001, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff_rot, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff_rot(i:size(dir2diff_rot):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ')'//&
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do

      call OPP%dir2diff_coeff_symmetry(dir2diff, .true., .false.)

      print *, cstr(' first angles but rotated for', 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ')'//&
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do

      do i = 1, OPP%LUT%dir_streams
        @mpiassertEqual(dir2diff_rot, dir2diff, eps, 'computed two opposite directions with BMC and then applied east_west symmetry on first one')
      end do

      angles = [89.999, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ')'//&
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do
      print *, ''

      !------------------------------
      angles = [-89.999, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff_rot, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff_rot(i:size(dir2diff_rot):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ')'//&
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do

      call OPP%dir2diff_coeff_symmetry(dir2diff, .true., .false.)

      print *, cstr(' first angles but rotated for', 'blue')
      do i = 1, OPP%LUT%dir_streams
        tmp = dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
        do j = 0, size(tmp) / N
          d0 = 1 + j * N; d1 = min((j + 1) * N, size(tmp))
          print *, 'dir2diff src '//toStr(i)//' : (', d0, '-', d1, ')'//&
            & colored_str_by_range(tmp(d0:d1), color_limits, colors)
        end do
      end do

      do i = 1, OPP%LUT%dir_streams
        @mpiassertEqual(dir2diff_rot, dir2diff, eps, 'computed two opposite directions with BMC and then applied east_west symmetry on first one')
      end do

      call OPP%destroy(ierr); call CHKERR(ierr)

      call PetscFinalize(ierr)
    end subroutine
    subroutine test_north_south_symmetry(OPP)
      class(t_optprop_cube) :: OPP
      real(ireals) :: vertices(24)
      real(ireals), parameter :: tauz = 1, w0 = 1, g = 1
      real(ireals), parameter :: dx = 1, dy = 1, dz = 2
      real(irealLUT) :: angles(2)
      logical, parameter :: ldir = .false.

      real(irealLUT), allocatable :: dir2dir(:), dir2diff(:)
      real(irealLUT), allocatable :: dir2dir_rot(:), dir2diff_rot(:)
      integer(iintegers) :: i

      call init_mpi_data_parameters(comm)

      call read_commandline_options(comm)

      call OPP%init(comm)

      call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

      allocate (dir2dir(OPP%LUT%dir_streams**2))
      allocate (dir2diff(OPP%LUT%diff_streams * OPP%LUT%dir_streams))
      allocate (dir2dir_rot(size(dir2dir)))
      allocate (dir2diff_rot(size(dir2diff)))

      angles = [0.001, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
      end do

      angles = [179.999, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff_rot, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff_rot(i:size(dir2diff_rot):OPP%LUT%dir_streams)
      end do

      call OPP%dir2diff_coeff_symmetry(dir2diff, .false., .true.)

      print *, cstr(' first angles but rotated for', 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
      end do

      do i = 1, OPP%LUT%dir_streams
        @mpiassertEqual(dir2diff_rot, dir2diff, eps, 'computed two opposite directions with BMC and then applied north_south symmetry on first one')
      end do

      !--------------------
      angles = [-0.001, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
      end do

      angles = [180.001, 30.]
      call OPP%get_coeff_bmc(vertices, tauz, w0, g, ldir, dir2diff_rot, angles)

      print *, cstr(' BMC coeff for angles '//toStr(angles), 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff_rot(i:size(dir2diff_rot):OPP%LUT%dir_streams)
      end do

      call OPP%dir2diff_coeff_symmetry(dir2diff, .false., .true.)

      print *, cstr(' first angles but rotated for', 'blue')
      do i = 1, OPP%LUT%dir_streams
        print *, 'dir2diff src', i, dir2diff(i:size(dir2diff):OPP%LUT%dir_streams)
      end do

      do i = 1, OPP%LUT%dir_streams
        @mpiassertEqual(dir2diff_rot, dir2diff, eps, 'computed two opposite directions with BMC and then applied north_south symmetry on first one')
      end do

      call OPP%destroy(ierr); call CHKERR(ierr)

      call PetscFinalize(ierr)
    end subroutine
  end subroutine

  @test(npes=[1])
  subroutine test_pprts_symmetry_roundtrip(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr

    real(irealLUT), allocatable :: dir2dir(:), dir2diff(:)
    integer(iintegers) :: i

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    allocate (t_optprop_3_10 :: OPP)
    call this_test(OPP)
    deallocate (OPP)

    allocate (t_optprop_3_16 :: OPP)
    call this_test(OPP)
    deallocate (OPP)

    allocate (t_optprop_3_24 :: OPP)
    call this_test(OPP)
    deallocate (OPP)

    allocate (t_optprop_8_10 :: OPP)
    call this_test(OPP)
    deallocate (OPP)

    allocate (t_optprop_8_16 :: OPP)
    call this_test(OPP)
    deallocate (OPP)

  contains
    subroutine this_test(OPP)
      class(t_optprop_cube) :: OPP

      call init_mpi_data_parameters(comm)

      call read_commandline_options(comm)

      call OPP%init(comm)

      allocate (dir2dir(OPP%LUT%dir_streams**2))
      allocate (dir2diff(OPP%LUT%diff_streams * OPP%LUT%dir_streams))

      do i = 1, ubound(dir2dir, 1)
        dir2dir(i) = real(i, kind=irealLUT)
      end do

      do i = 1, ubound(dir2diff, 1)
        dir2diff(i) = real(i, kind=irealLUT)
      end do

      call OPP%dir2dir_coeff_symmetry(dir2dir, .true., .true.)
      call OPP%dir2dir_coeff_symmetry(dir2dir, .true., .true.)

      do i = 1, ubound(dir2dir, 1)
        @mpiassertEqual(i, dir2dir(i), 'Coeff dir2dir not equal after switching two time north-south and east-west')
      end do

      call OPP%dir2diff_coeff_symmetry(dir2diff, .true., .true.)
      call OPP%dir2diff_coeff_symmetry(dir2diff, .true., .true.)

      do i = 1, ubound(dir2diff, 1)
        @mpiassertEqual(i, dir2diff(i), 'Coeff dir2diff not equal after switching two time north-south and east-west')
      end do

      deallocate (dir2dir)
      deallocate (dir2diff)
      call OPP%destroy(ierr); call CHKERR(ierr)

      call PetscFinalize(ierr)
    end subroutine
  end subroutine

  @test(npes=[1])
  subroutine test_pprts_symmetry_ex1(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers), parameter :: nxp = 5, nyp = 5, nv = 5
    real(ireals), parameter :: dx = 100, dy = dx
    real(ireals), parameter :: phi0 = 10, theta0 = 60
    real(ireals), parameter :: albedo = 0., dz = dx
    real(ireals), parameter :: incSolar = 1000
    real(ireals), parameter :: atolerance = .1

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call this_test(solver_3_10)

  contains
    subroutine this_test(solver)
      class(t_solver), intent(inout) :: solver
      real(ireals) :: dz1d(nv)

      real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g
      real(ireals), allocatable, dimension(:, :, :) :: fdir0, fdn0, fup0, fdiv0
      real(ireals), allocatable, dimension(:, :, :) :: fdir1, fdn1, fup1, fdiv1
      real(ireals), allocatable, dimension(:, :, :) :: fdir2, fdn2, fup2, fdiv2
      real(ireals), allocatable, dimension(:, :, :) :: fdir3, fdn3, fup3, fdiv3

      integer(iintegers) :: j
      integer(iintegers) :: cx, cy      ! global indices of cloud

      dz1d = dz

      call init_pprts(comm, nv, nxp, nyp, dx, dy, spherical_2_cartesian(phi0, theta0), solver, dz1d)

      allocate (kabs(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
      allocate (ksca(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
      allocate (g(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

      kabs = 1._ireals / nv / dz
      ksca = 1._ireals / nv / dz
      g = zero

      cx = int(real(nxp) / 2) + 1
      cy = int(real(nyp) / 2) + 1

      if (cx .le. (solver%C_one%xe + 1) .and. cx .gt. solver%C_one%xs) then
        if (cy .le. (solver%C_one%ye + 1) .and. cy .gt. solver%C_one%ys) then
          kabs(2, cx - solver%C_one%xs, cy - solver%C_one%ys) = 1 / dz
          ksca(2, cx - solver%C_one%xs, cy - solver%C_one%ys) = 1 / dz
          g(2, cx - solver%C_one%xs, cy - solver%C_one%ys) = .9
        end if
      end if

      call set_optical_properties(solver, albedo, kabs, ksca, g)
      call set_angles(solver, spherical_2_cartesian(10._ireals, theta0))
      call solve_pprts(solver, &
        & lthermal=.false., &
        & lsolar=.true., &
        & edirTOA=incSolar, &
        & opt_solution_uid=10_iintegers)

      call set_angles(solver, spherical_2_cartesian(190._ireals, theta0))
      call solve_pprts(solver, &
        & lthermal=.false., &
        & lsolar=.true., &
        & edirTOA=incSolar, &
        & opt_solution_uid=190_iintegers)

      call pprts_get_result_toZero(solver, fdn0, fup0, fdiv0, fdir0, opt_solution_uid=10_iintegers)
      call pprts_get_result_toZero(solver, fdn1, fup1, fdiv1, fdir1, opt_solution_uid=190_iintegers)

      call set_angles(solver, spherical_2_cartesian(100._ireals, theta0))
      call solve_pprts(solver, &
        & lthermal=.false., &
        & lsolar=.true., &
        & edirTOA=incSolar, &
        & opt_solution_uid=100_iintegers)

      call set_angles(solver, spherical_2_cartesian(280._ireals, theta0))
      call solve_pprts(solver, &
        & lthermal=.false., &
        & lsolar=.true., &
        & edirTOA=incSolar, &
        & opt_solution_uid=280_iintegers)

      call pprts_get_result_toZero(solver, fdn2, fup2, fdiv2, fdir2, opt_solution_uid=100_iintegers)
      call pprts_get_result_toZero(solver, fdn3, fup3, fdiv3, fdir3, opt_solution_uid=280_iintegers)

      if (myid .eq. 0) then
        do j = lbound(fdir0, 3), ubound(fdir0, 3)
          print *, j, 'edir0', fdir0(4, :, j)
        end do
        do j = lbound(fdir0, 3), ubound(fdir0, 3)
          print *, j, 'edir1', fdir1(4, :, j)
        end do
        do j = lbound(fdir0, 3), ubound(fdir0, 3)
          print *, j, 'edn0', fdn0(4, :, j)
        end do
        do j = lbound(fdir0, 3), ubound(fdir0, 3)
          print *, j, 'edn1', fdn1(4, :, j)
        end do
        fdir1(:, :, :) = fdir1(:, nxp:1:-1, nyp:1:-1)
        fdn1(:, :, :) = fdn1(:, nxp:1:-1, nyp:1:-1)
        fup1(:, :, :) = fup1(:, nxp:1:-1, nyp:1:-1)
        fdiv1(:, :, :) = fdiv1(:, nxp:1:-1, nyp:1:-1)
        @assertEqual(fdiv0, fdiv1, atolerance, '10 -> 190: divergence not symmetric for azimuth')
        @assertEqual(fdir0, fdir1, atolerance, '10 -> 190: Edirradiation not symmetric for azimuth ')
        @assertEqual(fdn0, fdn1, atolerance, '10 -> 190: Edn radiation not symmetric for azimuth ')
        @assertEqual(fup0, fup1, atolerance, '10 -> 190: Eup radiation not symmetric for azimuth ')

        fdir3(:, :, :) = fdir3(:, nxp:1:-1, nyp:1:-1)
        fdn3(:, :, :) = fdn3(:, nxp:1:-1, nyp:1:-1)
        fup3(:, :, :) = fup3(:, nxp:1:-1, nyp:1:-1)
        fdiv3(:, :, :) = fdiv3(:, nxp:1:-1, nyp:1:-1)
        @assertEqual(fdiv2, fdiv3, atolerance, '10 -> 190: divergence not symmetric for azimuth')
        @assertEqual(fdir2, fdir3, atolerance, '10 -> 190: Edirradiation not symmetric for azimuth ')
        @assertEqual(fdn2, fdn3, atolerance, '10 -> 190: Edn radiation not symmetric for azimuth ')
        @assertEqual(fup2, fup3, atolerance, '10 -> 190: Eup radiation not symmetric for azimuth ')
      end if
    end subroutine
  end subroutine
end module
