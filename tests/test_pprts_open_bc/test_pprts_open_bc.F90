module test_pprts_open_bc
  ! Pins the behaviour of the open boundary conditions (-pprts_open_bc) for direct radiation.
  !
  ! The horizontal DMDA is always periodic, open boundaries are emulated by
  !   * not sending direct radiation across the outer domain edges and
  !   * feeding the sunward edges with the side fluxes of a single column solve.
  ! Diffuse radiation is not covered here, we therefore use purely absorbing media.
  !
  ! Options:
  !   -pprts_open_bc     open boundaries in x and y
  !   -pprts_open_bc_x   open boundaries at the west/east edges only
  !   -pprts_open_bc_y   open boundaries at the south/north edges only
  !   -pprts_open_bc_2d  inflow of an edge cell is its own outflow (zero gradient across the edge), this is the default.
  !                      If set to false, the inflow is the side flux of a single column solve

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, imp_ireals, &
    & zero, one

  use mpi, only: MPI_SUM, MPI_INTEGER

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  use pfunit_mod

  use m_helper_functions, only: &
    & CHKERR, &
    & deg2rad, &
    & ind_1d_to_nd, &
    & insert_petsc_opt, &
    & spherical_2_cartesian, &
    & toStr

  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result, pprts_get_result_toZero

  use m_buildings, only: &
    & check_buildings_consistency, &
    & destroy_buildings, &
    & faceidx_by_cell_plus_offset, &
    & init_buildings, &
    & t_pprts_buildings

  implicit none

  integer(iintegers), parameter :: Nx = 8, Ny = 8, Nz = 8
  real(ireals), parameter :: dx = 100, dy = dx, dz = 100
  real(ireals), parameter :: theta0 = 60
  real(ireals), parameter :: incSolar = 1000
  real(ireals), parameter :: albedo = 0
  real(ireals), parameter :: kabs_clearsky = .2_ireals / (dz * Nz)
  real(ireals), parameter :: kabs_cloud = 5._ireals / dz
  integer(iintegers), parameter :: kcld = 2

  ! crater: terrain made of opaque building cells that rises towards the domain edges.
  ! The rim height differs on all four edges, i.e. the terrain height is not periodic.
  integer(iintegers), parameter :: Ncrater = 12 ! horizontal extent of the crater domain
  integer(iintegers), parameter :: rim_west = 4, rim_east = 2, rim_south = 3, rim_north = 1
  ! surroundings to embed the crater in a periodic domain,
  ! wide enough that no shadow of the continued terrain wraps around into the crater
  integer(iintegers), parameter :: Npad = 16
  ! and an even wider one to check that the result does not depend on the width of the surroundings
  integer(iintegers), parameter :: Npad_wide = 24

  ! sun azimuths, four that are grid aligned and all four diagonal quadrants
  real(ireals), parameter :: phis(8) = [real(ireals) :: 0, 90, 180, 270, 45, 135, 225, 315]

  ! set by the solve helpers if a solver did not pick up the requested -pprts_open_bc options.
  ! We must not assert inside of the helpers, see the note on the structure of the tests below
  logical :: loption_mismatch

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
    loption_mismatch = .false.
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call set_bc_options(.false., .false., .false.)
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  ! always set all of the options, they stay in the options database
  subroutine set_bc_options(lopen_x, lopen_y, l2d)
    logical, intent(in) :: lopen_x, lopen_y, l2d
    integer(mpiint) :: ierr
    call insert_petsc_opt('-pprts_open_bc '//merge('yes', 'no ', lopen_x .and. lopen_y), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_x '//merge('yes', 'no ', lopen_x), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_y '//merge('yes', 'no ', lopen_y), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_2d '//merge('yes', 'no ', l2d), ierr); call CHKERR(ierr)
  end subroutine

  subroutine check_bc_options(solver, lopen_x, lopen_y, l2d)
    class(t_solver), intent(in) :: solver
    logical, intent(in) :: lopen_x, lopen_y, l2d
    if (lopen_x .neqv. solver%lopen_bc_x) loption_mismatch = .true.
    if (lopen_y .neqv. solver%lopen_bc_y) loption_mismatch = .true.
    if ((lopen_x .or. lopen_y) .neqv. solver%lopen_bc) loption_mismatch = .true.
    if (l2d .neqv. solver%lopen_bc_2d) loption_mismatch = .true.
  end subroutine

  ! Beer-Lambert direct irradiance through a horizontally homogeneous clearsky atmosphere
  pure function edir_clearsky() result(edir)
    real(ireals) :: edir(Nz + 1)
    real(ireals) :: mu
    integer(iintegers) :: k
    mu = cos(deg2rad(theta0))
    do k = 1, Nz + 1
      edir(k) = incSolar * mu * exp(-kabs_clearsky * dz * real(k - 1, ireals) / mu)
    end do
  end function

  ! global index of the cloud column, i.e. the outermost column at the downwind domain edge(s)
  ! for grid aligned sun angles the cloud is put in the middle of the edge
  subroutine downwind_edge_column(phi0, icld, jcld)
    real(ireals), intent(in) :: phi0
    integer(iintegers), intent(out) :: icld, jcld
    real(ireals) :: sundir(3)
    real(ireals), parameter :: eps = 1e-3_ireals

    sundir = spherical_2_cartesian(phi0, theta0) ! points away from the sun

    if (sundir(1) .gt. eps) then
      icld = Nx
    else if (sundir(1) .lt. -eps) then
      icld = 1
    else
      icld = Nx / 2
    end if

    if (sundir(2) .gt. eps) then
      jcld = Ny
    else if (sundir(2) .lt. -eps) then
      jcld = 1
    else
      jcld = Ny / 2
    end if
  end subroutine

  ! solve for direct radiation, edir has global shape (Nz+1, Nx, Ny) and is only allocated on rank 0
  ! lopen_x, lopen_y default to lopen_bc
  subroutine solve_edir(comm, phi0, lopen_bc, lcloud, edir, abso, lopen_x, lopen_y, l2d)
    integer(mpiint), intent(in) :: comm
    real(ireals), intent(in) :: phi0
    logical, intent(in) :: lopen_bc, lcloud
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: edir
    real(ireals), allocatable, dimension(:, :, :), intent(out), optional :: abso ! dim (Nz, Nx, Ny), only on rank 0
    logical, intent(in), optional :: lopen_x, lopen_y, l2d

    class(t_solver), allocatable :: solver
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g
    real(ireals), allocatable, dimension(:, :, :) :: edn, eup, labso
    real(ireals) :: dz1d(Nz)
    integer(iintegers) :: icld, jcld, i, j
    logical :: lx, ly, l2
    integer(mpiint) :: ierr

    lx = lopen_bc
    ly = lopen_bc
    l2 = .false.
    if (present(lopen_x)) lx = lopen_x
    if (present(lopen_y)) ly = lopen_y
    if (present(l2d)) l2 = l2d
    call set_bc_options(lx, ly, l2)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    dz1d = dz
    call init_pprts(comm, Nz, Nx, Ny, dx, dy, spherical_2_cartesian(phi0, theta0), solver, dz1d=dz1d)
    call check_bc_options(solver, lx, ly, l2)

    associate (C => solver%C_one)
      allocate (kabs(C%zm, C%xm, C%ym), source=kabs_clearsky)
      allocate (ksca(C%zm, C%xm, C%ym), source=zero)
      allocate (g(C%zm, C%xm, C%ym), source=zero)

      if (lcloud) then
        call downwind_edge_column(phi0, icld, jcld)
        do j = C%ys, C%ye
          do i = C%xs, C%xe
            if (i + 1 .eq. icld .and. j + 1 .eq. jcld) then
              kabs(kcld, i - C%xs + 1, j - C%ys + 1) = kabs_cloud
            end if
          end do
        end do
      end if
    end associate

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=incSolar)
    call pprts_get_result_toZero(solver, edn, eup, labso, edir)
    if (present(abso)) then
      if (allocated(labso)) allocate (abso, source=labso)
    end if
    call destroy_pprts(solver, lfinalizepetsc=.false.)
  end subroutine

  ! terrain height in number of cells above ground, i,j are global indices starting at 1.
  ! Outside of the crater domain, the terrain continues with the height of the nearest edge column
  pure function crater_height(i, j) result(h)
    integer(iintegers), intent(in) :: i, j
    integer(iintegers) :: h, ic, jc
    ic = min(max(i, 1_iintegers), Ncrater)
    jc = min(max(j, 1_iintegers), Ncrater)
    h = 0
    h = max(h, rim_west - (ic - 1))
    h = max(h, rim_east - (Ncrater - ic))
    h = max(h, rim_south - (jc - 1))
    h = max(h, rim_north - (Ncrater - jc))
  end function

  ! solve for direct radiation in the crater which is surrounded by <pad> columns of terrain that continues the rim outwards.
  ! Results are cut to the crater domain, i.e. without the padding
  !   edir  (Nz+1, Ncrater, Ncrater) and abso (Nz, Ncrater, Ncrater) only on rank 0
  !   bedir (6, Nz, Ncrater, Ncrater) direct irradiance on the building faces, available on all ranks, -1 where there is no face
  !   bcount (6, Nz, Ncrater, Ncrater) number of ranks that provided a value for a building face, available on all ranks
  subroutine solve_crater(comm, phi0, lopen_bc, pad, edir, abso, bedir, bcount, l2d)
    integer(mpiint), intent(in) :: comm
    real(ireals), intent(in) :: phi0
    logical, intent(in) :: lopen_bc
    integer(iintegers), intent(in) :: pad
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: edir, abso
    real(ireals), intent(out) :: bedir(:, :, :, :)
    integer(mpiint), intent(out), optional :: bcount(:, :, :, :)
    logical, intent(in), optional :: l2d

    class(t_solver), allocatable :: solver
    type(t_pprts_buildings), allocatable :: buildings
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g
    real(ireals), allocatable, dimension(:, :, :) :: edn, eup, gabso, gedir
    real(ireals), allocatable, dimension(:, :, :) :: ledn, leup, labso, ledir
    real(ireals), allocatable :: lbedir(:, :, :, :)
    integer(mpiint), allocatable :: lbcount(:, :, :, :)
    real(ireals) :: dz1d(Nz)
    integer(iintegers) :: Nglob, Nfaces, m, i, j, k, h, faceid, idx(4)
    logical :: l2
    integer(mpiint) :: ierr

    l2 = .false.
    if (present(l2d)) l2 = l2d
    call set_bc_options(lopen_bc, lopen_bc, l2)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    Nglob = Ncrater + 2 * pad
    dz1d = dz
    call init_pprts(comm, Nz, Nglob, Nglob, dx, dy, spherical_2_cartesian(phi0, theta0), solver, dz1d=dz1d)
    call check_bc_options(solver, lopen_bc, lopen_bc, l2)

    associate (C => solver%C_one)
      allocate (kabs(C%zm, C%xm, C%ym), source=kabs_clearsky)
      allocate (ksca(C%zm, C%xm, C%ym), source=zero)
      allocate (g(C%zm, C%xm, C%ym), source=zero)

      Nfaces = 0
      do j = C%ys, C%ye
        do i = C%xs, C%xe
          Nfaces = Nfaces + 6 * crater_height(i + 1 - pad, j + 1 - pad)
        end do
      end do

      call init_buildings(buildings, [integer(iintegers) :: 6, C%zm, C%xm, C%ym], Nfaces, ierr); call CHKERR(ierr)

      m = 0
      do j = C%ys, C%ye
        do i = C%xs, C%xe
          h = crater_height(i + 1 - pad, j + 1 - pad)
          do k = Nz - h + 1, Nz
            do faceid = 1, 6
              m = m + 1
              buildings%iface(m) = faceidx_by_cell_plus_offset(buildings%da_offsets, k, i - C%xs + 1, j - C%ys + 1, faceid)
              buildings%albedo(m) = albedo
            end do
          end do
        end do
      end do
      call check_buildings_consistency(buildings, C%zm, C%xm, C%ym, ierr); call CHKERR(ierr)
    end associate

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=incSolar, opt_buildings=buildings)

    ! direct radiation on the building faces
    call pprts_get_result(solver, ledn, leup, labso, ledir, opt_buildings=buildings)

    allocate (lbedir(6, Nz, Ncrater, Ncrater), source=zero)
    allocate (lbcount(6, Nz, Ncrater, Ncrater), source=0_mpiint)
    associate (C => solver%C_one)
      do m = 1, size(buildings%iface)
        call ind_1d_to_nd(buildings%da_offsets, buildings%iface(m), idx)
        i = idx(3) + C%xs - pad
        j = idx(4) + C%ys - pad
        if (i .lt. 1 .or. i .gt. Ncrater .or. j .lt. 1 .or. j .gt. Ncrater) cycle ! terrain in the padding
        lbedir(idx(1), idx(2), i, j) = buildings%edir(m) + one ! shift by one so that we end up with -1 where there is no face
        lbcount(idx(1), idx(2), i, j) = lbcount(idx(1), idx(2), i, j) + 1_mpiint
      end do
    end associate
    call mpi_allreduce(lbedir, bedir, size(lbedir, kind=mpiint), imp_ireals, MPI_SUM, comm, ierr); call CHKERR(ierr)
    bedir = bedir - one
    if (present(bcount)) then
      call mpi_allreduce(lbcount, bcount, size(lbcount, kind=mpiint), MPI_INTEGER, MPI_SUM, comm, ierr); call CHKERR(ierr)
    end if

    call pprts_get_result_toZero(solver, edn, eup, gabso, gedir)
    if (allocated(gedir)) then
      allocate (edir, source=gedir(:, pad + 1:pad + Ncrater, pad + 1:pad + Ncrater))
      allocate (abso, source=gabso(:, pad + 1:pad + Ncrater, pad + 1:pad + Ncrater))
    end if

    call destroy_buildings(buildings, ierr); call CHKERR(ierr)
    call destroy_pprts(solver, lfinalizepetsc=.false.)
  end subroutine

  ! Note on the structure of the tests:
  !   results are gathered on rank 0 and can only be checked there.
  !   A failing pFUnit assert returns from the test, so we must not assert in between the (collective) solves
  !   or the other ranks deadlock. Hence, first do all solves, then check.
  !   For the same reason, the helpers only record in loption_mismatch if the solver did not pick up the open bc options.

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_clearsky_is_homogeneous(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: eps_1d = incSolar * 1e-2_ireals ! accuracy of the LUT vs Beer-Lambert
    real(ireals), parameter :: eps_periodic = incSolar * 1e-4_ireals

    real(ireals), allocatable, dimension(:, :, :) :: edir
    real(ireals), dimension(Nz + 1, Nx, Ny, size(phis)) :: edir_open, edir_open_2d, edir_periodic
    real(ireals) :: edir_1d(Nz + 1)
    integer(iintegers) :: iphi, i, j
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    edir_1d = edir_clearsky()
    edir_open = -one
    edir_periodic = -one

    edir_open_2d = -one
    do iphi = 1, size(phis)
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.false., edir=edir)
      if (allocated(edir)) edir_open(:, :, :, iphi) = edir
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.false., edir=edir, l2d=.true.)
      if (allocated(edir)) edir_open_2d(:, :, :, iphi) = edir
      call solve_edir(comm, phis(iphi), lopen_bc=.false., lcloud=.false., edir=edir)
      if (allocated(edir)) edir_periodic(:, :, :, iphi) = edir
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do iphi = 1, size(phis)
      print *, 'phi0', phis(iphi), 'edir 1D      ', edir_1d
      print *, 'phi0', phis(iphi), 'edir open min', minval(minval(edir_open(:, :, :, iphi), dim=3), dim=2)
      print *, 'phi0', phis(iphi), 'edir open max', maxval(maxval(edir_open(:, :, :, iphi), dim=3), dim=2)
    end do

    do iphi = 1, size(phis)
      do j = 1, Ny
        do i = 1, Nx
@assertEqual(edir_1d, edir_open(:, i, j, iphi), eps_1d, 'open bc clearsky edir should match Beer-Lambert, phi0 '//toStr(phis(iphi))//' column '//toStr(i)//','//toStr(j))
        end do
      end do
@assertEqual(edir_periodic(:, :, :, iphi), edir_open(:, :, :, iphi), eps_periodic, 'open bc clearsky edir should be the same as with periodic boundaries, phi0 '//toStr(phis(iphi)))
@assertEqual(edir_periodic(:, :, :, iphi), edir_open_2d(:, :, :, iphi), eps_periodic, '2D open bc clearsky edir should be the same as with periodic boundaries, phi0 '//toStr(phis(iphi)))
    end do
  end subroutine

  ! The absorption of the cells at the inflow edges depends on the side fluxes that enter through the open boundary
  @test(npes=[4, 2, 1])
  subroutine test_open_bc_clearsky_absorption_is_homogeneous(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), allocatable, dimension(:, :, :) :: edir, abso
    real(ireals), dimension(Nz, Nx, Ny, size(phis)) :: abso_open, abso_open_2d, abso_periodic
    real(ireals) :: eps
    integer(iintegers) :: iphi, j
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    abso_open = -one
    abso_open_2d = -one
    abso_periodic = -one

    do iphi = 1, size(phis)
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.false., edir=edir, abso=abso)
      if (allocated(abso)) abso_open(:, :, :, iphi) = abso
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.false., edir=edir, abso=abso, l2d=.true.)
      if (allocated(abso)) abso_open_2d(:, :, :, iphi) = abso
      call solve_edir(comm, phis(iphi), lopen_bc=.false., lcloud=.false., edir=edir, abso=abso)
      if (allocated(abso)) abso_periodic(:, :, :, iphi) = abso
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do iphi = 1, size(phis)
      print *, 'phi0', phis(iphi)
      do j = 1, Ny
        print *, 'abso [mW/m3] at k=4 open/periodic j', j, ':', int(abso_open(4, :, j, iphi) * 1e3_ireals), &
          & ':', int(abso_periodic(4, :, j, iphi) * 1e3_ireals)
      end do
    end do

    do iphi = 1, size(phis)
      eps = maxval(abso_periodic(:, :, :, iphi)) * 1e-2_ireals
@assertEqual(abso_periodic(:, :, :, iphi), abso_open(:, :, :, iphi), eps, 'open bc clearsky absorption should be the same as with periodic boundaries, phi0 '//toStr(phis(iphi)))
@assertEqual(abso_periodic(:, :, :, iphi), abso_open_2d(:, :, :, iphi), eps, '2D open bc clearsky absorption should be the same as with periodic boundaries, phi0 '//toStr(phis(iphi)))
    end do
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_cloud_shadow_does_not_wrap_around(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: eps = incSolar * 1e-4_ireals

    real(ireals), allocatable, dimension(:, :, :) :: edir, abso
    real(ireals), dimension(Nz + 1, Nx, Ny, size(phis)) :: edir_clear, edir_open, edir_open_2d, edir_periodic
    real(ireals), dimension(Nz, Nx, Ny, size(phis)) :: abso_clear, abso_open
    real(ireals) :: max_shadow_periodic, eps_abso
    integer(iintegers) :: iphi, i, j, icld, jcld
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    edir_clear = -one
    edir_open = -one
    edir_open_2d = -one
    edir_periodic = -one
    abso_clear = -one
    abso_open = -one

    do iphi = 1, size(phis)
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.false., edir=edir, abso=abso)
      if (allocated(edir)) edir_clear(:, :, :, iphi) = edir
      if (allocated(abso)) abso_clear(:, :, :, iphi) = abso
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.true., edir=edir, abso=abso)
      if (allocated(edir)) edir_open(:, :, :, iphi) = edir
      if (allocated(abso)) abso_open(:, :, :, iphi) = abso
      call solve_edir(comm, phis(iphi), lopen_bc=.true., lcloud=.true., edir=edir, l2d=.true.)
      if (allocated(edir)) edir_open_2d(:, :, :, iphi) = edir
      call solve_edir(comm, phis(iphi), lopen_bc=.false., lcloud=.true., edir=edir)
      if (allocated(edir)) edir_periodic(:, :, :, iphi) = edir
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do iphi = 1, size(phis)
      call downwind_edge_column(phis(iphi), icld, jcld)
      print *, 'phi0', phis(iphi), 'cloud column', icld, jcld
      do j = 1, Ny
        print *, 'srfc edir open/periodic j', j, ':', int(edir_open(Nz + 1, :, j, iphi)), ':', int(edir_periodic(Nz + 1, &
                                                                                                                 :, j, iphi))
      end do
    end do

    do iphi = 1, size(phis)
      call downwind_edge_column(phis(iphi), icld, jcld)

      ! the cloud itself has to cast a shadow, otherwise this test is void
@assertTrue(edir_open(kcld + 1, icld, jcld, iphi) .lt. edir_clear(kcld + 1, icld, jcld, iphi) * .5_ireals, 'expected a shadow below the cloud, phi0 '//toStr(phis(iphi)))
@assertTrue(edir_open_2d(kcld + 1, icld, jcld, iphi) .lt. edir_clear(kcld + 1, icld, jcld, iphi) * .5_ireals, '2D open bc: expected a shadow below the cloud, phi0 '//toStr(phis(iphi)))

      max_shadow_periodic = zero
      eps_abso = maxval(abso_clear(:, :, :, iphi)) * 1e-3_ireals
      do j = 1, Ny
        do i = 1, Nx
          if (i .eq. icld .and. j .eq. jcld) cycle
          max_shadow_periodic = max(max_shadow_periodic, maxval(edir_clear(:, i, j, iphi) - edir_periodic(:, i, j, iphi)))
@assertEqual(edir_clear(:, i, j, iphi), edir_open(:, i, j, iphi), eps, 'shadow of the cloud at the downwind edge leaked into the domain, phi0 '//toStr(phis(iphi))//' column '//toStr(i)//','//toStr(j))
@assertEqual(abso_clear(:, i, j, iphi), abso_open(:, i, j, iphi), eps_abso, 'shadow of the cloud at the downwind edge leaked into the absorption, phi0 '//toStr(phis(iphi))//' column '//toStr(i)//','//toStr(j))
@assertEqual(edir_clear(:, i, j, iphi), edir_open_2d(:, i, j, iphi), eps, '2D open bc: shadow of the cloud at the downwind edge leaked into the domain, phi0 '//toStr(phis(iphi))//' column '//toStr(i)//','//toStr(j))
        end do
      end do

      ! and make sure that the setup is sensitive to the boundary condition, i.e. the shadow wraps around with cyclic boundaries
@assertTrue(max_shadow_periodic .gt. incSolar * 1e-2_ireals, 'expected the shadow to wrap around with periodic boundaries, phi0 '//toStr(phis(iphi)))
    end do
  end subroutine

  ! -pprts_open_bc_x and -pprts_open_bc_y open the boundaries in one direction only, the other one stays periodic.
  ! With a cloud in the downwind corner of the domain, per direction either the shadow wraps around or it does not
  @test(npes=[4, 2, 1])
  subroutine test_open_bc_x_and_y_can_be_set_individually(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: eps = incSolar * 1e-4_ireals
    real(ireals), parameter :: phis_diag(4) = [real(ireals) :: 45, 135, 225, 315]
    integer(iintegers), parameter :: Nmodes = 6 ! periodic, open, x, y, x with 2D inflow, y with 2D inflow
    logical, parameter :: lxs(Nmodes) = [.false., .true., .true., .false., .true., .false.]
    logical, parameter :: lys(Nmodes) = [.false., .true., .false., .true., .false., .true.]
    logical, parameter :: l2ds(Nmodes) = [.false., .false., .false., .false., .true., .true.]

    real(ireals), allocatable, dimension(:, :, :) :: edir
    real(ireals), dimension(Nz + 1, Nx, Ny, size(phis_diag), Nmodes) :: e
    real(ireals) :: shadow_x(Nmodes), shadow_y(Nmodes), d
    integer(iintegers) :: iphi, imode, icld, jcld, i, j
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    e = -one
    do iphi = 1, size(phis_diag)
      do imode = 1, Nmodes
        call solve_edir(comm, phis_diag(iphi), lopen_bc=.false., lcloud=.true., edir=edir, &
          & lopen_x=lxs(imode), lopen_y=lys(imode), l2d=l2ds(imode))
        if (allocated(edir)) e(:, :, :, iphi, imode) = edir
      end do
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    @assertTrue(all(e .ge. zero), 'missing or negative edir')

    do iphi = 1, size(phis_diag)
      call downwind_edge_column(phis_diag(iphi), icld, jcld)

      ! The cloud sits in the downwind corner icld,jcld. The part of its shadow that leaves through
      !   the x edge comes back at the opposite x edge and travels on through the row jcld, at first
      !   the y edge comes back at the opposite y edge and travels on through the column icld, at first
      ! If x is open, the shadow can only come back through y. From there on it stays in the column icld
      ! because whatever leaves this column in x direction, leaves the domain. Same for y in the row jcld.
      shadow_x = zero
      shadow_y = zero
      do j = 1, Ny
        do i = 1, Nx
          do imode = 3, Nmodes
            d = maxval(abs(e(:, i, j, iphi, 2) - e(:, i, j, iphi, imode))) ! difference to open boundaries
            if (lxs(imode)) then
              if (i .ne. icld) then
@assertEqual(e(:, i, j, iphi, 2), e(:, i, j, iphi, imode), eps, 'open x: shadow wrapped around through the x edge, phi0 '//toStr(phis_diag(iphi))//' mode '//toStr(imode)//' column '//toStr(i)//','//toStr(j))
              else if (j .ne. jcld) then
                shadow_y(imode) = max(shadow_y(imode), d)
              end if
            else
              if (j .ne. jcld) then
@assertEqual(e(:, i, j, iphi, 2), e(:, i, j, iphi, imode), eps, 'open y: shadow wrapped around through the y edge, phi0 '//toStr(phis_diag(iphi))//' mode '//toStr(imode)//' column '//toStr(i)//','//toStr(j))
              else if (i .ne. icld) then
                shadow_x(imode) = max(shadow_x(imode), d)
              end if
            end if
          end do
        end do
      end do

      ! and the direction that is not open has to be periodic
      do imode = 3, Nmodes
        print *, 'phi0', phis_diag(iphi), 'mode', imode, &
          & 'max shadow that wrapped around in x', shadow_x(imode), 'y', shadow_y(imode)
        if (lxs(imode)) then
@assertTrue(shadow_y(imode) .gt. incSolar * 1e-2_ireals, 'open x: expected the shadow to wrap around in y, phi0 '//toStr(phis_diag(iphi))//' mode '//toStr(imode))
        else
@assertTrue(shadow_x(imode) .gt. incSolar * 1e-2_ireals, 'open y: expected the shadow to wrap around in x, phi0 '//toStr(phis_diag(iphi))//' mode '//toStr(imode))
        end if
      end do
@assertTrue(maxval(abs(e(:, :, :, iphi, 2) - e(:, :, :, iphi, 1))) .gt. incSolar * 1e-2_ireals, 'expected the shadow to wrap around with periodic boundaries, phi0 '//toStr(phis_diag(iphi)))
    end do
  end subroutine

  ! Open boundaries have to behave as if the edge columns of the domain, including their terrain, continue outwards forever.
  ! We check that with a crater, built from buildings, whose rim sits right at the domain edges:
  !   the open bc solution has to match the one where the crater is put in the middle of a much larger periodic domain.
  @test(npes=[4, 2, 1])
  subroutine test_open_bc_crater_matches_embedded_domain(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: eps = incSolar * 1e-3_ireals
    real(ireals), parameter :: eps_dir = 1e-3_ireals ! threshold on the horizontal components of the sun vector
    ! bounds on the deviations for suns that are not aligned with the grid, see the comment at the checks below
    real(ireals), parameter :: diag_max_diff = incSolar*.25_ireals
    real(ireals), parameter :: diag_mean_diff = incSolar * 5e-3_ireals
    real(ireals), parameter :: diag_frac_diff = .2_ireals

    real(ireals), allocatable, dimension(:, :, :) :: edir, abso
    real(ireals), dimension(Nz + 1, Ncrater, Ncrater, size(phis)) :: edir_open, edir_embedded, edir_embedded_wide, edir_periodic
    real(ireals), dimension(Nz + 1, Ncrater, Ncrater, size(phis)) :: edir_open_2d
    real(ireals), dimension(Nz, Ncrater, Ncrater, size(phis)) :: abso_open, abso_open_2d, abso_embedded
    real(ireals), dimension(6, Nz, Ncrater, Ncrater, size(phis)) :: bedir_open, bedir_embedded, bedir_embedded_wide, bedir_periodic
    real(ireals), dimension(6, Nz, Ncrater, Ncrater, size(phis)) :: bedir_open_2d
    integer(mpiint), dimension(6, Nz, Ncrater, Ncrater, size(phis)) :: bcount_open, bcount_embedded
    real(ireals) :: eps_abso, max_diff_periodic, mean_diff, max_diff, frac_diff, bmean_diff, bmax_diff, d
    real(ireals) :: edir_1d(Nz + 1), sundir(3)
    integer(iintegers) :: iphi, i, j, k, h, faceid, Nsensitive, Nsamples, Ndiff, Nbfaces, Nshadow, Nlit, Nwall
    integer(mpiint) :: comm, myid, expected_count
    logical :: lupwind_wall

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    edir_open = -one
    edir_embedded = -one
    edir_embedded_wide = -one
    edir_periodic = -one
    abso_open = -one
    abso_embedded = -one
    edir_open_2d = -one
    abso_open_2d = -one

    do iphi = 1, size(phis)
      call solve_crater(comm, phis(iphi), lopen_bc=.true., pad=0_iintegers, &
        & edir=edir, abso=abso, bedir=bedir_open_2d(:, :, :, :, iphi), l2d=.true.)
      if (allocated(edir)) edir_open_2d(:, :, :, iphi) = edir
      if (allocated(abso)) abso_open_2d(:, :, :, iphi) = abso
      call solve_crater(comm, phis(iphi), lopen_bc=.true., pad=0_iintegers, &
        & edir=edir, abso=abso, bedir=bedir_open(:, :, :, :, iphi), bcount=bcount_open(:, :, :, :, iphi))
      if (allocated(edir)) edir_open(:, :, :, iphi) = edir
      if (allocated(abso)) abso_open(:, :, :, iphi) = abso
      call solve_crater(comm, phis(iphi), lopen_bc=.false., pad=Npad, &
        & edir=edir, abso=abso, bedir=bedir_embedded(:, :, :, :, iphi), bcount=bcount_embedded(:, :, :, :, iphi))
      if (allocated(edir)) edir_embedded(:, :, :, iphi) = edir
      if (allocated(abso)) abso_embedded(:, :, :, iphi) = abso
      call solve_crater(comm, phis(iphi), lopen_bc=.false., pad=Npad_wide, &
        & edir=edir, abso=abso, bedir=bedir_embedded_wide(:, :, :, :, iphi))
      if (allocated(edir)) edir_embedded_wide(:, :, :, iphi) = edir
      call solve_crater(comm, phis(iphi), lopen_bc=.false., pad=0_iintegers, &
        & edir=edir, abso=abso, bedir=bedir_periodic(:, :, :, :, iphi))
      if (allocated(edir)) edir_periodic(:, :, :, iphi) = edir
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    print *, 'crater terrain height'
    do j = Ncrater, 1, -1
      print '(*(i4))', (crater_height(i, j), i=1, Ncrater)
    end do
    do iphi = 1, size(phis)
      print *, 'phi0', phis(iphi), 'edir on the terrain: open | embedded | periodic'
      do j = Ncrater, 1, -1
        print '(*(i4))', &
          & (int(edir_open(Nz - crater_height(i, j) + 1, i, j, iphi)), i=1, Ncrater), -1, &
          & (int(edir_embedded(Nz - crater_height(i, j) + 1, i, j, iphi)), i=1, Ncrater), -1, &
          & (int(edir_periodic(Nz - crater_height(i, j) + 1, i, j, iphi)), i=1, Ncrater)
      end do
    end do

    ! results have to be there and valid: the arrays were initialized with -1 and are only overwritten if we got a result
    @assertTrue(all(ieee_is_finite(edir_open)), 'open bc edir is not finite')
    @assertTrue(all(ieee_is_finite(abso_open)), 'open bc absorption is not finite')
    @assertTrue(all(ieee_is_finite(bedir_open)), 'open bc edir on building faces is not finite')
    @assertTrue(all(edir_open .ge. zero), 'missing or negative open bc edir')
    @assertTrue(all(edir_embedded .ge. zero), 'missing or negative embedded edir')
    @assertTrue(all(edir_embedded_wide .ge. zero), 'missing or negative embedded edir of the wider padding')
    @assertTrue(all(edir_periodic .ge. zero), 'missing or negative periodic edir')
    @assertTrue(all(abso_open .ge. zero), 'missing or negative open bc absorption')
    @assertTrue(all(ieee_is_finite(edir_open_2d)), '2D open bc edir is not finite')
    @assertTrue(all(ieee_is_finite(bedir_open_2d)), '2D open bc edir on building faces is not finite')
    @assertTrue(all(edir_open_2d .ge. zero), 'missing or negative 2D open bc edir')
    @assertTrue(all(abso_open_2d .ge. zero), 'missing or negative 2D open bc absorption')
    @assertTrue(all(abso_embedded .ge. zero), 'missing or negative embedded absorption')

    ! each terrain face has to be reported by exactly one rank, and there are no other faces
    do j = 1, Ncrater
      do i = 1, Ncrater
        h = crater_height(i, j)
        do k = 1, Nz
          expected_count = merge(1_mpiint, 0_mpiint, k .ge. Nz - h + 1)
          do faceid = 1, 6
@assertEqual(expected_count, bcount_open(faceid, k, i, j, 1), 'open bc: wrong number of results for building face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
@assertEqual(expected_count, bcount_embedded(faceid, k, i, j, 1), 'embedded: wrong number of results for building face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
          end do
        end do
      end do
    end do
    do iphi = 2, size(phis)
@assertTrue(all(bcount_open(:, :, :, :, iphi) .eq. bcount_open(:, :, :, :, 1)), 'building face ownership changed with the sun angle')
@assertTrue(all(bcount_embedded(:, :, :, :, iphi) .eq. bcount_embedded(:, :, :, :, 1)), 'building face ownership changed with the sun angle')
    end do

    ! the embedded crater is only a valid reference if the padding is wide enough, i.e. the result must not depend on it
    do iphi = 1, size(phis)
@assertEqual(edir_embedded(:, :, :, iphi), edir_embedded_wide(:, :, :, iphi), eps, 'embedded crater edir depends on the width of the padding, phi0 '//toStr(phis(iphi)))
@assertEqual(bedir_embedded(:, :, :, :, iphi), bedir_embedded_wide(:, :, :, :, iphi), eps, 'embedded crater edir on building faces depends on the width of the padding, phi0 '//toStr(phis(iphi)))
    end do

    ! the setup has to be sensitive to the boundary condition, i.e. with periodic boundaries the shadows of the rim wrap around.
    ! This does not hold for each azimuth, e.g. if the low rim is hidden behind the high rim on the opposite side
    Nsensitive = 0
    do iphi = 1, size(phis)
      max_diff_periodic = maxval(abs(edir_periodic(:, :, :, iphi) - edir_embedded(:, :, :, iphi)))
      print *, 'phi0', phis(iphi), 'max edir difference between periodic and embedded crater', max_diff_periodic
      print *, 'phi0', phis(iphi), 'max edir difference between open bc and embedded crater', &
        & maxval(abs(edir_open(:, :, :, iphi) - edir_embedded(:, :, :, iphi))), &
        & 'on terrain faces', maxval(abs(bedir_open(:, :, :, :, iphi) - bedir_embedded(:, :, :, :, iphi)))
      if (max_diff_periodic .gt. incSolar * 1e-2_ireals) Nsensitive = Nsensitive + 1
    end do
@assertTrue(Nsensitive .ge. size(phis) / 2, 'expected the shadows of the rim to wrap around with periodic boundaries for most azimuths')

    edir_1d = edir_clearsky()

    do iphi = 1, size(phis)
      sundir = spherical_2_cartesian(phis(iphi), theta0) ! points away from the sun

      ! the rim has to cast a shadow onto the terrain free crater floor while other parts of the floor are lit,
      ! otherwise this test is void
      Nshadow = 0
      Nlit = 0
      do j = 1, Ncrater
        do i = 1, Ncrater
          if (crater_height(i, j) .ne. 0) cycle
          if (edir_embedded(Nz + 1, i, j, iphi) .lt. edir_1d(Nz + 1)*.25_ireals) Nshadow = Nshadow + 1
          if (edir_embedded(Nz + 1, i, j, iphi) .gt. edir_1d(Nz + 1)*.75_ireals) Nlit = Nlit + 1
        end do
      end do
      print *, 'phi0', phis(iphi), 'crater floor cells in the shadow', Nshadow, 'lit', Nlit
      @assertTrue(Nshadow .gt. 0, 'expected the rim to cast a shadow onto the crater floor, phi0 '//toStr(phis(iphi)))
      @assertTrue(Nlit .gt. 0, 'expected parts of the crater floor to be lit, phi0 '//toStr(phis(iphi)))

      ! The inflow of an edge column is computed as if that column, including its terrain, is repeated in x and y,
      ! i.e. nothing enters through the outer walls of the terrain at the sunward edges. This holds for all sun angles
      Nwall = 0
      do j = 1, Ncrater
        do i = 1, Ncrater
          h = crater_height(i, j)
          do k = Nz - h + 1, Nz
            do faceid = 3, 6
              select case (faceid)
              case (3) ! x-low
                lupwind_wall = i .eq. 1 .and. sundir(1) .gt. eps_dir
              case (4) ! x-high
                lupwind_wall = i .eq. Ncrater .and. sundir(1) .lt. -eps_dir
              case (5) ! y-low
                lupwind_wall = j .eq. 1 .and. sundir(2) .gt. eps_dir
              case default ! y-high
                lupwind_wall = j .eq. Ncrater .and. sundir(2) .lt. -eps_dir
              end select
              if (.not. lupwind_wall) cycle
              Nwall = Nwall + 1
@assertEqual(zero, bedir_open(faceid, k, i, j, iphi), incSolar * 1e-8_ireals, 'expected no inflow into the outer wall of the terrain, phi0 '//toStr(phis(iphi))//' face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
@assertEqual(zero, bedir_open_2d(faceid, k, i, j, iphi), incSolar * 1e-8_ireals, '2D open bc: expected no inflow into the outer wall of the terrain, phi0 '//toStr(phis(iphi))//' face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
            end do
          end do
        end do
      end do
      @assertTrue(Nwall .gt. 0, 'expected terrain at the sunward edges, phi0 '//toStr(phis(iphi)))

      ! With the 2D inflow, the edge columns see their neighbours along the edge,
      ! i.e. the solution has to match the embedded crater for all sun angles
      eps_abso = maxval(abso_embedded(:, :, :, iphi)) * 1e-3_ireals
      max_diff = zero
      do j = 1, Ncrater
        do i = 1, Ncrater
          h = crater_height(i, j)
          do k = 1, Nz - h + 1
            max_diff = max(max_diff, abs(edir_embedded(k, i, j, iphi) - edir_open_2d(k, i, j, iphi)))
@assertEqual(edir_embedded(k, i, j, iphi), edir_open_2d(k, i, j, iphi), eps, '2D open bc crater edir differs from the embedded crater, phi0 '//toStr(phis(iphi))//' level '//toStr(k)//' column '//toStr(i)//','//toStr(j))
          end do
          do k = 1, Nz - h
@assertEqual(abso_embedded(k, i, j, iphi), abso_open_2d(k, i, j, iphi), eps_abso, '2D open bc crater absorption differs from the embedded crater, phi0 '//toStr(phis(iphi))//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
          end do
          do k = Nz - h + 1, Nz
            do faceid = 1, 6
@assertEqual(bedir_embedded(faceid, k, i, j, iphi), bedir_open_2d(faceid, k, i, j, iphi), eps, '2D open bc crater edir on building face differs from the embedded crater, phi0 '//toStr(phis(iphi))//' face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
            end do
          end do
        end do
      end do
      print *, 'phi0', phis(iphi), '2D open bc vs embedded crater above terrain: max', max_diff

      ! deviations from the embedded crater with the inflow from single column solves
      Nsamples = 0
      Ndiff = 0
      max_diff = zero
      mean_diff = zero
      Nbfaces = 0
      bmax_diff = zero
      bmean_diff = zero
      do j = 1, Ncrater
        do i = 1, Ncrater
          h = crater_height(i, j)
          do k = 1, Nz - h + 1 ! atmosphere above the terrain
            d = abs(edir_open(k, i, j, iphi) - edir_embedded(k, i, j, iphi))
            Nsamples = Nsamples + 1
            if (d .gt. eps) Ndiff = Ndiff + 1
            max_diff = max(max_diff, d)
            mean_diff = mean_diff + d
          end do
          do k = Nz - h + 1, Nz
            do faceid = 1, 6
              d = abs(bedir_open(faceid, k, i, j, iphi) - bedir_embedded(faceid, k, i, j, iphi))
              Nbfaces = Nbfaces + 1
              bmax_diff = max(bmax_diff, d)
              bmean_diff = bmean_diff + d
            end do
          end do
        end do
      end do
      mean_diff = mean_diff / real(Nsamples, ireals)
      bmean_diff = bmean_diff / real(Nbfaces, ireals)
      frac_diff = real(Ndiff, ireals) / real(Nsamples, ireals)
      print *, 'phi0', phis(iphi), 'open bc vs embedded crater above terrain: max', max_diff, 'mean', mean_diff, &
        & 'fraction off', frac_diff, 'on terrain faces: max', bmax_diff, 'mean', bmean_diff

      ! If the sun is not aligned with the grid, the inflow from the single column solves misses the shadows
      ! that the terrain beyond the neighbouring edge columns would cast,
      ! i.e. where terrain height varies along the inflow edge, we can only expect to roughly match the embedded crater.
      ! The bounds pin the current quality of the approximation: deviations have to stay local and bounded
      if (.not. modulo(nint(phis(iphi)), 90) .eq. 0) then
@assertTrue(max_diff .lt. diag_max_diff, 'open bc crater edir locally differs too much from the embedded crater, phi0 '//toStr(phis(iphi)))
@assertTrue(mean_diff .lt. diag_mean_diff, 'open bc crater edir differs too much from the embedded crater, phi0 '//toStr(phis(iphi)))
@assertTrue(frac_diff .lt. diag_frac_diff, 'open bc crater edir differs from the embedded crater in too many places, phi0 '//toStr(phis(iphi)))
@assertTrue(bmax_diff .lt. diag_max_diff, 'open bc crater edir on building faces locally differs too much from the embedded crater, phi0 '//toStr(phis(iphi)))
@assertTrue(bmean_diff .lt. diag_mean_diff, 'open bc crater edir on building faces differs too much from the embedded crater, phi0 '//toStr(phis(iphi)))
        cycle
      end if

      eps_abso = maxval(abso_embedded(:, :, :, iphi)) * 1e-3_ireals
      do j = 1, Ncrater
        do i = 1, Ncrater
          h = crater_height(i, j)
          ! atmosphere above the terrain
          do k = 1, Nz - h + 1
@assertEqual(edir_embedded(k, i, j, iphi), edir_open(k, i, j, iphi), eps, 'open bc crater edir differs from the embedded crater, phi0 '//toStr(phis(iphi))//' level '//toStr(k)//' column '//toStr(i)//','//toStr(j))
          end do
          do k = 1, Nz - h
@assertEqual(abso_embedded(k, i, j, iphi), abso_open(k, i, j, iphi), eps_abso, 'open bc crater absorption differs from the embedded crater, phi0 '//toStr(phis(iphi))//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
          end do
          ! direct radiation on the terrain faces
          do k = Nz - h + 1, Nz
            do faceid = 1, 6
@assertEqual(bedir_embedded(faceid, k, i, j, iphi), bedir_open(faceid, k, i, j, iphi), eps, 'open bc crater edir on building face differs from the embedded crater, phi0 '//toStr(phis(iphi))//' face '//toStr(faceid)//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
            end do
          end do
        end do
      end do
    end do
  end subroutine
end module
