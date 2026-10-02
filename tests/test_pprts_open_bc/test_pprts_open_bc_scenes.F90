module test_pprts_open_bc_scenes
  ! Further tests for the open boundary conditions (-pprts_open_bc) on small, deliberately irregular scenes:
  !   * rectangular domain (Nx /= Ny) with non uniform layer thickness (dz /= dx)
  !   * each column has its own absorption profile
  !   * some of the edge columns, including two corners, carry an elevated opaque slab made of buildings
  !
  ! All tests run twice, with the 2D inflow (-pprts_open_bc_2d, the default) and with the inflow from single column solves
  !
  ! Note on the structure of the tests:
  !   results are gathered on rank 0 and can only be checked there.
  !   A failing pFUnit assert returns from the test, so we must not assert in between the (collective) solves
  !   or the other ranks deadlock. Hence, first do all solves, then check.

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  use pfunit_mod

  use m_helper_functions, only: &
    & CHKERR, &
    & deg2rad, &
    & insert_petsc_opt, &
    & spherical_2_cartesian, &
    & toStr

  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
  use m_pprts, only: init_pprts, set_angles, set_optical_properties, solve_pprts, pprts_get_result_toZero

  use m_buildings, only: &
    & check_buildings_consistency, &
    & destroy_buildings, &
    & faceidx_by_cell_plus_offset, &
    & init_buildings, &
    & t_pprts_buildings

  implicit none

  integer(iintegers), parameter :: Nz = 6
  real(ireals), parameter :: dx = 100
  ! all layers are solved in 3D
  real(ireals), parameter :: dz_3d(Nz) = [real(ireals) :: 80, 120, 100, 100, 90, 110]
  ! the two uppermost layers exceed the aspect ratio of -twostr_ratio (default 2) and are solved in 1D
  real(ireals), parameter :: dz_mixed(Nz) = [real(ireals) :: 250, 300, 100, 100, 90, 110]
  real(ireals), parameter :: tau_clearsky = .3_ireals ! vertically integrated optical thickness of the homogeneous scenes
  integer(iintegers), parameter :: kslab_top = 3, kslab_bot = 4 ! layers that are occupied by the slabs

  type t_cfg
    integer(iintegers) :: Nx = 6, Ny = 5         ! global domain size
    real(ireals) :: dz(Nz) = dz_3d
    real(ireals) :: phi0 = 0, theta0 = 60
    real(ireals) :: edirTOA = 1000
    real(ireals) :: albedo = 0
    real(ireals) :: w0 = 0, g = 0                ! single scattering albedo and asymmetry parameter
    real(ireals) :: kabs_scale = 1
    logical :: lopen_bc = .true.
    logical :: l2d = .false.                     ! -pprts_open_bc_2d, i.e. zero gradient inflow instead of single column solves
    logical :: lheterogeneous = .true.           ! each column gets its own absorption profile, otherwise homogeneous
    logical :: lslabs = .true.                   ! put opaque slabs into some of the edge columns
    integer(iintegers) :: icol = -1, jcol = -1   ! if set, replicate this global column all over the domain
    character(len=8) :: solvername = '3_10'
    integer(iintegers), allocatable :: nxproc(:), nyproc(:) ! if allocated, force the domain decomposition
  end type

  type t_res
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! global, only on rank 0
    logical, allocatable :: l1d(:)
    integer(iintegers) :: Ndir_streams = -1
  end type

  ! set by the solve helpers if a solver did not pick up the requested -pprts_open_bc options
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
    call set_open_bc_option(.false., .false.)
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  ! absorption coefficient of the cell in layer k of the global column i,j (indices start at 1)
  pure function scene_kabs(cfg, k, i, j) result(kabs)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: k, i, j
    real(ireals) :: kabs
    kabs = cfg%kabs_scale * tau_clearsky / sum(cfg%dz)
    if (cfg%lheterogeneous) kabs = kabs * (one + real(modulo(3 * i + 5 * j + k, 4_iintegers), ireals) / 2)
  end function

  ! does the global column i,j carry a slab?
  ! Slabs sit on all four edges, in two of the four corners, and one edge column next to a slab is always free
  pure function scene_has_slab(cfg, i, j) result(lslab)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: i, j
    logical :: lslab
    lslab = .false.
    if (.not. cfg%lslabs) return
    lslab = (i .eq. 1 .and. j .eq. 1) &
      & .or. (i .eq. cfg%Nx .and. j .eq. cfg%Ny) &
      & .or. (i .eq. 1 .and. j .eq. cfg%Ny - 1) &
      & .or. (i .eq. cfg%Nx .and. j .eq. 2) &
      & .or. (i .eq. 3 .and. j .eq. 1) &
      & .or. (i .eq. cfg%Nx - 2 .and. j .eq. cfg%Ny)
  end function

  ! global column that provides the properties for the local column i,j (indices start at 1)
  pure subroutine source_column(cfg, i, j, isrc, jsrc)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: i, j
    integer(iintegers), intent(out) :: isrc, jsrc
    if (cfg%icol .gt. 0) then
      isrc = cfg%icol
      jsrc = cfg%jcol
    else
      isrc = i
      jsrc = j
    end if
  end subroutine

  ! always set all of the options, they stay in the options database
  subroutine set_open_bc_option(lopen_bc, l2d)
    logical, intent(in) :: lopen_bc, l2d
    integer(mpiint) :: ierr
    call insert_petsc_opt('-pprts_open_bc '//merge('yes', 'no ', lopen_bc), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_x '//merge('yes', 'no ', lopen_bc), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_y '//merge('yes', 'no ', lopen_bc), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_2d '//merge('yes', 'no ', l2d), ierr); call CHKERR(ierr)
  end subroutine

  subroutine init_scene(comm, cfg, solver)
    integer(mpiint), intent(in) :: comm
    type(t_cfg), intent(in) :: cfg
    class(t_solver), allocatable, intent(inout) :: solver
    integer(mpiint) :: ierr

    call set_open_bc_option(cfg%lopen_bc, cfg%l2d)
    call allocate_pprts_solver_from_commandline(solver, trim(cfg%solvername), ierr); call CHKERR(ierr)

    if (allocated(cfg%nxproc)) then
      call init_pprts(comm, Nz, cfg%Nx, cfg%Ny, dx, dx, spherical_2_cartesian(cfg%phi0, cfg%theta0), solver, &
        & dz1d=cfg%dz, nxproc=cfg%nxproc, nyproc=cfg%nyproc)
    else
      call init_pprts(comm, Nz, cfg%Nx, cfg%Ny, dx, dx, spherical_2_cartesian(cfg%phi0, cfg%theta0), solver, dz1d=cfg%dz)
    end if
    if (cfg%lopen_bc .neqv. solver%lopen_bc) loption_mismatch = .true.
    if (cfg%lopen_bc .neqv. solver%lopen_bc_x) loption_mismatch = .true.
    if (cfg%lopen_bc .neqv. solver%lopen_bc_y) loption_mismatch = .true.
    if (cfg%l2d .neqv. solver%lopen_bc_2d) loption_mismatch = .true.
  end subroutine

  ! set the optical properties and buildings of the scene, solve and gather the results on rank 0.
  ! The sun angles are the ones the solver currently has, i.e. from init_pprts or set_angles
  subroutine run_scene(solver, cfg, res, res_again)
    class(t_solver), allocatable, intent(inout) :: solver
    type(t_cfg), intent(in) :: cfg
    type(t_res), intent(out) :: res
    type(t_res), intent(out), optional :: res_again ! results gathered a second time from the same solution

    type(t_pprts_buildings), allocatable :: buildings
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g
    integer(iintegers) :: Nfaces, m, i, j, k, isrc, jsrc, faceid
    logical :: lbuildings
    integer(mpiint) :: ierr

    ! if there are slabs anywhere in the domain, has to be the same on all ranks
    lbuildings = cfg%lslabs
    if (cfg%icol .gt. 0) lbuildings = scene_has_slab(cfg, cfg%icol, cfg%jcol)

    associate (C => solver%C_one)
      allocate (kabs(C%zm, C%xm, C%ym))
      allocate (ksca(C%zm, C%xm, C%ym))
      allocate (g(C%zm, C%xm, C%ym), source=cfg%g)

      Nfaces = 0
      do j = C%ys, C%ye
        do i = C%xs, C%xe
          call source_column(cfg, i + 1, j + 1, isrc, jsrc)
          do k = 1, Nz
            kabs(k, i - C%xs + 1, j - C%ys + 1) = scene_kabs(cfg, k, isrc, jsrc)
          end do
          if (scene_has_slab(cfg, isrc, jsrc)) Nfaces = Nfaces + 6 * (kslab_bot - kslab_top + 1)
        end do
      end do
      ksca = kabs * cfg%w0 / (one - cfg%w0)

      if (lbuildings) then
        call init_buildings(buildings, [integer(iintegers) :: 6, C%zm, C%xm, C%ym], Nfaces, ierr); call CHKERR(ierr)
        m = 0
        do j = C%ys, C%ye
          do i = C%xs, C%xe
            call source_column(cfg, i + 1, j + 1, isrc, jsrc)
            if (.not. scene_has_slab(cfg, isrc, jsrc)) cycle
            do k = kslab_top, kslab_bot
              do faceid = 1, 6
                m = m + 1
                buildings%iface(m) = faceidx_by_cell_plus_offset(buildings%da_offsets, k, i - C%xs + 1, j - C%ys + 1, faceid)
                buildings%albedo(m) = cfg%albedo
              end do
            end do
          end do
        end do
        call check_buildings_consistency(buildings, C%zm, C%xm, C%ym, ierr); call CHKERR(ierr)
      end if
    end associate

    call set_optical_properties(solver, cfg%albedo, kabs, ksca, g)
    if (lbuildings) then
      call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=cfg%edirTOA, opt_buildings=buildings)
    else
      call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=cfg%edirTOA)
    end if

    call pprts_get_result_toZero(solver, res%edn, res%eup, res%abso, res%edir)
    if (present(res_again)) call pprts_get_result_toZero(solver, res_again%edn, res_again%eup, res_again%abso, res_again%edir)

    allocate (res%l1d(Nz))
    res%l1d = solver%atm%l1d
    res%Ndir_streams = solver%dirtop%dof + solver%dirside%dof * 2

    if (lbuildings) then
      call destroy_buildings(buildings, ierr); call CHKERR(ierr)
    end if
  end subroutine

  ! solve the scene with a fresh solver
  subroutine solve_scene(comm, cfg, res)
    integer(mpiint), intent(in) :: comm
    type(t_cfg), intent(in) :: cfg
    type(t_res), intent(out) :: res
    class(t_solver), allocatable :: solver

    call init_scene(comm, cfg, solver)
    call run_scene(solver, cfg, res)
    call destroy_pprts(solver, lfinalizepetsc=.false.)
  end subroutine

  ! The inflow from the single column solves is a fixed boundary condition, the explicit solver is done after one sweep.
  ! The 2D inflow is part of the iterative solution and only as good as the convergence criteria of the solver
  pure function get_tol_scale(l2d) result(tol_scale)
    logical, intent(in) :: l2d
    real(ireals) :: tol_scale
    tol_scale = merge(10._ireals, one, l2d)
  end function

  ! global columns along the sunward edges of the domain whose complete inflow is given by the open boundary condition:
  !   * if the sun is aligned with the grid, all columns of the sunward edge
  !   * otherwise only the sunward corner, the other edge columns also get radiation from their neighbours
  subroutine inflow_only_columns(cfg, icols, jcols)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), allocatable, intent(out) :: icols(:), jcols(:)
    real(ireals), parameter :: eps = 1e-3_ireals
    real(ireals) :: sundir(3)
    integer(iintegers) :: iedge, jedge, i

    sundir = spherical_2_cartesian(cfg%phi0, cfg%theta0) ! points away from the sun

    iedge = -1
    if (sundir(1) .gt. eps) iedge = 1
    if (sundir(1) .lt. -eps) iedge = cfg%Nx
    jedge = -1
    if (sundir(2) .gt. eps) jedge = 1
    if (sundir(2) .lt. -eps) jedge = cfg%Ny

    if (iedge .gt. 0 .and. jedge .gt. 0) then
      icols = [iedge]
      jcols = [jedge]
    else if (iedge .gt. 0) then
      icols = [(iedge, i=1, cfg%Ny)]
      jcols = [(i, i=1, cfg%Ny)]
    else if (jedge .gt. 0) then
      icols = [(i, i=1, cfg%Nx)]
      jcols = [(jedge, i=1, cfg%Nx)]
    else
      allocate (icols(0), jcols(0))
    end if
  end subroutine

  ! The open boundaries feed each sunward edge column with the side fluxes it would get if this very column,
  ! including its buildings, is repeated all over. For the columns that get all of their inflow from the boundary condition,
  ! the solution therefore has to be the same as in a periodic domain that is filled with copies of that column.
  ! The reference solves do not use the open boundary code at all
  subroutine check_open_bc_edge_columns_match_replicated_column(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d
    real(ireals) :: tol_scale

    real(ireals), parameter :: phis(12) = [real(ireals) :: 0, 90, 180, 270, 45, 135, 225, 315, 20, 110, 200, 290]
    real(ireals), parameter :: thetas(2) = [real(ireals) :: 60, 25]
    integer(iintegers), parameter :: Nmax = 6 ! max number of columns to check per sun angle

    type(t_cfg) :: cfg, cfg_ref
    type(t_res) :: res
    real(ireals), dimension(Nz + 1, Nmax, size(phis), size(thetas)) :: edir_open, edir_ref
    real(ireals), dimension(Nz, Nmax, size(phis), size(thetas)) :: abso_open, abso_ref
    integer(iintegers), dimension(Nmax, size(phis), size(thetas)) :: ic, jc
    integer(iintegers) :: Ncols(size(phis), size(thetas))
    integer(iintegers), allocatable :: icols(:), jcols(:)
    integer(iintegers) :: iphi, itheta, m, k, Nslab, Nfree
    real(ireals) :: eps, eps_abso, mu
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()
    cfg%l2d = l2d
    tol_scale = get_tol_scale(l2d)

    edir_open = -one
    edir_ref = -one
    abso_open = -one
    abso_ref = -one

    do itheta = 1, size(thetas)
      do iphi = 1, size(phis)
        cfg%phi0 = phis(iphi)
        cfg%theta0 = thetas(itheta)
        cfg%lopen_bc = .true.
        call inflow_only_columns(cfg, icols, jcols)
        Ncols(iphi, itheta) = size(icols)
        ic(1:size(icols), iphi, itheta) = icols
        jc(1:size(icols), iphi, itheta) = jcols

        call solve_scene(comm, cfg, res)
        do m = 1, size(icols)
          if (allocated(res%edir)) edir_open(:, m, iphi, itheta) = res%edir(:, icols(m), jcols(m))
          if (allocated(res%abso)) abso_open(:, m, iphi, itheta) = res%abso(:, icols(m), jcols(m))
        end do

        do m = 1, size(icols)
          cfg_ref = cfg
          cfg_ref%lopen_bc = .false.
          cfg_ref%icol = icols(m)
          cfg_ref%jcol = jcols(m)
          call solve_scene(comm, cfg_ref, res)
          ! any column will do, pick one in the middle of the domain
          if (allocated(res%edir)) edir_ref(:, m, iphi, itheta) = res%edir(:, cfg%Nx / 2, cfg%Ny / 2)
          if (allocated(res%abso)) abso_ref(:, m, iphi, itheta) = res%abso(:, cfg%Nx / 2, cfg%Ny / 2)
        end do
      end do
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    Nslab = 0
    Nfree = 0
    do itheta = 1, size(thetas)
      mu = cos(deg2rad(thetas(itheta)))
      eps = cfg%edirTOA * mu * 1e-4_ireals * tol_scale
      do iphi = 1, size(phis)
        @assertTrue(Ncols(iphi, itheta) .ge. 1, 'expected at least one column to check')
        do m = 1, Ncols(iphi, itheta)
          msg = 'phi0 '//toStr(phis(iphi))//' theta0 '//toStr(thetas(itheta))// &
            & ' column '//toStr(ic(m, iphi, itheta))//','//toStr(jc(m, iphi, itheta))
          print *, msg, ' slab ', scene_has_slab(cfg, ic(m, iphi, itheta), jc(m, iphi, itheta)), &
            & 'edir open', edir_open(:, m, iphi, itheta), 'ref', edir_ref(:, m, iphi, itheta)

          @assertTrue(all(ieee_is_finite(edir_open(:, m, iphi, itheta))), 'open bc edir is not finite, '//msg)
          @assertTrue(all(edir_ref(:, m, iphi, itheta) .ge. zero), 'missing or negative reference edir, '//msg)
          @assertTrue(all(abso_ref(:, m, iphi, itheta) .ge. zero), 'missing or negative reference absorption, '//msg)

@assertEqual(edir_ref(:, m, iphi, itheta), edir_open(:, m, iphi, itheta), eps, 'open bc edir differs from the replicated column, '//msg)
          eps_abso = maxval(abso_ref(:, m, iphi, itheta)) * 1e-3_ireals * tol_scale
@assertEqual(abso_ref(:, m, iphi, itheta), abso_open(:, m, iphi, itheta), eps_abso, 'open bc absorption differs from the replicated column, '//msg)

          ! the comparison is only worth something if there is light, it must not be dark above the slabs or in free columns
          do k = 1, kslab_top
@assertTrue(edir_open(k, m, iphi, itheta) .gt. cfg%edirTOA * mu * .5_ireals, 'expected light above the slab, '//msg//' level '//toStr(k))
          end do
          if (scene_has_slab(cfg, ic(m, iphi, itheta), jc(m, iphi, itheta))) then
            Nslab = Nslab + 1
            ! nothing must enter through the open boundary below a slab that continues outwards forever
            do k = kslab_bot + 1, Nz + 1
@assertEqual(zero, edir_open(k, m, iphi, itheta), cfg%edirTOA * 1e-8_ireals, 'expected no light below the slab, '//msg//' level '//toStr(k))
            end do
          else
            Nfree = Nfree + 1
           @assertTrue(edir_open(Nz + 1, m, iphi, itheta) .gt. cfg%edirTOA * mu * .2_ireals, 'expected light at the surface, '//msg)
          end if
        end do
      end do
    end do
    @assertTrue(Nslab .ge. 8, 'expected to check columns with slabs')
    @assertTrue(Nfree .ge. 8, 'expected to check columns without slabs')
  end subroutine

  ! In a horizontally homogeneous atmosphere, open and periodic boundaries have to give the same answer.
  ! Run this for various zenith and azimuth angles that are off the symmetry axes and with 1D layers on top of 3D layers.
  ! The purely absorbing atmosphere additionally allows to check against Beer-Lambert and to close the energy budget per layer
  subroutine check_open_bc_homogeneous_angles_layers_and_solvers(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d

    real(ireals), parameter :: phis(5) = [real(ireals) :: 20, 110, 200, 290, 90]
    real(ireals), parameter :: thetas(3) = [real(ireals) :: 0, 30, 75]
    ! only the LUTs of the 3_10 solver are part of the commonly available ones, add further solvers here if you have their LUTs
    character(len=8), parameter :: solvers(1) = [character(len=8) :: '3_10']
    integer(iintegers), parameter :: Ndir_streams(1) = [3]
    integer(iintegers), parameter :: Ndz = 2
    real(ireals), parameter :: eps_1d = 5e-2_ireals ! relative accuracy of the transport coefficients vs Beer-Lambert

    type(t_cfg) :: cfg
    type(t_res) :: res
    real(ireals), allocatable, dimension(:, :, :, :, :, :, :) :: edir_open, edir_periodic, abso_open, abso_periodic
    logical :: l1d(Nz, Ndz, size(solvers))
    integer(iintegers) :: Nstreams(Ndz, size(solvers))
    real(ireals) :: edir_1d(Nz + 1), mu, eps, eps_abso, tau
    integer(iintegers) :: iphi, itheta, idz, isolver, i, j, k
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()
    cfg%l2d = l2d

    allocate (edir_open(Nz + 1, cfg%Nx, cfg%Ny, size(phis), size(thetas), Ndz, size(solvers)), source=-one)
    allocate (edir_periodic, source=edir_open)
    allocate (abso_open(Nz, cfg%Nx, cfg%Ny, size(phis), size(thetas), Ndz, size(solvers)), source=-one)
    allocate (abso_periodic, source=abso_open)

    cfg%lheterogeneous = .false.
    cfg%lslabs = .false.

    do isolver = 1, size(solvers)
      cfg%solvername = solvers(isolver)
      do idz = 1, Ndz
        cfg%dz = merge(dz_3d, dz_mixed, idz .eq. 1)
        do itheta = 1, size(thetas)
          cfg%theta0 = thetas(itheta)
          do iphi = 1, size(phis)
            cfg%phi0 = phis(iphi)

            cfg%lopen_bc = .true.
            call solve_scene(comm, cfg, res)
            if (allocated(res%edir)) edir_open(:, :, :, iphi, itheta, idz, isolver) = res%edir
            if (allocated(res%abso)) abso_open(:, :, :, iphi, itheta, idz, isolver) = res%abso
            l1d(:, idz, isolver) = res%l1d
            Nstreams(idz, isolver) = res%Ndir_streams

            cfg%lopen_bc = .false.
            call solve_scene(comm, cfg, res)
            if (allocated(res%edir)) edir_periodic(:, :, :, iphi, itheta, idz, isolver) = res%edir
            if (allocated(res%abso)) abso_periodic(:, :, :, iphi, itheta, idz, isolver) = res%abso
          end do
        end do
      end do
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do isolver = 1, size(solvers)
      do idz = 1, Ndz
        cfg%dz = merge(dz_3d, dz_mixed, idz .eq. 1)

        ! make sure that we actually ran what we think we ran
@assertEqual(Ndir_streams(isolver), Nstreams(idz, isolver), 'unexpected number of direct streams for solver '//trim(solvers(isolver)))
        if (idz .eq. 1) then
          @assertFalse(any(l1d(:, idz, isolver)), 'did not expect 1D layers')
        else
          @assertTrue(all(l1d(1:2, idz, isolver)), 'expected the two uppermost layers to be 1D layers')
          @assertFalse(any(l1d(3:Nz, idz, isolver)), 'expected the lower layers to be 3D layers')
        end if

        do itheta = 1, size(thetas)
          mu = cos(deg2rad(thetas(itheta)))
          eps = cfg%edirTOA * mu * 1e-4_ireals
          tau = zero
          edir_1d(1) = cfg%edirTOA * mu
          do k = 1, Nz
            tau = tau + scene_kabs(cfg, k, 1_iintegers, 1_iintegers) * cfg%dz(k)
            edir_1d(k + 1) = cfg%edirTOA * mu * exp(-tau / mu)
          end do

          do iphi = 1, size(phis)
         msg = 'solver '//trim(solvers(isolver))//' dz '//toStr(idz)//' theta0 '//toStr(thetas(itheta))//' phi0 '//toStr(phis(iphi))
            associate ( &
                & eo => edir_open(:, :, :, iphi, itheta, idz, isolver), &
                & ep => edir_periodic(:, :, :, iphi, itheta, idz, isolver), &
                & ao => abso_open(:, :, :, iphi, itheta, idz, isolver), &
                & ap => abso_periodic(:, :, :, iphi, itheta, idz, isolver))

              print *, msg, ' edir 1D', edir_1d, 'open min', minval(minval(eo, dim=3), dim=2), 'max', maxval(maxval(eo, &
                                                                                                                    dim=3), dim=2)

              @assertTrue(all(ieee_is_finite(eo)), 'open bc edir is not finite, '//msg)
              @assertTrue(all(ieee_is_finite(ao)), 'open bc absorption is not finite, '//msg)
              @assertTrue(all(ep .ge. zero), 'missing or negative periodic edir, '//msg)
              @assertTrue(all(ap .ge. zero), 'missing or negative periodic absorption, '//msg)

              @assertEqual(ep, eo, eps, 'open bc edir should be the same as with periodic boundaries, '//msg)
              eps_abso = maxval(ap) * 1e-3_ireals
              @assertEqual(ap, ao, eps_abso, 'open bc absorption should be the same as with periodic boundaries, '//msg)

              do j = 1, size(eo, dim=3)
                do i = 1, size(eo, dim=2)
@assertEqual(edir_1d, eo(:, i, j), cfg%edirTOA * mu * eps_1d, 'open bc edir should match Beer-Lambert, '//msg//' column '//toStr(i)//','//toStr(j))
                  ! what is lost from the direct beam has to show up as absorption
                  do k = 1, Nz
@assertEqual(eo(k, i, j) - eo(k + 1, i, j), ao(k, i, j) * cfg%dz(k), eps, 'open bc absorption does not close the energy budget, '//msg//' layer '//toStr(k)//' column '//toStr(i)//','//toStr(j))
                  end do
                end do
              end do
            end associate
          end do
        end do
      end do
    end do
  end subroutine

  ! Diffuse radiation stays periodic but its source terms are computed from the direct radiation,
  ! including the direct radiation that enters through the open boundaries.
  ! In a horizontally homogeneous atmosphere with scattering and a reflecting surface,
  ! open and periodic boundaries therefore have to give the same diffuse fluxes
  subroutine check_open_bc_homogeneous_scattering_matches_periodic(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d

    real(ireals), parameter :: phis(6) = [real(ireals) :: 0, 90, 20, 110, 200, 290]

    type(t_cfg) :: cfg
    type(t_res) :: res
    real(ireals), dimension(Nz + 1, 6, 5, size(phis)) :: edir_open, edir_periodic, edn_open, edn_periodic, eup_open, eup_periodic
    real(ireals), dimension(Nz, 6, 5, size(phis)) :: abso_open, abso_periodic
    real(ireals) :: eps, eps_abso
    integer(iintegers) :: iphi
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()
    cfg%l2d = l2d

    edir_open = -one; edir_periodic = -one
    edn_open = -one; edn_periodic = -one
    eup_open = -one; eup_periodic = -one
    abso_open = -one; abso_periodic = -one

    cfg%lheterogeneous = .false.
    cfg%lslabs = .false.
    cfg%kabs_scale = 3
    cfg%w0 = .7_ireals
    cfg%g = .5_ireals
    cfg%albedo = .3_ireals

    do iphi = 1, size(phis)
      cfg%phi0 = phis(iphi)

      cfg%lopen_bc = .true.
      call solve_scene(comm, cfg, res)
      if (allocated(res%edir)) edir_open(:, :, :, iphi) = res%edir
      if (allocated(res%edn)) edn_open(:, :, :, iphi) = res%edn
      if (allocated(res%eup)) eup_open(:, :, :, iphi) = res%eup
      if (allocated(res%abso)) abso_open(:, :, :, iphi) = res%abso

      cfg%lopen_bc = .false.
      call solve_scene(comm, cfg, res)
      if (allocated(res%edir)) edir_periodic(:, :, :, iphi) = res%edir
      if (allocated(res%edn)) edn_periodic(:, :, :, iphi) = res%edn
      if (allocated(res%eup)) eup_periodic(:, :, :, iphi) = res%eup
      if (allocated(res%abso)) abso_periodic(:, :, :, iphi) = res%abso
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    eps = cfg%edirTOA * 1e-4_ireals
    do iphi = 1, size(phis)
      msg = 'phi0 '//toStr(phis(iphi))
      print *, msg, ' open edir', edir_open(:, 1, 1, iphi), 'edn', edn_open(:, 1, 1, iphi), 'eup', eup_open(:, 1, 1, iphi)
      print *, msg, ' max diff open vs periodic: edir', maxval(abs(edir_open(:, :, :, iphi) - edir_periodic(:, :, :, iphi))), &
        & 'edn', maxval(abs(edn_open(:, :, :, iphi) - edn_periodic(:, :, :, iphi))), &
        & 'eup', maxval(abs(eup_open(:, :, :, iphi) - eup_periodic(:, :, :, iphi)))

      ! the diffuse radiation field has to be substantial, otherwise this test is void
      @assertTrue(minval(eup_periodic(:, :, :, iphi)) .gt. cfg%edirTOA * 1e-2_ireals, 'expected upwelling diffuse radiation, '//msg)
@assertTrue(minval(edn_periodic(Nz + 1, :, :, iphi)) .gt. cfg%edirTOA * 1e-2_ireals, 'expected diffuse radiation at the surface, '//msg)
      @assertTrue(all(abso_periodic(:, :, :, iphi) .gt. zero), 'missing periodic absorption, '//msg)

@assertEqual(edir_periodic(:, :, :, iphi), edir_open(:, :, :, iphi), eps, 'open bc edir should be the same as with periodic boundaries, '//msg)
@assertEqual(edn_periodic(:, :, :, iphi), edn_open(:, :, :, iphi), eps, 'open bc edn should be the same as with periodic boundaries, '//msg)
@assertEqual(eup_periodic(:, :, :, iphi), eup_open(:, :, :, iphi), eps, 'open bc eup should be the same as with periodic boundaries, '//msg)
      eps_abso = maxval(abso_periodic(:, :, :, iphi)) * 1e-3_ireals
@assertEqual(abso_periodic(:, :, :, iphi), abso_open(:, :, :, iphi), eps_abso, 'open bc absorption should be the same as with periodic boundaries, '//msg)
    end do
  end subroutine

  ! The solution must not depend on how the domain is split between the ranks,
  ! in particular not on which ranks own the sunward edges or if a rank touches opposite edges at the same time
  subroutine check_open_bc_is_independent_of_domain_decomposition(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d
    real(ireals) :: tol_scale

    real(ireals), parameter :: phis(6) = [real(ireals) :: 20, 110, 200, 290, 90, 180]
    integer(iintegers), parameter :: Nx = 8, Ny = 6
    integer(iintegers), parameter :: Nlayouts = 3

    type(t_cfg) :: cfg
    type(t_res) :: res
    real(ireals), dimension(Nz + 1, Nx, Ny, size(phis), 0:Nlayouts) :: edir
    real(ireals), dimension(Nz, Nx, Ny, size(phis), 0:Nlayouts) :: abso
    real(ireals) :: eps, eps_abso
    integer(iintegers) :: iphi, ilayout
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid, numnodes

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()
    numnodes = this%getNumProcesses()
    cfg%l2d = l2d
    tol_scale = get_tol_scale(l2d)

    edir = -one
    abso = -one

    do ilayout = 0, Nlayouts
      cfg%Nx = Nx
      cfg%Ny = Ny
      if (allocated(cfg%nxproc)) deallocate (cfg%nxproc)
      if (allocated(cfg%nyproc)) deallocate (cfg%nyproc)
      ! layout 0 is the default decomposition
      if (numnodes .eq. 4) then
        select case (ilayout)
        case (1) ! stripes along x, the two ranks in the middle do not touch the x edges
          cfg%nxproc = [integer(iintegers) :: 3, 1, 2, 2]
          cfg%nyproc = [integer(iintegers) :: Ny]
        case (2) ! stripes along y
          cfg%nxproc = [integer(iintegers) :: Nx]
          cfg%nyproc = [integer(iintegers) :: 1, 2, 2, 1]
        case (3) ! uneven 2x2
          cfg%nxproc = [integer(iintegers) :: 5, 3]
          cfg%nyproc = [integer(iintegers) :: 2, 4]
        end select
      else
        select case (ilayout)
        case (1)
          cfg%nxproc = [integer(iintegers) :: 5, 3]
          cfg%nyproc = [integer(iintegers) :: Ny]
        case (2)
          cfg%nxproc = [integer(iintegers) :: Nx]
          cfg%nyproc = [integer(iintegers) :: 2, 4]
        case (3)
          cfg%nxproc = [integer(iintegers) :: 1, 7]
          cfg%nyproc = [integer(iintegers) :: Ny]
        end select
      end if

      do iphi = 1, size(phis)
        cfg%phi0 = phis(iphi)
        call solve_scene(comm, cfg, res)
        if (allocated(res%edir)) edir(:, :, :, iphi, ilayout) = res%edir
        if (allocated(res%abso)) abso(:, :, :, iphi, ilayout) = res%abso
      end do
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    eps = cfg%edirTOA * 1e-5_ireals * tol_scale
    do iphi = 1, size(phis)
      @assertTrue(all(edir(:, :, :, iphi, 0) .ge. zero), 'missing or negative edir')
      @assertTrue(all(abso(:, :, :, iphi, 0) .ge. zero), 'missing or negative absorption')
      ! the slabs have to cast shadows, i.e. the scene is not trivial
      @assertTrue(minval(edir(Nz + 1, :, :, iphi, 0)) .lt. maxval(edir(Nz + 1, :, :, iphi, 0)) * .5_ireals, 'expected shadows')

      eps_abso = maxval(abso(:, :, :, iphi, 0)) * 1e-4_ireals * tol_scale
      do ilayout = 1, Nlayouts
        msg = 'phi0 '//toStr(phis(iphi))//' layout '//toStr(ilayout)//' on '//toStr(numnodes)//' ranks'
        print *, msg, ' max diff edir', maxval(abs(edir(:, :, :, iphi, ilayout) - edir(:, :, :, iphi, 0))), &
          & 'abso', maxval(abs(abso(:, :, :, iphi, ilayout) - abso(:, :, :, iphi, 0))), 'eps_abso', eps_abso
        @assertEqual(edir(:, :, :, iphi, 0), edir(:, :, :, iphi, ilayout), eps, 'edir depends on the domain decomposition, '//msg)
@assertEqual(abso(:, :, :, iphi, 0), abso(:, :, :, iphi, ilayout), eps_abso, 'absorption depends on the domain decomposition, '//msg)
      end do
    end do
  end subroutine

  ! The inflow at the open boundaries is kept along with the solution.
  ! Reusing a solver for a sequence of different suns, optical properties, buildings and incoming fluxes
  ! must give the same results as a fresh solver each time, i.e. nothing stale is left from the previous solves
  subroutine check_open_bc_solver_reuse_matches_fresh_solver(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d
    real(ireals) :: tol_scale

    integer(iintegers), parameter :: Nsteps = 7

    class(t_solver), allocatable :: solver
    type(t_cfg) :: cfgs(Nsteps)
    type(t_res) :: reused(Nsteps), fresh(Nsteps), again
    real(ireals) :: eps, eps_abso
    integer(iintegers) :: istep
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    cfgs(1)%l2d = l2d
    tol_scale = get_tol_scale(l2d)
    cfgs(1)%phi0 = 20
    cfgs(2) = cfgs(1); cfgs(2)%phi0 = 200       ! sun moves to the opposite quadrant, i.e. the inflow edges change
    cfgs(3) = cfgs(2); cfgs(3)%kabs_scale = 2   ! other optical properties
    cfgs(4) = cfgs(3); cfgs(4)%edirTOA = 2 * cfgs(3)%edirTOA ! brighter sun
    cfgs(5) = cfgs(4); cfgs(5)%lslabs = .false. ! buildings vanish
    cfgs(5)%phi0 = 290; cfgs(5)%theta0 = 30
    cfgs(6) = cfgs(5); cfgs(6)%phi0 = 90        ! only one inflow edge left
    cfgs(7) = cfgs(1)                           ! and back to where we started

    call init_scene(comm, cfgs(1), solver)
    do istep = 1, Nsteps
      call set_angles(solver, spherical_2_cartesian(cfgs(istep)%phi0, cfgs(istep)%theta0))
      if (istep .eq. 1) then
        call run_scene(solver, cfgs(istep), reused(istep), res_again=again)
      else
        call run_scene(solver, cfgs(istep), reused(istep))
      end if
    end do
    call destroy_pprts(solver, lfinalizepetsc=.false.)

    do istep = 1, Nsteps
      call solve_scene(comm, cfgs(istep), fresh(istep))
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do istep = 1, Nsteps
      msg = 'step '//toStr(istep)
      @assertTrue(allocated(reused(istep)%edir), 'missing edir of the reused solver, '//msg)
      @assertTrue(allocated(reused(istep)%abso), 'missing absorption of the reused solver, '//msg)
      @assertTrue(allocated(fresh(istep)%edir), 'missing edir of the fresh solver, '//msg)
      @assertTrue(allocated(fresh(istep)%abso), 'missing absorption of the fresh solver, '//msg)
      @assertTrue(all(ieee_is_finite(reused(istep)%edir)), 'edir is not finite, '//msg)
      @assertTrue(maxval(fresh(istep)%edir(Nz + 1, :, :)) .gt. zero, 'expected light at the surface, '//msg)

      eps = cfgs(istep)%edirTOA * 1e-5_ireals * tol_scale
      eps_abso = maxval(fresh(istep)%abso) * 1e-4_ireals * tol_scale
      print *, msg, ' max diff reused vs fresh solver: edir', maxval(abs(reused(istep)%edir - fresh(istep)%edir)), &
        & 'abso', maxval(abs(reused(istep)%abso - fresh(istep)%abso)), 'eps_abso', eps_abso
      @assertEqual(fresh(istep)%edir, reused(istep)%edir, eps, 'edir of a reused solver differs from a fresh one, '//msg)
      @assertEqual(fresh(istep)%abso, reused(istep)%abso, eps_abso, 'absorption of a reused solver differs from a fresh one, '//msg)
    end do

    ! retrieving the results a second time does not change them, e.g. by scaling the stored fluxes twice
    @assertEqual(reused(1)%edir, again%edir, 'edir changed when retrieving the results a second time')
    @assertEqual(reused(1)%abso, again%abso, 'absorption changed when retrieving the results a second time')

    ! the steps have to differ from each other, otherwise we would not notice stale data
    do istep = 2, Nsteps - 1
@assertTrue(maxval(abs(fresh(istep)%edir - fresh(istep - 1)%edir)) .gt. cfgs(1)%edirTOA * 1e-2_ireals, 'expected the solution to change in step '//toStr(istep))
    end do

    ! the solution is linear in the incoming flux
    eps = cfgs(4)%edirTOA * 1e-5_ireals * tol_scale
    eps_abso = maxval(reused(4)%abso) * 1e-4_ireals * tol_scale
    @assertEqual(2 * reused(3)%edir, reused(4)%edir, eps, 'edir is not linear in edirTOA')
    @assertEqual(2 * reused(3)%abso, reused(4)%abso, eps_abso, 'absorption is not linear in edirTOA')

    ! back at the first configuration
@assertEqual(reused(1)%edir, reused(Nsteps)%edir, cfgs(1)%edirTOA * 1e-5_ireals * tol_scale, 'edir differs after returning to the first configuration')
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_edge_columns_match_replicated_column(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_edge_columns_match_replicated_column(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_edge_columns_match_replicated_column_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_edge_columns_match_replicated_column(this, l2d=.true.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_homogeneous_angles_layers_and_solvers(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_homogeneous_angles_layers_and_solvers(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_homogeneous_angles_layers_and_solvers_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_homogeneous_angles_layers_and_solvers(this, l2d=.true.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_homogeneous_scattering_matches_periodic(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_homogeneous_scattering_matches_periodic(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_homogeneous_scattering_matches_periodic_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_homogeneous_scattering_matches_periodic(this, l2d=.true.)
  end subroutine

  @test(npes=[4, 2])
  subroutine test_open_bc_is_independent_of_domain_decomposition(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_is_independent_of_domain_decomposition(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2])
  subroutine test_open_bc_is_independent_of_domain_decomposition_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_is_independent_of_domain_decomposition(this, l2d=.true.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_solver_reuse_matches_fresh_solver(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_solver_reuse_matches_fresh_solver(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_solver_reuse_matches_fresh_solver_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_solver_reuse_matches_fresh_solver(this, l2d=.true.)
  end subroutine
end module
