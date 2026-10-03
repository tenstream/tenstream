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
    logical :: lopen_x = .true., lopen_y = .true. ! if lopen_bc, which of the boundaries are open
    logical :: lthermal = .false.                ! thermal instead of solar radiation
    integer(iintegers) :: ivar = 1, jvar = 1     ! set to 0 to have no variations of the scene along x or y
    integer(iintegers) :: pad = 0                ! embed the scene in a larger domain, edge columns continue outwards
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
    call set_open_bc_option(.false., .false., .false.)
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
    if (cfg%lheterogeneous) kabs = kabs * (one + real(modulo(3 * i * cfg%ivar + 5 * j * cfg%jvar + k, 4_iintegers), ireals) / 2)
  end function

  ! planck emission [W/m2] at level k (1..Nz+1) of the global column i,j, the surface is warmer than the air above
  pure function scene_planck(cfg, k, i, j) result(planck)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: k, i, j
    real(ireals) :: planck
    planck = 100 + 20 * real(k, ireals)
    if (cfg%lheterogeneous) planck = planck + 15 * real(modulo(i * cfg%ivar + 2 * j * cfg%jvar, 3_iintegers), ireals)
    if (k .gt. Nz + 1) planck = planck + 30
  end function

  ! does the global column i,j carry a slab?
  ! Slabs sit on all four edges, in two of the four corners, and one edge column next to a slab is always free.
  ! If the scene must not vary along x or y, the slab is a bar across the domain
  pure function scene_has_slab(cfg, i, j) result(lslab)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: i, j
    logical :: lslab
    lslab = .false.
    if (.not. cfg%lslabs) return
    if (cfg%ivar .eq. 0 .and. cfg%jvar .eq. 0) return
    if (cfg%ivar .eq. 0) then
      lslab = j .eq. 2
    else if (cfg%jvar .eq. 0) then
      lslab = i .eq. 3
    else
      lslab = (i .eq. 1 .and. j .eq. 1) &
        & .or. (i .eq. cfg%Nx .and. j .eq. cfg%Ny) &
        & .or. (i .eq. 1 .and. j .eq. cfg%Ny - 1) &
        & .or. (i .eq. cfg%Nx .and. j .eq. 2) &
        & .or. (i .eq. 3 .and. j .eq. 1) &
        & .or. (i .eq. cfg%Nx - 2 .and. j .eq. cfg%Ny)
    end if
  end function

  ! column of the scene that provides the properties for the global column i,j of the domain (indices start at 1)
  pure subroutine source_column(cfg, i, j, isrc, jsrc)
    type(t_cfg), intent(in) :: cfg
    integer(iintegers), intent(in) :: i, j
    integer(iintegers), intent(out) :: isrc, jsrc
    if (cfg%icol .gt. 0) then
      isrc = cfg%icol
      jsrc = cfg%jcol
    else ! in the padding, the edge columns continue outwards
      isrc = min(max(i - cfg%pad, 1_iintegers), cfg%Nx)
      jsrc = min(max(j - cfg%pad, 1_iintegers), cfg%Ny)
    end if
  end subroutine

  ! always set all of the options, they stay in the options database
  subroutine set_open_bc_option(lopen_x, lopen_y, l2d)
    logical, intent(in) :: lopen_x, lopen_y, l2d
    integer(mpiint) :: ierr
    call insert_petsc_opt('-pprts_open_bc '//merge('yes', 'no ', lopen_x .and. lopen_y), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_x '//merge('yes', 'no ', lopen_x), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_y '//merge('yes', 'no ', lopen_y), ierr); call CHKERR(ierr)
    call insert_petsc_opt('-pprts_open_bc_2d '//merge('yes', 'no ', l2d), ierr); call CHKERR(ierr)
  end subroutine

  subroutine init_scene(comm, cfg, solver)
    integer(mpiint), intent(in) :: comm
    type(t_cfg), intent(in) :: cfg
    class(t_solver), allocatable, intent(inout) :: solver
    logical :: lx, ly
    integer(mpiint) :: ierr

    lx = cfg%lopen_bc .and. cfg%lopen_x
    ly = cfg%lopen_bc .and. cfg%lopen_y
    call set_open_bc_option(lx, ly, cfg%l2d)
    call allocate_pprts_solver_from_commandline(solver, trim(cfg%solvername), ierr); call CHKERR(ierr)

    if (allocated(cfg%nxproc)) then
      call init_pprts(comm, Nz, cfg%Nx, cfg%Ny, dx, dx, spherical_2_cartesian(cfg%phi0, cfg%theta0), solver, &
        & dz1d=cfg%dz, nxproc=cfg%nxproc, nyproc=cfg%nyproc)
    else
      call init_pprts(comm, Nz, cfg%Nx + 2 * cfg%pad, cfg%Ny + 2 * cfg%pad, dx, dx, &
        & spherical_2_cartesian(cfg%phi0, cfg%theta0), solver, dz1d=cfg%dz)
    end if
    if ((lx .or. ly) .neqv. solver%lopen_bc) loption_mismatch = .true.
    if (lx .neqv. solver%lopen_bc_x) loption_mismatch = .true.
    if (ly .neqv. solver%lopen_bc_y) loption_mismatch = .true.
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
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, planck
    real(ireals), allocatable, dimension(:, :) :: planck_srfc
    integer(iintegers) :: Nfaces, m, i, j, k, isrc, jsrc, faceid
    logical :: lbuildings
    integer(mpiint) :: ierr

    ! if there are slabs anywhere in the domain, has to be the same on all ranks
    lbuildings = cfg%lslabs .and. .not. cfg%lthermal .and. (cfg%ivar .ne. 0 .or. cfg%jvar .ne. 0)
    if (cfg%icol .gt. 0) lbuildings = scene_has_slab(cfg, cfg%icol, cfg%jcol)

    associate (C => solver%C_one)
      allocate (kabs(C%zm, C%xm, C%ym))
      allocate (ksca(C%zm, C%xm, C%ym))
      allocate (g(C%zm, C%xm, C%ym), source=cfg%g)
      if (cfg%lthermal) allocate (planck(C%zm + 1, C%xm, C%ym), planck_srfc(C%xm, C%ym))

      Nfaces = 0
      do j = C%ys, C%ye
        do i = C%xs, C%xe
          call source_column(cfg, i + 1, j + 1, isrc, jsrc)
          do k = 1, Nz
            kabs(k, i - C%xs + 1, j - C%ys + 1) = scene_kabs(cfg, k, isrc, jsrc)
          end do
          if (cfg%lthermal) then
            do k = 1, Nz + 1
              planck(k, i - C%xs + 1, j - C%ys + 1) = scene_planck(cfg, k, isrc, jsrc)
            end do
            planck_srfc(i - C%xs + 1, j - C%ys + 1) = scene_planck(cfg, Nz + 2, isrc, jsrc)
          end if
          if (lbuildings .and. scene_has_slab(cfg, isrc, jsrc)) Nfaces = Nfaces + 6 * (kslab_bot - kslab_top + 1)
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

    if (cfg%lthermal) then
      call set_optical_properties(solver, cfg%albedo, kabs, ksca, g, planck, planck_srfc)
      call solve_pprts(solver, lthermal=.true., lsolar=.false., edirTOA=cfg%edirTOA)
      call pprts_get_result_toZero(solver, res%edn, res%eup, res%abso)
      if (present(res_again)) call pprts_get_result_toZero(solver, res_again%edn, res_again%eup, res_again%abso)
    else
      call set_optical_properties(solver, cfg%albedo, kabs, ksca, g)
      if (lbuildings) then
        call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=cfg%edirTOA, opt_buildings=buildings)
      else
        call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=cfg%edirTOA)
      end if
      call pprts_get_result_toZero(solver, res%edn, res%eup, res%abso, res%edir)
      if (present(res_again)) call pprts_get_result_toZero(solver, res_again%edn, res_again%eup, res_again%abso, res_again%edir)
    end if

    if (cfg%pad .gt. 0) then ! cut out the scene
      call cut(res%edn); call cut(res%eup); call cut(res%abso); call cut(res%edir)
    end if

    allocate (res%l1d(Nz))
    res%l1d = solver%atm%l1d
    res%Ndir_streams = solver%dirtop%dof + solver%dirside%dof * 2

    if (lbuildings) then
      call destroy_buildings(buildings, ierr); call CHKERR(ierr)
    end if
  contains
    subroutine cut(arr)
      real(ireals), allocatable, intent(inout) :: arr(:, :, :)
      real(ireals), allocatable :: tmp(:, :, :)
      if (.not. allocated(arr)) return
      allocate (tmp, source=arr(:, cfg%pad + 1:cfg%pad + cfg%Nx, cfg%pad + 1:cfg%pad + cfg%Ny))
      call move_alloc(tmp, arr)
    end subroutine
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
    real(ireals), dimension(Nz + 1, Nx, Ny, size(phis), 0:Nlayouts) :: edir, edn, eup
    real(ireals), dimension(Nz, Nx, Ny, size(phis), 0:Nlayouts) :: abso
    real(ireals) :: eps, eps_abso
    ! the diffuse solver iterates differently for each decomposition, results agree up to its convergence criteria
    real(ireals), parameter :: eps_diff = 1e-1_ireals
    integer(iintegers) :: iphi, ilayout
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid, numnodes

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()
    numnodes = this%getNumProcesses()
    cfg%l2d = l2d
    tol_scale = get_tol_scale(l2d)

    edir = -one
    edn = -one
    eup = -one
    abso = -one

    ! with scattering and a reflecting surface, to also cover the diffuse open boundaries
    cfg%w0 = .5_ireals
    cfg%g = .5_ireals
    cfg%albedo = .2_ireals

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
        if (allocated(res%edn)) edn(:, :, :, iphi, ilayout) = res%edn
        if (allocated(res%eup)) eup(:, :, :, iphi, ilayout) = res%eup
        if (allocated(res%abso)) abso(:, :, :, iphi, ilayout) = res%abso
      end do
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    eps = cfg%edirTOA * 1e-5_ireals * tol_scale
    do iphi = 1, size(phis)
      @assertTrue(all(edir(:, :, :, iphi, 0) .ge. zero), 'missing or negative edir')
      @assertTrue(all(edn(:, :, :, iphi, 0) .ge. zero), 'missing or negative edn')
      @assertTrue(maxval(eup(:, :, :, iphi, 0)) .gt. 10._ireals, 'expected upwelling diffuse radiation')
      @assertTrue(all(abso(:, :, :, iphi, 0) .ge. zero), 'missing or negative absorption')
      ! the slabs have to cast shadows, i.e. the scene is not trivial
      @assertTrue(minval(edir(Nz + 1, :, :, iphi, 0)) .lt. maxval(edir(Nz + 1, :, :, iphi, 0)) * .5_ireals, 'expected shadows')

      eps_abso = maxval(abso(:, :, :, iphi, 0)) * 1e-3_ireals
      do ilayout = 1, Nlayouts
        msg = 'phi0 '//toStr(phis(iphi))//' layout '//toStr(ilayout)//' on '//toStr(numnodes)//' ranks'
        print *, msg, ' max diff edir', maxval(abs(edir(:, :, :, iphi, ilayout) - edir(:, :, :, iphi, 0))), &
          & 'abso', maxval(abs(abso(:, :, :, iphi, ilayout) - abso(:, :, :, iphi, 0))), 'eps_abso', eps_abso
        @assertEqual(edir(:, :, :, iphi, 0), edir(:, :, :, iphi, ilayout), eps, 'edir depends on the domain decomposition, '//msg)
        print *, msg, ' max diff edn', maxval(abs(edn(:, :, :, iphi, ilayout) - edn(:, :, :, iphi, 0))), &
          & 'eup', maxval(abs(eup(:, :, :, iphi, ilayout) - eup(:, :, :, iphi, 0)))
        @assertEqual(edn(:, :, :, iphi, 0), edn(:, :, :, iphi, ilayout), eps_diff, 'edn depends on the domain decomposition, '//msg)
        @assertEqual(eup(:, :, :, iphi, 0), eup(:, :, :, iphi, ilayout), eps_diff, 'eup depends on the domain decomposition, '//msg)
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

  ! Open boundaries for diffuse radiation: what enters an edge cell through the domain edge is what leaves it in that direction.
  ! If the scene does not vary along x, neither does the radiation field and periodic boundaries in x are exact.
  ! Opening the x boundaries then must not change anything, no matter how the scene looks like along y. Same for y.
  ! This holds for direct, diffuse and thermal radiation and for the absorption
  subroutine check_open_bc_along_invariant_direction_matches_periodic(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d

    real(ireals), parameter :: phis(6) = [real(ireals) :: 20, 110, 200, 290, 0, 90]
    integer(iintegers), parameter :: Ncases = 2 * (size(phis) + 1) ! per direction: solar for each sun azimuth, plus thermal

    type(t_cfg) :: cfg, cfgs(Ncases)
    type(t_res) :: open (Ncases), periodic(Ncases)
    real(ireals) :: eps, eps_abso, variation
    integer(iintegers) :: icase, idir, iphi
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    cfg%l2d = l2d
    cfg%kabs_scale = 3
    cfg%w0 = .7_ireals
    cfg%g = .5_ireals
    cfg%albedo = .3_ireals

    icase = 0
    do idir = 1, 2
      do iphi = 1, size(phis) + 1
        icase = icase + 1
        cfgs(icase) = cfg
        if (idir .eq. 1) then ! no variations along x, open in x
          cfgs(icase)%ivar = 0
          cfgs(icase)%lopen_y = .false.
        else
          cfgs(icase)%jvar = 0
          cfgs(icase)%lopen_x = .false.
        end if
        if (iphi .le. size(phis)) then
          cfgs(icase)%phi0 = phis(iphi)
        else
          cfgs(icase)%lthermal = .true.
        end if
      end do
    end do

    do icase = 1, Ncases
      cfg = cfgs(icase)
      call solve_scene(comm, cfg, open (icase))
      cfg%lopen_bc = .false.
      call solve_scene(comm, cfg, periodic(icase))
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do icase = 1, Ncases
      associate (c => cfgs(icase), o => open (icase), p => periodic(icase))
        msg = 'open '//merge('x', 'y', c%lopen_x)//merge(' thermal', ' solar  ', c%lthermal)//' phi0 '//toStr(c%phi0)

        ! the direct inflow from single column solves does not know about the neighbours along the edge,
        ! i.e. it only gives the invariant solution if the sun is aligned with the grid
        if (.not. l2d .and. .not. c%lthermal .and. modulo(nint(c%phi0), 90) .ne. 0) cycle

        @assertTrue(allocated(o%edn) .and. allocated(o%eup) .and. allocated(o%abso), 'missing results, '//msg)
        @assertTrue(allocated(p%edn) .and. allocated(p%eup) .and. allocated(p%abso), 'missing results, '//msg)
     @assertTrue(all(ieee_is_finite(o%edn)) .and. all(ieee_is_finite(o%eup)) .and. all(ieee_is_finite(o%abso)), 'not finite, '//msg)

        ! the diffuse radiation field has to be substantial and has to vary along the edge, otherwise this test is void
        @assertTrue(minval(p%eup(1, :, :)) .gt. 10._ireals, 'expected upwelling diffuse radiation at the top of the domain, '//msg)
        if (c%lopen_x) then
          variation = maxval(maxval(p%eup(Nz + 1, :, :), dim=1) - minval(p%eup(Nz + 1, :, :), dim=1))
          @assertEqual(zero, variation, 1e-1_ireals, 'expected no variations along x in the periodic solution, '//msg)
          variation = maxval(maxval(p%eup(:, 1, :), dim=2) - minval(p%eup(:, 1, :), dim=2))
        else
          variation = maxval(maxval(p%eup(Nz + 1, :, :), dim=2) - minval(p%eup(Nz + 1, :, :), dim=2))
          @assertEqual(zero, variation, 1e-1_ireals, 'expected no variations along y in the periodic solution, '//msg)
          variation = maxval(maxval(p%eup(:, :, 1), dim=2) - minval(p%eup(:, :, 1), dim=2))
        end if
        @assertTrue(variation .gt. one, 'expected the upwelling radiation to vary along the open edge, '//msg)

        eps = 1e-1_ireals
        eps_abso = maxval(abs(p%abso)) * 1e-3_ireals
        print *, msg, ' max diff open vs periodic: edn', maxval(abs(o%edn - p%edn)), 'eup', maxval(abs(o%eup - p%eup)), &
          & 'abso', maxval(abs(o%abso - p%abso)), 'eps_abso', eps_abso, 'variation along the edge', variation

        if (.not. c%lthermal) then
          @assertEqual(p%edir, o%edir, eps, 'edir changed by opening the boundaries along the invariant direction, '//msg)
        end if
        @assertEqual(p%edn, o%edn, eps, 'edn changed by opening the boundaries along the invariant direction, '//msg)
        @assertEqual(p%eup, o%eup, eps, 'eup changed by opening the boundaries along the invariant direction, '//msg)
        @assertEqual(p%abso, o%abso, eps_abso, 'absorption changed by opening the boundaries along the invariant direction, '//msg)
      end associate
    end do
  end subroutine

  ! Diffuse radiation travels in all directions. Other than for direct radiation, the open boundaries are therefore
  ! only an approximation to a scene whose edge columns continue outwards forever: with -pprts_open_bc_2d a layer of
  ! ghost cells continues the edge columns and beyond them we assume zero gradient, otherwise zero gradient right at the edge.
  ! It has to be much closer to that than periodic boundaries though. Pin the quality of the approximation
  subroutine check_open_bc_diffuse_is_close_to_embedded_domain(this, l2d)
    class(MpiTestMethod), intent(inout) :: this
    logical, intent(in) :: l2d ! use -pprts_open_bc_2d

    real(ireals), parameter :: phis(4) = [real(ireals) :: 20, 110, 200, 290]
    integer(iintegers), parameter :: Ncases = size(phis) + 1 ! solar for each sun azimuth, plus thermal
    integer(iintegers), parameter :: pad = 8

    type(t_cfg) :: cfg, cfgs(Ncases)
    type(t_res) :: open (Ncases), periodic(Ncases), embedded(Ncases)
    real(ireals) :: err_open, err_periodic, rmse_open, rmse_periodic
    integer(iintegers) :: icase
    character(len=:), allocatable :: msg
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    cfg%l2d = l2d
    cfg%kabs_scale = 3
    cfg%w0 = .7_ireals
    cfg%g = .5_ireals
    cfg%albedo = .3_ireals

    do icase = 1, Ncases
      cfgs(icase) = cfg
      if (icase .le. size(phis)) then
        cfgs(icase)%phi0 = phis(icase)
      else
        cfgs(icase)%lthermal = .true.
      end if
    end do

    do icase = 1, Ncases
      cfg = cfgs(icase)
      call solve_scene(comm, cfg, open (icase))
      cfg%lopen_bc = .false.
      call solve_scene(comm, cfg, periodic(icase))
      cfg%pad = pad
      call solve_scene(comm, cfg, embedded(icase))
    end do

    @assertFalse(loption_mismatch, 'solver did not pick up the -pprts_open_bc options')
    if (myid .ne. 0) return

    do icase = 1, Ncases
      associate (c => cfgs(icase), o => open (icase), p => periodic(icase), e => embedded(icase))
        msg = merge('thermal', 'solar  ', c%lthermal)//' phi0 '//toStr(c%phi0)
        @assertTrue(allocated(o%edn) .and. allocated(p%edn) .and. allocated(e%edn), 'missing results, '//msg)
        @assertTrue(all(shape(o%edn) .eq. shape(e%edn)), 'embedded result has the wrong shape, '//msg)
     @assertTrue(all(ieee_is_finite(o%edn)) .and. all(ieee_is_finite(o%eup)) .and. all(ieee_is_finite(o%abso)), 'not finite, '//msg)

        err_open = max(maxval(abs(o%edn - e%edn)), maxval(abs(o%eup - e%eup)))
        err_periodic = max(maxval(abs(p%edn - e%edn)), maxval(abs(p%eup - e%eup)))
        rmse_open = sqrt((sum((o%edn - e%edn)**2) + sum((o%eup - e%eup)**2)) / real(2 * size(e%edn), ireals))
        rmse_periodic = sqrt((sum((p%edn - e%edn)**2) + sum((p%eup - e%eup)**2)) / real(2 * size(e%edn), ireals))
        print *, msg, ' diffuse fluxes vs embedded: max err open', err_open, 'periodic', err_periodic, &
          & 'rmse open', rmse_open, 'periodic', rmse_periodic, 'mean eup', sum(e%eup) / real(size(e%eup), ireals)

        @assertTrue(rmse_periodic .gt. one, 'expected periodic boundaries to differ from the embedded domain, '//msg)
@assertTrue(rmse_open .lt. rmse_periodic * merge(.2_ireals, .5_ireals, l2d), 'open bc diffuse fluxes are not closer to the embedded domain than periodic ones, '//msg)
@assertTrue(err_open .lt. err_periodic, 'open bc diffuse fluxes locally differ more from the embedded domain than periodic ones, '//msg)

        ! the absorption of the edge cells needs the fluxes through the domain edges
        rmse_open = sqrt(sum((o%abso - e%abso)**2) / real(size(e%abso), ireals))
        rmse_periodic = sqrt(sum((p%abso - e%abso)**2) / real(size(e%abso), ireals))
        print *, msg, ' absorption vs embedded: rmse open', rmse_open, 'periodic', rmse_periodic, &
          & 'max err open', maxval(abs(o%abso - e%abso)), 'periodic', maxval(abs(p%abso - e%abso))
        ! with the inflow from single column solves, the absorption of direct radiation is further off
@assertTrue(rmse_open .lt. rmse_periodic * merge(.25_ireals, .75_ireals, l2d), 'open bc absorption is not closer to the embedded domain than periodic one, '//msg)
@assertTrue(maxval(abs(o%abso - e%abso)) .lt. maxval(abs(p%abso - e%abso)), 'open bc absorption locally differs more from the embedded domain than periodic one, '//msg)
      end associate
    end do
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

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_along_invariant_direction_matches_periodic(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_along_invariant_direction_matches_periodic(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_along_invariant_direction_matches_periodic_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_along_invariant_direction_matches_periodic(this, l2d=.true.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_diffuse_is_close_to_embedded_domain(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_diffuse_is_close_to_embedded_domain(this, l2d=.false.)
  end subroutine

  @test(npes=[4, 2, 1])
  subroutine test_open_bc_diffuse_is_close_to_embedded_domain_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    call check_open_bc_diffuse_is_close_to_embedded_domain(this, l2d=.true.)
  end subroutine
end module
