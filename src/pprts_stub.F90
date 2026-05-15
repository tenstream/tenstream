! Stub module providing m_pprts public interface when PETSc is not available.
! Implements init_pprts, set_angles, and the result/gather calls using plain arrays.
! PETSc-dependent calls (solve_pprts, set_optical_properties, etc.) still abort.
module m_pprts
  use mpi, only: MPI_COMM_RANK, MPI_COMM_SIZE
  use m_data_parameters, only: &
    & iintegers, ireals, irealLUT, mpiint, &
    & init_mpi_data_parameters, &
    & zero, one, &
    & i0, i1, i2
  use m_helper_functions, only: &
    & approx, &
    & cartesian_2_spherical, &
    & CHKERR, &
    & CHKWARN, &
    & deg2rad, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_max, &
    & imp_allreduce_sum, &
    & rad2deg, &
    & spherical_2_cartesian, &
    & toStr
  use m_pprts_base, only: &
    & atmk, &
    & setup_coord_native, &
    & t_atmosphere, &
    & t_coord, &
    & t_dof, &
    & t_solver, &
    & t_solver_1_2, &
    & t_solver_2str, &
    & t_solver_3_10, &
    & t_solver_3_16, &
    & t_solver_3_24, &
    & t_solver_3_30, &
    & t_solver_3_6, &
    & t_solver_8_10, &
    & t_solver_8_16, &
    & t_solver_8_18, &
    & t_solver_disort, &
    & t_solver_mcdmda, &
    & t_solver_rayli, &
    & t_state_container, &
    & t_suninfo
  use m_optprop, only: &
    & t_optprop_1_2, &
    & t_optprop_3_6, &
    & t_optprop_3_10, &
    & t_optprop_3_10_ann, &
    & t_optprop_3_16, &
    & t_optprop_3_24, &
    & t_optprop_3_30, &
    & t_optprop_8_10, &
    & t_optprop_8_16, &
    & t_optprop_8_18
  use m_tenstream_options, only: &
    & luse_eddington, &
    & options_max_solution_err, &
    & options_max_solution_time, &
    & read_commandline_options, &
    & twostr_ratio
  use m_buildings, only: t_pprts_buildings
  implicit none
  private

  public :: &
    & init_pprts, &
    & set_optical_properties, set_global_optical_properties, &
    & solve_pprts, set_angles, pprts_get_result, &
    & pprts_get_result_toZero, &
    & gather_all_toZero, &
    & gather_all_to_all, &
    & scale_flx

  integer(iintegers), parameter :: minimal_dimension = 3

contains

  subroutine init_pprts(icomm, Nz, Nx, Ny, dx, dy, sundir, solver, &
      & dz1d, dz3d, nxproc, nyproc, collapseindex, solvername)
    integer(mpiint), intent(in) :: icomm
    integer(iintegers), intent(in) :: Nz, Nx, Ny
    real(ireals), intent(in) :: dx, dy, sundir(:)
    class(t_solver), intent(inout) :: solver
    real(ireals), optional, intent(in) :: dz1d(:), dz3d(:, :, :)
    integer(iintegers), optional, intent(in) :: nxproc(:), nyproc(:), collapseindex
    character(len=*), optional, intent(in) :: solvername

    integer(iintegers) :: Nz_in, Nz_solver, Nx_eff, Ny_eff, k, i, j
    logical :: luse_ann, lflg
    integer(mpiint) :: ierr

    if (solver%linitialized) then
      print *, 'init_pprts: solver already initialized'
      ierr = 1; call CHKERR(ierr)
      return
    end if

    call init_mpi_data_parameters(icomm)
    solver%comm = icomm
    call MPI_COMM_RANK(solver%comm, solver%myid, ierr); call CHKERR(ierr)
    call MPI_COMM_SIZE(solver%comm, solver%numnodes, ierr); call CHKERR(ierr)

    call read_commandline_options(solver%comm)

    luse_ann = .false.
    call get_petsc_opt(solver%prefix, '-pprts_use_ANN', luse_ann, lflg, ierr); call CHKERR(ierr)

    solver%difftop%area_divider = 1
    solver%diffside%area_divider = 1
    solver%dirtop%area_divider = 1
    solver%dirside%area_divider = 1

    select type (solver)
    class is (t_solver_2str)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(0))
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(0))

    class is (t_solver_disort)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(0))
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(0))

    class is (t_solver_mcdmda)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(4))
      solver%diffside%is_inward = [.false., .true.]
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = [.true.]

    class is (t_solver_rayli)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(4))
      solver%diffside%is_inward = [.false., .true., .false., .true.]
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = [.true.]

    class is (t_solver_1_2)
      allocate (t_optprop_1_2 :: solver%OPP)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(0))
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(0))

    class is (t_solver_3_6)
      allocate (t_optprop_3_6 :: solver%OPP)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(2))
      solver%diffside%is_inward = [.false., .true.]
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = [.true.]

    class is (t_solver_3_10)
      if (luse_ann) then
        allocate (t_optprop_3_10_ann :: solver%OPP)
      else
        allocate (t_optprop_3_10 :: solver%OPP)
      end if
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(4))
      solver%diffside%is_inward = [.false., .true., .false., .true.]
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = [.true.]
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = [.true.]

    class is (t_solver_8_10)
      allocate (t_optprop_8_10 :: solver%OPP)
      allocate (solver%difftop%is_inward(2))
      solver%difftop%is_inward = [.false., .true.]
      allocate (solver%diffside%is_inward(4))
      solver%diffside%is_inward = [.false., .true., .false., .true.]
      allocate (solver%dirtop%is_inward(4))
      solver%dirtop%is_inward = .true.
      solver%dirtop%area_divider = 4
      allocate (solver%dirside%is_inward(2))
      solver%dirside%is_inward = .true.
      solver%dirside%area_divider = 2

    class is (t_solver_3_16)
      allocate (t_optprop_3_16 :: solver%OPP)
      allocate (solver%difftop%is_inward(8), source= &
                [.false., .true., .false., .true., .false., .true., .false., .true.])
      allocate (solver%diffside%is_inward(4))
      solver%diffside%is_inward = [.false., .true., .false., .true.]
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = .true.
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = .true.

    class is (t_solver_3_24)
      allocate (t_optprop_3_24 :: solver%OPP)
      allocate (solver%difftop%is_inward(8), source= &
                [.false., .true., .false., .true., .false., .true., .false., .true.])
      allocate (solver%diffside%is_inward(8), source= &
                [.true., .false., .true., .false., .true., .false., .true., .false.])
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = .true.
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = .true.

    class is (t_solver_3_30)
      allocate (t_optprop_3_30 :: solver%OPP)
      allocate (solver%difftop%is_inward(10), source= &
                [.false., .true., .false., .true., .false., .true., .false., .true., .false., .true.])
      allocate (solver%diffside%is_inward(10), source= &
                [.true., .false., .true., .false., .true., .false., .true., .false., .true., .false.])
      allocate (solver%dirtop%is_inward(1))
      solver%dirtop%is_inward = .true.
      allocate (solver%dirside%is_inward(1))
      solver%dirside%is_inward = .true.

    class is (t_solver_8_16)
      allocate (t_optprop_8_16 :: solver%OPP)
      allocate (solver%difftop%is_inward(8), source= &
                [.false., .true., .false., .true., .false., .true., .false., .true.])
      allocate (solver%diffside%is_inward(4), source=[.false., .true., .false., .true.])
      allocate (solver%dirtop%is_inward(4), source=.true.)
      solver%dirtop%area_divider = 4
      allocate (solver%dirside%is_inward(2), source=.true.)
      solver%dirside%area_divider = 2

    class is (t_solver_8_18)
      allocate (t_optprop_8_18 :: solver%OPP)
      allocate (solver%difftop%is_inward(10), source= &
                [.false., .true., .false., .true., .false., .true., .false., .true., .false., .true.])
      allocate (solver%diffside%is_inward(4), source=[.false., .true., .false., .true.])
      allocate (solver%dirtop%is_inward(4), source=.true.)
      solver%dirtop%area_divider = 4
      allocate (solver%dirside%is_inward(2), source=.true.)
      solver%dirside%area_divider = 2

    class default
      call CHKERR(1_mpiint, 'init_pprts: unexpected solver type')
    end select

    solver%difftop%dof = size(solver%difftop%is_inward)
    solver%diffside%dof = size(solver%diffside%is_inward)
    solver%dirtop%dof = size(solver%dirtop%is_inward)
    solver%dirside%dof = size(solver%dirside%is_inward)

    solver%difftop%streams = solver%difftop%dof / 2
    solver%diffside%streams = solver%diffside%dof / 2
    solver%dirtop%streams = solver%dirtop%dof
    solver%dirside%streams = solver%dirside%dof

    allocate (solver%solutions(-1000:1000))
    solver%lenable_solutions_err_estimates = &
      options_max_solution_err .gt. zero .and. options_max_solution_time .gt. zero

    if (.not. approx(dx, dy)) &
      call CHKERR(1_mpiint, 'dx and dy currently have to be the same '//toStr(dx)//' vs '//toStr(dy))

    call get_petsc_opt('', '-pprts_open_bc', solver%lopen_bc, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, '-pprts_open_bc', solver%lopen_bc, lflg, ierr); call CHKERR(ierr)

    ! Atmosphere layers count
    Nz_in = Nz

    ! Solver grid may be smaller when collapsing upper layers
    Nz_solver = Nz_in
    if (present(collapseindex)) then
      if (collapseindex .gt. i1) Nz_solver = Nz_in - collapseindex + i1
    end if

    ! Without explicit decomposition, enforce minimum grid size
    if (present(nxproc) .and. present(nyproc)) then
      Nx_eff = Nx
      Ny_eff = Ny
    else
      Nx_eff = max(minimal_dimension, Nx)
      Ny_eff = max(minimal_dimension, Ny)
    end if

    ! Set up all coordinate objects using native MPI Cartesian decomposition
    call setup_coord_native(icomm, Nz_solver + 1, Nx_eff, Ny_eff, &
      & solver%difftop%dof + 2 * solver%diffside%dof, solver%C_diff, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_solver + 1, Nx_eff, Ny_eff, &
      & solver%dirtop%dof + 2 * solver%dirside%dof, solver%C_dir, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_solver, Nx_eff, Ny_eff, &
      & i1, solver%C_one, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_solver + 1, Nx_eff, Ny_eff, &
      & i1, solver%C_one1, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_solver + 1, Nx_eff, Ny_eff, &
      & i2, solver%C_two1, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_in, Nx_eff, Ny_eff, &
      & i1, solver%C_one_atm, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_in + 1, Nx_eff, Ny_eff, &
      & i1, solver%C_one_atm1, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_in + 1, Nx_eff, Ny_eff, &
      & i1, solver%C_one_atm1_box, nxproc, nyproc)

    call setup_coord_native(icomm, i1, Nx_eff, Ny_eff, &
      & i1, solver%Csrfc_one, nxproc, nyproc)

    call setup_coord_native(icomm, Nz_in + 1, Nx_eff + 1, Ny_eff + 1, &
      & i1, solver%Cvert_one_atm1)

    call setup_atm()

    if (present(solvername)) then
      solver%solvername = trim(solver%solvername)//trim(solvername)
    else
      solver%solvername = ''
    end if

    solver%linitialized = .true.

    call set_angles(solver, sundir)

    select type (solver)
    class is (t_solver_2str)
      return
    class is (t_solver_disort)
      return
    class is (t_solver_rayli)
      return
    class is (t_solver_mcdmda)
      return
    end select

    if (.not. luse_eddington) then
      allocate (t_optprop_1_2 :: solver%OPP1d)
      call solver%OPP1d%init(solver%comm)
    end if

    call solver%OPP%init(solver%comm)

  contains
    subroutine setup_atm()
      real(ireals) :: N1dlayers, N1dlayers_max
      integer(iintegers) :: count3d, countall

      if (.not. allocated(solver%atm)) allocate (solver%atm)
      solver%atm%dx = dx
      solver%atm%dy = dy

      if (.not. allocated(solver%atm%dz)) then
        allocate (solver%atm%dz( &
          & solver%C_one_atm%zs:solver%C_one_atm%ze, &
          & solver%C_one_atm%xs:solver%C_one_atm%xe, &
          & solver%C_one_atm%ys:solver%C_one_atm%ye))
      end if

      if (present(dz1d)) then
        do j = solver%C_one_atm%ys, solver%C_one_atm%ye
          do i = solver%C_one_atm%xs, solver%C_one_atm%xe
            solver%atm%dz(:, i, j) = dz1d
          end do
        end do
      end if
      if (present(dz3d)) then
        if (any(shape(dz3d) .ne. shape(solver%atm%dz))) then
          call CHKERR(1_mpiint, 'dz3d shape mismatch: '// &
            & toStr(shape(dz3d))//' vs '//toStr(shape(solver%atm%dz)))
        end if
        solver%atm%dz = dz3d
      end if
      if (.not. present(dz1d) .and. .not. present(dz3d)) then
        call CHKERR(1_mpiint, 'init_pprts: must supply dz1d or dz3d')
      end if

      solver%sun%luse_topography = .false. ! topography requires PETSc Vecs

      ! Determine 1D layers (no -pprts_1d_height support without PETSc)
      if (.not. allocated(solver%atm%l1d)) then
        allocate (solver%atm%l1d(solver%C_one_atm%zs:solver%C_one_atm%ze))
      end if
      solver%atm%l1d = .false.

      select type (solver)
      class is (t_solver_2str)
        solver%atm%l1d = .true.
      class is (t_solver_disort)
        solver%atm%l1d = .true.
      class default
        associate (C => solver%C_one_atm, dz => solver%atm%dz)
          solver%atm%l1d(C%ze) = any((dz(C%ze, C%xs:C%xe, C%ys:C%ye) / solver%atm%dx) .gt. twostr_ratio)
          do k = C%ze - 1, C%zs, -1
            if (any((dz(k, C%xs:C%xe, C%ys:C%ye) / solver%atm%dx) .gt. twostr_ratio)) then
              solver%atm%l1d(C%zs:k) = .true.
              exit
            end if
          end do
        end associate
      end select

      if (present(collapseindex)) then
        solver%atm%lcollapse = collapseindex .gt. i1
        solver%atm%icollapse = collapseindex
        if (solver%atm%lcollapse) then
          ierr = count(.not. solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs)))
          if (ierr .ne. 0_mpiint) then
            call CHKWARN(ierr, 'Non-1D layers in collapse region -- collapse index: '//toStr(collapseindex))
          end if
          solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs)) = .true.
        end if
      end if

      N1dlayers = count(solver%atm%l1d)
      call imp_allreduce_max(solver%comm, N1dlayers, N1dlayers_max)
      if (N1dlayers_max .gt. 0) then
        ierr = int(N1dlayers - N1dlayers_max, mpiint)
        if (ierr .ne. 0_mpiint) then
          call CHKWARN(ierr, 'Nr of 1D layers does not match across ranks')
        end if
        solver%atm%l1d(solver%C_one_atm%zs:solver%C_one_atm%zs + int(N1dlayers_max, iintegers) - i1) = .true.
      end if

      call imp_allreduce_sum(solver%comm, count(.not. solver%atm%l1d, kind=iintegers), count3d)
      call imp_allreduce_sum(solver%comm, size(solver%atm%l1d, kind=iintegers), countall)
      solver%atm%unconstrained_fraction = real(count3d, ireals) / real(countall, ireals)

    end subroutine setup_atm
  end subroutine init_pprts

  subroutine set_angles(solver, sundir)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: sundir(:)

    real(ireals) :: round_sun_angles, round_sun_theta, round_sun_phi
    logical :: lflg
    integer(mpiint) :: ierr

    if (.not. solver%linitialized) then
      print *, 'set_angles: solver not yet initialized'
      ierr = 1; call CHKERR(ierr)
    end if

    call cartesian_2_spherical(sundir, solver%sun%phi, solver%sun%theta, ierr)
    solver%sun%mu = max(cos(deg2rad(solver%sun%theta)), zero)

    round_sun_angles = -1._ireals
    call get_petsc_opt('', '-pprts_round_sun_angles', round_sun_angles, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, '-pprts_round_sun_angles', round_sun_angles, lflg, ierr); call CHKERR(ierr)
    round_sun_theta = round_sun_angles
    round_sun_phi = round_sun_angles
    call get_petsc_opt('', '-pprts_round_sun_theta', round_sun_theta, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, '-pprts_round_sun_theta', round_sun_theta, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-pprts_round_sun_phi', round_sun_phi, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, '-pprts_round_sun_phi', round_sun_phi, lflg, ierr); call CHKERR(ierr)

    if (round_sun_theta .gt. 0._ireals) &
      solver%sun%theta = nint(solver%sun%theta / round_sun_theta) * round_sun_theta
    if (round_sun_phi .gt. 0._ireals) &
      solver%sun%phi = nint(solver%sun%phi / round_sun_phi) * round_sun_phi

    call get_petsc_opt(solver%prefix, '-pprts_force_zenith', solver%sun%theta, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, '-pprts_force_azimuth', solver%sun%phi, lflg, ierr); call CHKERR(ierr)

    solver%sun%sundir = spherical_2_cartesian(solver%sun%phi, solver%sun%theta)

    if (solver%sun%theta .ge. 90._ireals) solver%sun%theta = -one

    solver%sun%costheta = max(cos(deg2rad(solver%sun%theta)), zero)

    ! Reduce phi to symmetry range [0,90]
    solver%sun%symmetry_phi = sym_rot_phi(solver%sun%phi)

    if (sin(deg2rad(solver%sun%phi)) .gt. zero) then
      solver%sun%xinc = i0
    else
      solver%sun%xinc = i1
    end if

    if (cos(deg2rad(solver%sun%phi)) .lt. zero) then
      solver%sun%yinc = i1
    else
      solver%sun%yinc = i0
    end if

  contains
    pure elemental function sym_rot_phi(phi)
      real(ireals) :: sym_rot_phi
      real(ireals), intent(in) :: phi
      sym_rot_phi = acos(cos(deg2rad(phi)))
      sym_rot_phi = min(90._ireals, max(0._ireals, rad2deg(asin(sin(sym_rot_phi)))))
    end function
  end subroutine set_angles

  subroutine set_optical_properties(solver, albedo, kabs, ksca, g, &
      & planck, planck_srfc, albedo_2d, ldelta_scaling)
    class(t_solver) :: solver
    real(ireals), intent(in) :: albedo
    real(ireals), intent(in), dimension(:, :, :), optional :: kabs, ksca, g, planck
    real(ireals), intent(in), dimension(:, :), optional :: planck_srfc, albedo_2d
    logical, intent(in), optional :: ldelta_scaling
    call CHKERR(1_mpiint, 'set_optical_properties requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine set_global_optical_properties(solver, global_albedo, global_kabs, global_ksca, global_g, global_planck)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(inout), optional :: global_albedo
    real(ireals), intent(inout), dimension(:, :, :), allocatable, optional :: global_kabs, global_ksca, global_g, global_planck
    call CHKERR(1_mpiint, 'set_global_optical_properties requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  recursive subroutine solve_pprts(solver, lthermal, lsolar, edirTOA, &
      & opt_solution_uid, opt_solution_time, opt_buildings)
    class(t_solver), intent(inout) :: solver
    logical, intent(in) :: lthermal, lsolar
    real(ireals), intent(in) :: edirTOA
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    real(ireals), optional, intent(in) :: opt_solution_time
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings
    call CHKERR(1_mpiint, 'solve_pprts requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine pprts_get_result(solver, redn, reup, rabso, redir, opt_solution_uid, opt_buildings)
    class(t_solver) :: solver
    real(ireals), dimension(:, :, :), intent(inout), allocatable :: redn, reup, rabso
    real(ireals), dimension(:, :, :), intent(inout), allocatable, optional :: redir
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    type(t_pprts_buildings), optional, intent(inout) :: opt_buildings
    call CHKERR(1_mpiint, 'pprts_get_result requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine pprts_get_result_toZero(solver, gedn, geup, gabso, gedir, opt_solution_uid)
    class(t_solver) :: solver
    real(ireals), intent(inout), dimension(:, :, :), allocatable :: gedn, geup, gabso
    real(ireals), intent(inout), dimension(:, :, :), allocatable, optional :: gedir
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    call CHKERR(1_mpiint, 'pprts_get_result_toZero requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine gather_all_toZero(C, inp, outp)
    type(t_coord), intent(in) :: C
    real(ireals), intent(in) :: inp(:, :, :)
    real(ireals), intent(inout), allocatable :: outp(:, :, :)
    call CHKERR(1_mpiint, 'gather_all_toZero requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine gather_all_to_all(C, inp, outp)
    type(t_coord), intent(in) :: C
    real(ireals), intent(in), allocatable :: inp(:, :, :)
    real(ireals), intent(inout), allocatable :: outp(:, :, :)
    call CHKERR(1_mpiint, 'gather_all_to_all requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

  subroutine scale_flx(solver, solution, lWm2)
    class(t_solver), intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution
    logical, intent(in) :: lWm2
    call CHKERR(1_mpiint, 'scale_flx requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

end module
