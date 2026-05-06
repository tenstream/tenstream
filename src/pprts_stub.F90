! Stub module providing m_pprts public interface when PETSc is not available.
! All subroutines abort at runtime with a clear error message.
module m_pprts
  use m_data_parameters, only: iintegers, ireals, irealLUT, mpiint
  use m_helper_functions, only: CHKERR
  use m_pprts_base, only: t_solver, t_coord, t_state_container
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
    call CHKERR(1_mpiint, 'init_pprts requires PETSc -- rebuild with -DWITH_PETSC=ON')
  end subroutine

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

  subroutine set_angles(solver, sundir)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: sundir(:)
    call CHKERR(1_mpiint, 'set_angles requires PETSc -- rebuild with -DWITH_PETSC=ON')
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
