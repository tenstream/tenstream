module m_plex_rt_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: CHKERR, CHKWARN, get_arg, get_petsc_opt
  use m_data_parameters, only: ireals, iintegers, mpiint, default_str_len
  use m_pprts_base, only: t_state_container, t_solver_log_events
  use m_plex_grid, only: t_plexgrid
  use m_optprop, only: t_optprop_wedge

  implicit none

  private
  public :: t_plex_solver, &
            t_plex_solver_2str, &
            t_plex_solver_disort, &
            t_plex_solver_5_8, &
            t_plex_solver_rectilinear_5_8, &
            t_plex_solver_18_8, &
            t_plex_solver_rayli, &
            t_dof, &
            allocate_plexrt_solver_from_commandline

  type t_dof
    integer(iintegers) :: dof, area_divider, streams
  end type

  type, abstract :: t_plex_solver
    type(t_dof) :: difftop, diffside, dirtop, dirside
    integer(iintegers) :: dirdof, diffdof
    type(t_plexgrid), allocatable :: plex
    class(t_optprop_wedge), allocatable :: OPP

    real(ireals), allocatable :: dir2dir(:, :)      ! in each cell [0:dirdof**2-1     , pStart..pEnd-1]
    real(ireals), allocatable :: dir2diff(:, :)     ! in each cell [0:dirdof*diffdof-1, pStart..pEnd-1]
    real(ireals), allocatable :: diff2diff(:, :)    ! in each cell [0:diffdof**2-1    , pStart..pEnd-1]

    type(tVec), allocatable :: kabs, ksca, g       ! in each cell [pStart..pEnd-1]
    type(tVec), allocatable :: albedo              ! on each surface face [defined on plex%srfc_boundary_dm]

    type(tVec), allocatable :: plck ! Planck Radiation in each vertical level [W] [fStart..pEnd-1]

    type(t_state_container), allocatable :: solutions(:)

    type(tVec), allocatable :: dirsrc, diffsrc
    type(tMat), allocatable :: Mdir
    type(tMat), allocatable :: Mdiff
    type(tKSP), allocatable :: ksp_solar_dir
    type(tKSP), allocatable :: ksp_solar_diff
    type(tKSP), allocatable :: ksp_thermal_diff

    type(tVec), allocatable :: dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2
    type(tVec), allocatable :: diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2
    type(tIS), allocatable :: IS_diff_in_out_dof

    logical :: lenable_solutions_err_estimates = .true.

    type(t_solver_log_events) :: logs
  end type

  type, extends(t_plex_solver) :: t_plex_solver_2str
  end type
  type, extends(t_plex_solver) :: t_plex_solver_disort
  end type
  type, extends(t_plex_solver) :: t_plex_solver_5_8
  end type
  type, extends(t_plex_solver) :: t_plex_solver_rectilinear_5_8
  end type
  type, extends(t_plex_solver) :: t_plex_solver_18_8
  end type
  type, extends(t_plex_solver) :: t_plex_solver_rayli
  end type

contains

  subroutine allocate_plexrt_solver_from_commandline(plexrt_solver, default_solver)
    class(t_plex_solver), intent(inout), allocatable :: plexrt_solver
    character(len=*), intent(in), optional :: default_solver

    logical :: lflg
    character(len=default_str_len) :: solver_str
    integer(mpiint) :: ierr

    if (allocated(plexrt_solver)) then
      call CHKWARN(1_mpiint, 'called allocate_plexrt_solver_from_commandline on an already allocated solver...'// &
                   'have you been trying to change the solver type on the fly?'// &
                   'this is not possible, please destroy the old one and create a new one')
      return
    end if

    solver_str = get_arg('none', trim(default_solver))
    call get_petsc_opt(PETSC_NULL_CHARACTER, '-solver', solver_str, lflg, ierr); call CHKERR(ierr)

    select case (solver_str)
    case ('2str')
      allocate (t_plex_solver_2str :: plexrt_solver)

    case ('disort')
      allocate (t_plex_solver_disort :: plexrt_solver)

    case ('5_8')
      allocate (t_plex_solver_5_8 :: plexrt_solver)

    case ('rectilinear_5_8')
      allocate (t_plex_solver_rectilinear_5_8 :: plexrt_solver)

    case ('18_8')
      allocate (t_plex_solver_18_8 :: plexrt_solver)

    case ('rayli')
      allocate (t_plex_solver_rayli :: plexrt_solver)

    case default
      print *, 'error, have to provide solver type as argument, e.g. call with'
      print *, '-solver 2str'
      print *, '-solver disort'
      print *, '-solver 5_8'
      print *, '-solver rectilinear_5_8'
      print *, '-solver 18_8'
      print *, '-solver rayli'
      call CHKERR(1_mpiint, 'have to provide solver type')
    end select

  end subroutine

end module
