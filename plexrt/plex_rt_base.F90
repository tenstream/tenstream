module m_plex_rt_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: CHKERR, CHKWARN, get_arg, itoa
  use m_data_parameters, only : ireals, iintegers, mpiint, default_str_len, &
    zero, one, i0
  use m_pprts_base, only : t_solver_log_events, t_state_container
  use m_plex_grid, only: t_plexgrid
  use m_optprop, only : t_optprop_wedge

  implicit none

  private
  public :: t_plex_solver, &
    t_plex_solver_5_8, &
    t_plex_solver_rectilinear_5_8, &
    t_plex_solver_18_8, &
    t_field, &
    allocate_plexrt_solver_from_commandline, &
    t_state_container_plexrt, prepare_solution, destroy_solution, &
    setup_field_IS

  type, extends(t_state_container) :: t_state_container_plexrt
    type(tVec), allocatable :: flx,abso

    ! save state of solution vectors... they are either in [W](false) or [W/m**2](true)
    logical :: lWm2=.False.
    logical :: lsolar_rad=.False.

    !save error statistics
    real(ireals),allocatable :: ksp_residual_history(:)

    integer(iintegers) :: Niter=-1
  end type

  type t_field
    integer(iintegers) :: id, dof, area_divider, streams
    character(len=default_str_len) :: fieldname
    type(tIS), allocatable :: subIS ! encodes the global idx of the subproblem
    type(tDM), allocatable :: subdm

    ! idx of the participating points
    type(tIS), allocatable :: cellIS, topfaceIS, botfaceIS

    type(tDM), allocatable :: dofdm  ! hold section for in/out dof ISes
    type(tVec),allocatable :: in_out_dof ! dim cellIS[cStart..cEnd-1]
  end type

  type, abstract :: t_plex_solver
    type(t_field), allocatable :: fields(:)
    type(t_plexgrid), allocatable :: plex
    class(t_optprop_wedge), allocatable :: OPP

    type(tVec), allocatable :: kabs, ksca, g       ! in each cell [pStart..pEnd-1]
    type(tVec), allocatable :: albedo              ! on each surface face [defined on plex%srfc_boundary_dm]

    type(tVec), allocatable :: plck ! Planck Radiation in each vertical level [W] [fStart..pEnd-1]

    type(t_state_container_plexrt) :: solutions(-1000:1000)

    type(tVec), allocatable :: dirsrc, diffsrc
    type(tMat), allocatable :: Mdir
    type(tMat), allocatable :: Mdiff
    type(tKSP), allocatable :: ksp_solar_dir
    type(tKSP), allocatable :: ksp_solar_diff
    type(tKSP), allocatable :: ksp_thermal_diff

    type(tVec), allocatable :: scalevec_Wm2_to_W, scalevec_W_to_Wm2
    !type(tIS),  allocatable :: IS_diff_in_out_dof

    logical :: lenable_solutions_err_estimates=.True.

    type(t_solver_log_events) :: logs

    type(tDM), allocatable :: mergedm
  end type

  type, extends(t_plex_solver) :: t_plex_solver_5_8
  end type
  type, extends(t_plex_solver) :: t_plex_solver_rectilinear_5_8
  end type
  type, extends(t_plex_solver) :: t_plex_solver_18_8
  end type

  contains

    subroutine allocate_plexrt_solver_from_commandline(plexrt_solver, default_solver)
      class(t_plex_solver), intent(inout), allocatable :: plexrt_solver
      character(len=*), intent(in), optional :: default_solver

      logical :: lflg
      character(len=default_str_len) :: solver_str
      integer(mpiint) :: ierr

      if(allocated(plexrt_solver)) then
        call CHKWARN(1_mpiint, 'called allocate_plexrt_solver_from_commandline on an already allocated solver...'//&
          'have you been trying to change the solver type on the fly?'// &
          'this is not possible, please destroy the old one and create a new one')
        return
      endif

      solver_str = get_arg('none', trim(default_solver))
      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-solver', solver_str, lflg, ierr) ; call CHKERR(ierr)

      select case (solver_str)
      case('5_8')
        allocate(t_plex_solver_5_8::plexrt_solver)

      case('rectilinear_5_8')
        allocate(t_plex_solver_rectilinear_5_8::plexrt_solver)

      case('18_8')
        allocate(t_plex_solver_18_8::plexrt_solver)

      case default
        print *,'error, have to provide solver type as argument, e.g. call with'
        print *,'-solver 5_8'
        print *,'-solver rectilinear_5_8'
        print *,'-solver 18_8'
        call CHKERR(1_mpiint, 'have to provide solver type')
      end select

    end subroutine

    subroutine prepare_solution(dm, abso_dm, solution, uid)
      type(tDM), intent(in) :: dm, abso_dm
      type(t_state_container_plexrt), intent(inout) :: solution
      integer(iintegers), optional, intent(in) :: uid
      integer(mpiint) :: ierr

      if(solution%lset) call CHKERR(1_mpiint, 'solution has already been prepared before')

      solution%lchanged = .True.
      solution%lWm2 = .True.

      solution%uid = get_arg(i0, uid)

      allocate(solution%flx)
      call DMCreateGlobalVector(dm, solution%flx, ierr)  ; call CHKERR(ierr)
      call PetscObjectSetName(solution%flx,'initialized_flx_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%flx, zero, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%flx, dm, "-show_init_flx", ierr); call CHKERR(ierr)

      allocate(solution%abso)
      call DMCreateGlobalVector(abso_dm, solution%abso, ierr)  ; call CHKERR(ierr)
      call PetscObjectSetName(solution%abso,'initialized_abso_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-show_init_abso", ierr); call CHKERR(ierr)

      solution%lset = .True.
    end subroutine
    subroutine destroy_solution(solution)
      type(t_state_container_plexrt), intent(inout) :: solution
      integer(mpiint) :: ierr
      if( solution%lset ) then
        if(allocated(solution%flx)) then
          call VecDestroy(solution%flx, ierr) ;call CHKERR(ierr)
          deallocate(solution%flx)
        endif
        if(allocated(solution%abso)) then
          call VecDestroy(solution%abso     , ierr) ;call CHKERR(ierr)
          deallocate(solution%abso)
        endif

        if(allocated(solution%ksp_residual_history)) &
          deallocate(solution%ksp_residual_history)

        solution%lset = .False.
      endif
    end subroutine

    subroutine setup_field_IS(pntStart, pntEnd, mask, fIS, optname)
      integer(iintegers), intent(in) :: pntStart, pntEnd
      logical :: mask(pntStart:)
      type(tIS), allocatable, intent(inout) :: fIS
      character(len=*), intent(in), optional :: optname
      integer(iintegers) :: i, j, numidx
      integer(iintegers), allocatable :: idx(:)

      integer(mpiint) :: ierr

      if(allocated(fIS)) call CHKERR(1_mpiint, 'fIS already allocated?')
      allocate(fIS)

      call CHKERR(int(size(mask)-pntEnd+pntStart,mpiint), &
        'input sizes dont match: pntStart/End '//itoa(pntStart)//'/'//itoa(pntEnd)// &
        ' ( '//itoa(pntEnd-pntStart)//' ) '// &
        'vs size(mask): '//itoa(size(mask)))

      numidx = count(mask)
      allocate(idx(max(1_iintegers,numidx)), source=-1)

      j = 1
      do i = pntStart, pntEnd-1
        if(mask(i)) then
          idx(j) = i
          j = j+1
        endif
      enddo

      call ISCreateGeneral(PETSC_COMM_SELF, numidx, idx, PETSC_COPY_VALUES, fIS, ierr); call CHKERR(ierr)

      if(present(optname)) then
        call PetscObjectSetName(fIS, trim(optname), ierr); call CHKERR(ierr)
      endif
    end subroutine
end module
