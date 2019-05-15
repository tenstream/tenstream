module m_pprts_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    zero, one, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_helper_functions, only : CHKWARN, CHKERR, get_arg, itoa

  use m_optprop, only: t_optprop_cube

  implicit none

  public :: t_solver, t_solver_1_2, t_solver_3_6, t_solver_3_10, &
    t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18, &
    allocate_pprts_solver_from_commandline, &
    t_coord, t_suninfo, &
    t_state_container, destroy_solution, &
    t_dof, t_solver_log_events, setup_log_events

  type t_coord
    integer(iintegers)      :: xs,xe                   ! local domain start and end indices
    integer(iintegers)      :: ys,ye                   ! local domain start and end indices
    integer(iintegers)      :: zs,ze                   ! local domain start and end indices
    integer(iintegers)      :: xm,ym,zm                ! size of local domain
    integer(iintegers)      :: gxs,gys,gzs             ! domain indices including ghost points
    integer(iintegers)      :: gxe,gye,gze             !
    integer(iintegers)      :: gxm,gym,gzm             ! size of local domain including ghosts

    integer(iintegers)      :: glob_xm,glob_ym,glob_zm ! global domain size

    integer(iintegers)      :: dof,dim                 ! degrees of freedom of Petsc Domain, dimension of dmda
    type(tDM)               :: da                      ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:)            ! all 3d neighbours((x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
    integer(mpiint)         :: comm                    ! mpi communicatior for this DMDA
  end type

  type t_atmosphere
    real(ireals), allocatable , dimension(:,:,:) :: planck, kabs, ksca, g
    real(ireals), allocatable , dimension(:,:,:) :: a11, a12, a21, a22, a13, a23, a33
    real(ireals), allocatable , dimension(:,:,:) :: dz
    logical     , allocatable , dimension(:,:,:) :: l1d
    real(ireals), allocatable , dimension(:,:)   :: albedo
    real(ireals), allocatable , dimension(:,:)   :: Btop, Bbot ! TOA layer planck emissions, special case memory for icollapse
    real(ireals)                                 :: dx,dy
    integer(iintegers)                           :: icollapse=1
    logical                                      :: lcollapse = .False.
  end type

  !type t_sunangles
  !  real(ireals)        :: symmetry_phi
  !  integer(iintegers)  :: yinc,xinc
  !  real(ireals)        :: theta, phi, costheta, sintheta
  !end type

  type t_suninfo
    !type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    real(ireals), allocatable, dimension(:,:,:) :: & ! C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye
      symmetry_phi, theta, phi, costheta, sintheta
    integer(iintegers), allocatable, dimension(:,:,:) :: xinc, yinc
    logical :: luse_topography=.False.
  end type

  type t_state_container
    integer(iintegers)  :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    type(tVec), allocatable    :: edir,ediff,abso

    logical             :: lset        = .False. ! initialized?
    logical             :: lsolar_rad  = .False. ! direct radiation calculated?
    logical             :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](false) or [W/m**2](true)
    logical             :: lWm2_dir=.False. , lWm2_diff=.False.

    !save error statistics
    real(ireals)        :: time   (30) = -one
    real(ireals)        :: maxnorm(30) = zero
    real(ireals)        :: twonorm(30) = zero
    real(ireals),allocatable :: dir_ksp_residual_history(:)
    real(ireals),allocatable :: diff_ksp_residual_history(:)
  end type

  type t_dof
    integer(iintegers) :: dof, area_divider, streams
    logical, allocatable :: is_inward(:)
  end type

  type t_solver_log_events
    PetscLogStage :: stage_solve
    PetscLogEvent :: set_optprop
    PetscLogEvent :: setup_dir_src
    PetscLogEvent :: compute_edir
    PetscLogEvent :: solve_Mdir
    PetscLogEvent :: setup_Mdir
    PetscLogEvent :: setup_diff_src
    PetscLogEvent :: compute_ediff
    PetscLogEvent :: solve_Mdiff
    PetscLogEvent :: setup_Mdiff
    PetscLogEvent :: compute_absorption

    PetscLogEvent :: solve_twostream
    PetscLogEvent :: solve_mcrts
    PetscLogEvent :: get_coeff_dir2dir
    PetscLogEvent :: get_coeff_dir2diff
    PetscLogEvent :: get_coeff_diff2diff
    PetscLogEvent :: get_result

    PetscLogEvent :: scatter_to_Zero
    PetscLogEvent :: scale_flx
    PetscLogEvent :: compute_orientation
    PetscLogEvent :: orient_face_normals
    PetscLogEvent :: debug_output
  end type

  type, abstract :: t_solver
    character(len=default_str_len)     :: solvername='' ! name to prefix e.g. log stages. If you create more than one solver, make sure that it has a unique name
    integer(mpiint)                    :: comm, myid, numnodes     ! mpi communicator, my rank and number of ranks in comm
    type(t_coord), allocatable         :: C_dir, C_diff, C_one, C_one1, C_one_atm, C_one_atm1
    type(t_atmosphere),allocatable     :: atm
    type(t_suninfo)                    :: sun
    type(tMat),allocatable             :: Mdir,Mdiff
    type(tKSP), allocatable            :: ksp_solar_dir
    type(tKSP), allocatable            :: ksp_solar_diff
    type(tKSP), allocatable            :: ksp_thermal_diff
    class(t_optprop_cube), allocatable :: OPP

    type(t_dof)                     :: difftop, diffside, dirtop, dirside

    logical                         :: lenable_solutions_err_estimates=.True.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
    type(tVec),allocatable          :: incSolar, b
    type(tVec),allocatable          :: dir_scalevec_Wm2_to_W, diff_scalevec_Wm2_to_W
    type(tVec),allocatable          :: dir_scalevec_W_to_Wm2, diff_scalevec_W_to_Wm2
    type(tVec),allocatable          :: abso_scalevec

    logical                         :: linitialized=.False.
    type(t_state_container)         :: solutions(-1000:1000)
    type(t_solver_log_events)       :: logs
  end type

  type, extends(t_solver) :: t_solver_1_2
  end type
  type, extends(t_solver) :: t_solver_3_6
  end type
  type, extends(t_solver) :: t_solver_3_10
  end type
  type, extends(t_solver) :: t_solver_8_10
  end type
  type, extends(t_solver) :: t_solver_3_16
  end type
  type, extends(t_solver) :: t_solver_8_16
  end type
  type, extends(t_solver) :: t_solver_8_18
  end type


  contains
    subroutine prepare_solution(edir_dm, ediff_dm, abso_dm, lsolar, solution, uid)
      type(tDM), intent(in) :: edir_dm, ediff_dm, abso_dm
      logical, intent(in) :: lsolar
      type(t_state_container), intent(inout) :: solution
      integer(iintegers), optional, intent(in) :: uid
      integer(mpiint) :: ierr

      if(solution%lset) call CHKERR(1_mpiint, 'solution has already been prepared before')

      solution%lset = .True.
      solution%lsolar_rad = lsolar

      solution%lchanged = .True.
      solution%lWm2_dir = .True.
      solution%lWm2_diff= .True.

      solution%uid = get_arg(i0, uid)

      if(solution%lsolar_rad) then
        allocate(solution%edir)
        call DMCreateGlobalVector(edir_dm, solution%edir, ierr)  ; call CHKERR(ierr)
        call PetscObjectSetName(solution%edir,'initialized_edir_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_init_edir", ierr); call CHKERR(ierr)
      endif

      allocate(solution%ediff)
      call DMCreateGlobalVector(ediff_dm, solution%ediff, ierr)  ; call CHKERR(ierr)
      call PetscObjectSetName(solution%ediff,'initialized_ediff_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_init_ediff", ierr); call CHKERR(ierr)

      allocate(solution%abso)
      call DMCreateGlobalVector(abso_dm, solution%abso, ierr)  ; call CHKERR(ierr)
      call PetscObjectSetName(solution%abso,'initialized_abso_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-show_init_abso", ierr); call CHKERR(ierr)

    end subroutine
    subroutine destroy_solution(solution)
      type(t_state_container), intent(inout) :: solution
      integer(mpiint) :: ierr
      if( solution%lset ) then
        if(solution%lsolar_rad) then
          if(allocated(solution%edir)) then
            call VecDestroy(solution%edir , ierr) ;call CHKERR(ierr)
            deallocate(solution%edir)
          endif
          solution%lsolar_rad = .False.
        endif

        if(allocated(solution%ediff)) then
          call VecDestroy(solution%ediff    , ierr) ;call CHKERR(ierr)
          deallocate(solution%ediff)
        endif
        if(allocated(solution%abso)) then
          call VecDestroy(solution%abso     , ierr) ;call CHKERR(ierr)
          deallocate(solution%abso)
        endif

        if(allocated(solution%dir_ksp_residual_history)) &
          deallocate(solution%dir_ksp_residual_history)
        if(allocated(solution%diff_ksp_residual_history)) &
          deallocate(solution%diff_ksp_residual_history)

        solution%lset = .False.
      endif
    end subroutine

    subroutine setup_log_events(logs, solvername)
      type(t_solver_log_events), intent(inout) :: logs
      character(len=*), optional :: solvername
      PetscClassId, parameter :: cid=0
      integer(mpiint) :: ierr
      character(len=default_str_len) :: s

      s = get_arg('pprts.', solvername)

      call PetscLogStageRegister(trim(s)//'solve_stage', logs%stage_solve, ierr); call CHKERR(ierr)

      call PetscLogEventRegister(trim(s)//'set_optprop', cid, logs%set_optprop, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'setup_dir_src', cid, logs%setup_dir_src, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'comp_Edir', cid, logs%compute_Edir, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'solve_Mdir', cid, logs%solve_Mdir, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'setup_Mdir', cid, logs%setup_Mdir, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'setup_diff_src', cid, logs%setup_diff_src, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'comp_Ediff', cid, logs%compute_Ediff, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'solve_Mdiff', cid, logs%solve_Mdiff, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'setup_Mdiff', cid, logs%setup_Mdiff, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'compute_absorption', cid, logs%compute_absorption, ierr); call CHKERR(ierr)

      call PetscLogEventRegister(trim(s)//'solve_twostr', cid, logs%solve_twostream, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'solve_mcrts', cid, logs%solve_mcrts, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'dir2dir', cid, logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'dir2diff', cid, logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'diff2diff', cid, logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'get_result', cid, logs%get_result, ierr); call CHKERR(ierr)

      call PetscLogEventRegister(trim(s)//'scttr2Zero', cid, logs%scatter_to_Zero, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'scale_flx', cid, logs%scale_flx, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'compute_orientation', cid, logs%compute_orientation, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'orient_face_normals', cid, logs%orient_face_normals, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'debug_output', cid, logs%debug_output, ierr); call CHKERR(ierr)
    end subroutine

  subroutine allocate_pprts_solver_from_commandline(pprts_solver, default_solver)
    class(t_solver), intent(inout), allocatable :: pprts_solver
    character(len=*), intent(in), optional :: default_solver

    logical :: lflg
    character(len=default_str_len) :: solver_str
    integer(mpiint) :: ierr

    if(allocated(pprts_solver)) then
      call CHKWARN(1_mpiint, 'called allocate_pprts_solver_from_commandline on an already allocated solver...'//&
        'have you been trying to change the solver type on the fly?'// &
        'this is not possible, please destroy the old one and create a new one')
      return
    endif

    solver_str = get_arg('none', trim(default_solver))
    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-solver', solver_str, lflg, ierr) ; call CHKERR(ierr)

    select case (solver_str)
      case('1_2')
        allocate(t_solver_1_2::pprts_solver)

      case('3_6')
        allocate(t_solver_3_6::pprts_solver)

      case('3_10')
        allocate(t_solver_3_10::pprts_solver)

      case('8_10')
        allocate(t_solver_8_10::pprts_solver)

      case('3_16')
        allocate(t_solver_3_16::pprts_solver)

      case('8_16')
        allocate(t_solver_8_16::pprts_solver)

      case('8_18')
        allocate(t_solver_8_18::pprts_solver)

      case default
        print *,'error, have to provide solver type as argument, e.g. call with'
        print *,'-solver 1_2'
        print *,'-solver 3_6'
        print *,'-solver 3_10'
        print *,'-solver 8_10'
        print *,'-solver 3_16'
        print *,'-solver 8_16'
        print *,'-solver 8_18'
        call CHKERR(1_mpiint, 'have to provide solver type')
    end select

  end subroutine
end module
