module m_pprts_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & ireals, irealLUT, iintegers, mpiint, &
    & zero, one, pi, &
    & i0, i1, i2, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKWARN, CHKERR, &
    & cstr, &
    & deallocate_allocatable, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_min, &
    & is_inrange, &
    & rel_approx, &
    & toStr

  use m_petsc_helpers, only: &
    & getvecpointer, restorevecpointer, is_local_vec

  use m_optprop, only: t_optprop_cube

  implicit none

  private
  public :: &
    & allocate_pprts_solver_from_commandline, &
    & atmk, &
    & compute_gradient, &
    & destroy_pprts, &
    & destroy_solution, &
    & determine_ksp_tolerances, &
    & get_coeff, &
    & get_solution_uid, &
    & interpolate_cell_values_to_vertices, &
    & prepare_solution, &
    & print_solution, &
    & set_dmda_cell_coordinates, &
    & setup_incSolar, &
    & setup_log_events, &
    & solver_to_str, &
    & t_atmosphere, &
    & t_coord, &
    & t_dof, &
    & t_mat_permute_info, &
    & t_pprts_shell_ctx, &
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
    & t_solver_log_events, &
    & t_solver_mcdmda, &
    & t_solver_rayli, &
    & t_state_container, &
    & t_suninfo

  type t_coord
    integer(iintegers) :: xs, xe                   ! local domain start and end indices
    integer(iintegers) :: ys, ye                   ! local domain start and end indices
    integer(iintegers) :: zs, ze                   ! local domain start and end indices
    integer(iintegers) :: xm, ym, zm                ! size of local domain
    integer(iintegers) :: gxs, gys, gzs             ! domain indices including ghost points
    integer(iintegers) :: gxe, gye, gze             !
    integer(iintegers) :: gxm, gym, gzm             ! size of local domain including ghosts

    integer(iintegers) :: glob_xm, glob_ym, glob_zm ! global domain size

    integer(iintegers) :: dof, dim                 ! degrees of freedom of Petsc Domain, dimension of dmda
    type(tDM) :: da                      ! The Domain Decomposition Object
    integer(mpiint) :: comm                    ! mpi communicatior for this DMDA
    integer(mpiint), allocatable :: neighbors(:)        ! all 3d neighbours((x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self
  end type

  type t_atmosphere
    real(ireals), allocatable, dimension(:, :, :) :: planck, kabs, ksca, g
    real(ireals), allocatable, dimension(:, :, :) :: a11, a12, a21, a22, a13, a23, a33
    real(ireals), allocatable, dimension(:, :, :) :: dz
    logical, allocatable, dimension(:) :: l1d
    real(ireals) :: unconstrained_fraction ! fraction of 3D voxels
    real(ireals), allocatable, dimension(:, :) :: albedo
    real(ireals), allocatable, dimension(:, :) :: Btop, Bbot ! TOA layer planck emissions, special case memory for icollapse
    real(ireals), allocatable, dimension(:, :) :: Bsrfc      ! Srfc planck emissions
    real(ireals) :: dx, dy
    integer(iintegers) :: icollapse = 1
    logical :: lcollapse = .false.
    type(tVec), allocatable :: hhl   ! vertical level heights, local vec on C_one_atm1_box
    type(tVec), allocatable :: vert_heights ! vertical level heights, local vec on Cvert_one_atm1
    type(tVec), allocatable :: hgrad ! horizontal gradient of heights, C_two1
  end type

  type t_suninfo
    !type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    real(ireals) :: sundir(3)
    real(ireals) :: &
      symmetry_phi, theta, phi, costheta, sintheta
    integer(iintegers) :: xinc, yinc
    logical :: luse_topography = .false.
  end type

  type t_state_container
    integer(iintegers) :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    type(tVec), allocatable :: edir, ediff, abso

    logical :: lset = .false. ! initialized?
    logical :: lsolar_rad = .false. ! direct radiation calculated?
    logical :: lthermal_rad = .false. ! thermal radiation calculated?
    logical :: lchanged = .true.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](false) or [W/m**2](true)
    logical :: lWm2_dir = .false., lWm2_diff = .false.

    !save error statistics
    real(ireals) :: time(30) = -one
    real(ireals) :: maxnorm(30) = zero
    real(ireals) :: dir_ksp_residual_history(100) = -one
    real(ireals) :: diff_ksp_residual_history(100) = -one

    integer(iintegers) :: Niter_dir = -1, Niter_diff = -1
  end type

  type t_dof
    integer(iintegers) :: dof, area_divider, streams
    logical, allocatable :: is_inward(:)
  end type

  type t_solver_log_events
    PetscLogEvent :: set_optprop
    PetscLogEvent :: setup_dir_src
    PetscLogEvent :: compute_edir
    PetscLogEvent :: solve_Mdir
    PetscLogEvent :: setup_Mdir
    PetscLogEvent :: permute_mat_dir
    PetscLogEvent :: permute_mat_gen_dir
    PetscLogEvent :: permute_mat_diff
    PetscLogEvent :: permute_mat_gen_diff
    PetscLogEvent :: setup_diff_src
    PetscLogEvent :: compute_ediff
    PetscLogEvent :: solve_Mdiff
    PetscLogEvent :: setup_Mdiff
    PetscLogEvent :: compute_absorption

    PetscLogEvent :: solve_twostream
    PetscLogEvent :: solve_schwarzschild
    PetscLogEvent :: solve_rayli
    PetscLogEvent :: solve_mcrts
    PetscLogEvent :: solve_disort
    PetscLogEvent :: rayli_tracing
    PetscLogEvent :: get_coeff_dir2dir
    PetscLogEvent :: get_coeff_dir2diff
    PetscLogEvent :: get_coeff_diff2diff
    PetscLogEvent :: get_result

    PetscLogEvent :: scatter_to_Zero
    PetscLogEvent :: scale_flx
    PetscLogEvent :: compute_orientation
    PetscLogEvent :: orient_face_normals
    PetscLogEvent :: debug_output
    PetscLogEvent :: setup_initial_guess
  end type

  type t_pprts_shell_ctx
    class(t_solver), pointer :: solver
  end type

  type t_mat_permute_info
    type(tIS) :: is ! col/row permutations
    logical :: rev_x, rev_y, rev_z
    logical :: zlast, switch_xy
  end type

  type, abstract :: t_solver
    character(len=default_str_len) :: prefix = '' ! name to prefix options
    character(len=default_str_len) :: solvername = '' ! name to prefix e.g. log stages. If you create more than one solver, make sure that it has a unique name
    integer(mpiint) :: comm, myid, numnodes     ! mpi communicator, my rank and number of ranks in comm
    logical :: lopen_bc = .false. ! switch if boundaries are cyclic or open
    type(t_coord), allocatable :: C_dir
    type(t_coord), allocatable :: C_diff
    type(t_coord), allocatable :: C_one
    type(t_coord), allocatable :: C_one1
    type(t_coord), allocatable :: C_one_atm
    type(t_coord), allocatable :: C_one_atm1
    type(t_coord), allocatable :: C_one_atm1_box
    type(t_coord), allocatable :: C_two1
    type(t_coord), allocatable :: Cvert_one_atm1
    type(t_coord), allocatable :: Csrfc_one
    type(t_atmosphere), allocatable :: atm
    type(t_suninfo) :: sun
    type(tMat), allocatable :: Mdir, Mdiff, Mth
    type(tMat), allocatable :: Mdir_perm, Mdiff_perm, Mth_perm
    type(t_mat_permute_info), allocatable :: perm_dir, perm_diff ! col/row permutations
    type(tKSP), allocatable :: ksp_solar_dir
    type(tKSP), allocatable :: ksp_solar_diff
    type(tKSP), allocatable :: ksp_thermal_diff
    class(t_optprop_cube), allocatable :: OPP
    class(t_optprop_cube), allocatable :: OPP1d

    type(t_dof) :: difftop, diffside, dirtop, dirside
    real(ireals), allocatable, dimension(:, :, :, :) :: dir2dir, dir2diff, diff2diff

    logical :: lenable_solutions_err_estimates = .true.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
    type(tVec), allocatable :: incSolar, b
    type(tVec), allocatable :: dir_scalevec_Wm2_to_W, diff_scalevec_Wm2_to_W
    type(tVec), allocatable :: dir_scalevec_W_to_Wm2, diff_scalevec_W_to_Wm2
    type(tVec), allocatable :: abso_scalevec

    logical :: linitialized = .false.
    type(t_solver_log_events) :: logs
    type(t_pprts_shell_ctx) :: shell_ctx
    type(t_state_container), allocatable :: solutions(:)
  end type

  type, extends(t_solver) :: t_solver_2str
  end type
  type, extends(t_solver) :: t_solver_disort
  end type
  type, extends(t_solver) :: t_solver_rayli
  end type
  type, extends(t_solver) :: t_solver_mcdmda
  end type
  type, extends(t_solver) :: t_solver_1_2
  end type
  type, extends(t_solver) :: t_solver_3_6
  end type
  type, extends(t_solver) :: t_solver_3_10
  end type
  type, extends(t_solver) :: t_solver_3_16
  end type
  type, extends(t_solver) :: t_solver_3_24
  end type
  type, extends(t_solver) :: t_solver_3_30
  end type
  type, extends(t_solver) :: t_solver_8_10
  end type
  type, extends(t_solver) :: t_solver_8_16
  end type
  type, extends(t_solver) :: t_solver_8_18
  end type

contains
  subroutine prepare_solution(edir_dm, ediff_dm, abso_dm, lsolar, lthermal, solution, uid)
    type(tDM), intent(in) :: edir_dm, ediff_dm, abso_dm
    logical, intent(in) :: lsolar, lthermal
    type(t_state_container), intent(inout) :: solution
    integer(iintegers), optional, intent(in) :: uid
    integer(mpiint) :: ierr

    if (solution%lset) call CHKERR(1_mpiint, 'solution has already been prepared before')

    solution%lsolar_rad = lsolar
    solution%lthermal_rad = lthermal

    solution%lchanged = .true.
    solution%lWm2_dir = .true.
    solution%lWm2_diff = .true.

    solution%uid = get_arg(i0, uid)

    if (solution%lsolar_rad) then
      allocate (solution%edir)
      call DMCreateGlobalVector(edir_dm, solution%edir, ierr); call CHKERR(ierr)
      call PetscObjectSetName(solution%edir, 'initialized_edir_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, edir_dm, "-show_init_edir", ierr); call CHKERR(ierr)
    end if

    allocate (solution%ediff)
    call DMCreateGlobalVector(ediff_dm, solution%ediff, ierr); call CHKERR(ierr)
    call PetscObjectSetName(solution%ediff, 'initialized_ediff_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%ediff, ediff_dm, "-show_init_ediff", ierr); call CHKERR(ierr)

    allocate (solution%abso)
    call DMCreateGlobalVector(abso_dm, solution%abso, ierr); call CHKERR(ierr)
    call PetscObjectSetName(solution%abso, 'initialized_abso_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
    call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%abso, abso_dm, "-show_init_abso", ierr); call CHKERR(ierr)

    solution%lset = .true.
  end subroutine
  subroutine destroy_solution(solution)
    type(t_state_container), intent(inout) :: solution
    if (solution%lset) then
      call deallocate_allocatable(solution%edir)
      solution%lsolar_rad = .false.
      solution%lthermal_rad = .false.

      call deallocate_allocatable(solution%ediff)
      call deallocate_allocatable(solution%abso)

      solution%lset = .false.
    end if
  end subroutine
  subroutine print_solution(solution)
    type(t_state_container), intent(inout) :: solution
    integer(mpiint) :: ierr
    character(len=30) :: header
    header = cstr('Solution('//toStr(solution%uid)//') ', 'blue')
    print *, trim(header)//'is initialized?', solution%lset
    if (.not. solution%lset) return
    print *, trim(header)//'is a solar solution?', solution%lsolar_rad
    print *, trim(header)//'has changed?', solution%lchanged
    print *, trim(header)//'direct  radiation is in W/m2?', solution%lWm2_dir
    print *, trim(header)//'diffuse radiation is in W/m2?', solution%lWm2_diff

    call investigate_vec(solution%edir, 'Edir ')
    call investigate_vec(solution%ediff, 'Ediff')
    call investigate_vec(solution%abso, 'Abso ')
  contains
    subroutine investigate_vec(v, title)
      type(tVec), allocatable, intent(in) :: v
      character(len=*), intent(in) :: title
      real(ireals) :: n2
      n2 = -1
      if (allocated(v)) then
        call VecNorm(v, NORM_2, n2, ierr); call CHKERR(ierr)
      end if
      print *, trim(header)//title// &
        ' (alloc='//toStr(allocated(v))//' 2-norm = '//toStr(n2)
    end subroutine
  end subroutine

  subroutine setup_log_events(logs, solvername)
    type(t_solver_log_events), intent(inout) :: logs
    character(len=*), optional :: solvername
    PetscClassId :: cid
    integer(mpiint) :: ierr
    character(len=default_str_len) :: s

    s = get_arg('pprts.', solvername)
    call PetscClassIdRegister(trim(s), cid, ierr); call CHKERR(ierr)

    call PetscLogEventRegister(trim(s)//'set_optprop', cid, logs%set_optprop, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'setup_dir_src', cid, logs%setup_dir_src, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'comp_Edir', cid, logs%compute_Edir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_Mdir', cid, logs%solve_Mdir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'setup_Mdir', cid, logs%setup_Mdir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'perm_mat_dir', cid, logs%permute_mat_dir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'perm_mat_gen_dir', cid, logs%permute_mat_gen_dir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'perm_mat_diff', cid, logs%permute_mat_diff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'perm_mat_gen_diff', cid, logs%permute_mat_gen_diff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'setup_diff_src', cid, logs%setup_diff_src, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'comp_Ediff', cid, logs%compute_Ediff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_Mdiff', cid, logs%solve_Mdiff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'setup_Mdiff', cid, logs%setup_Mdiff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'comp_abso', cid, logs%compute_absorption, ierr); call CHKERR(ierr)

    call PetscLogEventRegister(trim(s)//'solve_twostr', cid, logs%solve_twostream, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_schwarz', cid, logs%solve_schwarzschild, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_rayli', cid, logs%solve_rayli, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_disort', cid, logs%solve_disort, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'solve_mcrts', cid, logs%solve_mcrts, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'rayli_tracing', cid, logs%rayli_tracing, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'dir2dir', cid, logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'dir2diff', cid, logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'diff2diff', cid, logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'get_result', cid, logs%get_result, ierr); call CHKERR(ierr)

    call PetscLogEventRegister(trim(s)//'scttr2Zero', cid, logs%scatter_to_Zero, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'scale_flx', cid, logs%scale_flx, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'compute_orientation', cid, logs%compute_orientation, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'orient_face_normals', cid, logs%orient_face_normals, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'debug_output', cid, logs%debug_output, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'setup_initial_guess', cid, logs%setup_initial_guess, ierr); call CHKERR(ierr)
  end subroutine

  subroutine allocate_pprts_solver_from_commandline(pprts_solver, default_solver, ierr, prefix)
    class(t_solver), intent(inout), allocatable :: pprts_solver
    character(len=*), intent(in) :: default_solver
    integer(mpiint), intent(out) :: ierr
    character(len=*), intent(in), optional :: prefix

    logical :: lflg
    character(len=default_str_len) :: solver_str, pref

    ierr = 0

    if (allocated(pprts_solver)) then
      call CHKWARN(1_mpiint, 'called allocate_pprts_solver_from_commandline on an already allocated solver...'// &
                   'have you been trying to change the solver type on the fly?'// &
                   'this is not possible, please destroy the old one and create a new one')
      return
    end if

    solver_str = get_arg('none', trim(default_solver))
    pref = get_arg(PETSC_NULL_CHARACTER, prefix)
    call get_petsc_opt(pref, '-solver', solver_str, lflg, ierr); call CHKERR(ierr)

    select case (solver_str)
    case ('2str')
      allocate (t_solver_2str :: pprts_solver)

    case ('disort')
      allocate (t_solver_disort :: pprts_solver)

    case ('rayli')
      allocate (t_solver_rayli :: pprts_solver)

    case ('mcdmda')
      allocate (t_solver_mcdmda :: pprts_solver)

    case ('1_2')
      allocate (t_solver_1_2 :: pprts_solver)

    case ('3_6')
      allocate (t_solver_3_6 :: pprts_solver)

    case ('3_10')
      allocate (t_solver_3_10 :: pprts_solver)

    case ('8_10')
      allocate (t_solver_8_10 :: pprts_solver)

    case ('3_16')
      allocate (t_solver_3_16 :: pprts_solver)

    case ('3_24')
      allocate (t_solver_3_24 :: pprts_solver)

    case ('3_30')
      allocate (t_solver_3_30 :: pprts_solver)

    case ('8_16')
      allocate (t_solver_8_16 :: pprts_solver)

    case ('8_18')
      allocate (t_solver_8_18 :: pprts_solver)

    case default
      print *, 'error, have to provide solver type as argument, e.g. call with'
      print *, '-'//trim(pref)//'solver 2str'
      print *, '-'//trim(pref)//'solver disort'
      print *, '-'//trim(pref)//'solver 1_2'
      print *, '-'//trim(pref)//'solver 3_6'
      print *, '-'//trim(pref)//'solver 3_10'
      print *, '-'//trim(pref)//'solver 3_16'
      print *, '-'//trim(pref)//'solver 3_24'
      print *, '-'//trim(pref)//'solver 8_10'
      print *, '-'//trim(pref)//'solver 8_16'
      ierr = 1
      call CHKERR(ierr, 'have to provide valid solver type, '// &
        & 'given <'//trim(solver_str)//'> (prefix='//trim(pref)//')')
    end select

    pprts_solver%prefix = pref
  end subroutine

  function solver_to_str(solver) result(s)
    character(len=default_str_len) :: s
    class(t_solver), intent(in) :: solver

    select type (solver)
    type is (t_solver_2str)
      s = '2str'

    type is (t_solver_disort)
      s = 'disort'

    type is (t_solver_rayli)
      s = 'rayli'

    type is (t_solver_mcdmda)
      s = 'mcdmda'

    type is (t_solver_1_2)
      s = '1_2'

    type is (t_solver_3_6)
      s = '3_6'

    type is (t_solver_3_10)
      s = '3_10'

    type is (t_solver_8_10)
      s = '8_10'

    type is (t_solver_3_16)
      s = '3_16'

    type is (t_solver_3_24)
      s = '3_24'

    type is (t_solver_3_30)
      s = '3_30'

    type is (t_solver_8_16)
      s = '8_16'

    type is (t_solver_8_18)
      s = '8_18'

    class default
      s = 'NaN'
      call CHKERR(1_mpiint, 'could not determine description str for solver')
    end select
  end function

  subroutine destroy_coord(C)
    type(t_coord), allocatable, intent(inout) :: C
    integer(mpiint) :: ierr
    if (allocated(C)) then
      call DMDestroy(C%da, ierr); call CHKERR(ierr)
      deallocate (C)
    end if
  end subroutine
  subroutine destroy_atm(atm)
    type(t_atmosphere), intent(inout) :: atm
    call deallocate_allocatable(atm%planck)
    call deallocate_allocatable(atm%kabs)
    call deallocate_allocatable(atm%ksca)
    call deallocate_allocatable(atm%g)
    call deallocate_allocatable(atm%a11)
    call deallocate_allocatable(atm%a12)
    call deallocate_allocatable(atm%a21)
    call deallocate_allocatable(atm%a22)
    call deallocate_allocatable(atm%a13)
    call deallocate_allocatable(atm%a23)
    call deallocate_allocatable(atm%a33)
    call deallocate_allocatable(atm%dz)
    call deallocate_allocatable(atm%l1d)
    call deallocate_allocatable(atm%albedo)
    call deallocate_allocatable(atm%Btop)
    call deallocate_allocatable(atm%Bbot)
    call deallocate_allocatable(atm%Bsrfc)
    call deallocate_allocatable(atm%hhl)
    call deallocate_allocatable(atm%vert_heights)
    call deallocate_allocatable(atm%hgrad)
  end subroutine
  subroutine destroy_pprts(solver, lfinalizepetsc)
    class(t_solver) :: solver
    logical, optional :: lfinalizepetsc
    logical :: lfinalize, lpetsc_is_initialized
    integer(iintegers) :: uid
    integer(mpiint) :: ierr

    if (solver%linitialized) then
      call deallocate_allocatable(solver%ksp_solar_dir)
      call deallocate_allocatable(solver%ksp_solar_diff)
      call deallocate_allocatable(solver%ksp_thermal_diff)

      call deallocate_allocatable(solver%incSolar)
      call deallocate_allocatable(solver%b)

      call destroy_matrices(solver)

      if (allocated(solver%solutions)) then
        do uid = lbound(solver%solutions, 1), ubound(solver%solutions, 1)
          call destroy_solution(solver%solutions(uid))
        end do
        deallocate (solver%solutions)
      end if

      call deallocate_allocatable(solver%dir_scalevec_Wm2_to_W)
      call deallocate_allocatable(solver%diff_scalevec_Wm2_to_W)
      call deallocate_allocatable(solver%dir_scalevec_W_to_Wm2)
      call deallocate_allocatable(solver%diff_scalevec_W_to_Wm2)
      call deallocate_allocatable(solver%abso_scalevec)

      if (allocated(solver%atm)) then
        call destroy_atm(solver%atm)
        deallocate (solver%atm)
      end if

      if (allocated(solver%perm_dir)) then
        call ISDestroy(solver%perm_dir%is, ierr); call CHKERR(ierr)
        deallocate (solver%perm_dir)
      end if
      if (allocated(solver%perm_diff)) then
        call ISDestroy(solver%perm_diff%is, ierr); call CHKERR(ierr)
        deallocate (solver%perm_diff)
      end if

      if (allocated(solver%OPP)) then
        call solver%OPP%destroy(ierr); call CHKERR(ierr)
        deallocate (solver%OPP)
      end if
      if (allocated(solver%OPP1d)) then
        call solver%OPP1d%destroy(ierr); call CHKERR(ierr)
        deallocate (solver%OPP1d)
      end if

      call deallocate_allocatable(solver%dir2dir)
      call deallocate_allocatable(solver%dir2diff)
      call deallocate_allocatable(solver%diff2diff)

      call destroy_coord(solver%C_dir)
      call destroy_coord(solver%C_diff)
      call destroy_coord(solver%C_one)
      call destroy_coord(solver%C_one1)
      call destroy_coord(solver%C_one_atm)
      call destroy_coord(solver%C_one_atm1)
      call destroy_coord(solver%C_one_atm1_box)
      call destroy_coord(solver%C_two1)
      call destroy_coord(solver%Cvert_one_atm1)
      call destroy_coord(solver%Csrfc_one)

      call deallocate_allocatable(solver%difftop%is_inward)
      call deallocate_allocatable(solver%diffside%is_inward)
      call deallocate_allocatable(solver%dirtop%is_inward)
      call deallocate_allocatable(solver%dirside%is_inward)

      solver%comm = -1
      solver%linitialized = .false.
    end if

    lfinalize = get_arg(.false., lfinalizepetsc)
    if (lfinalize) then
      call PetscInitialized(lpetsc_is_initialized, ierr); call CHKERR(ierr)
      if (lpetsc_is_initialized) call PetscFinalize(ierr); call CHKERR(ierr)
    end if
  end subroutine

  subroutine destroy_matrices(solver)
    class(t_solver) :: solver

    call deallocate_allocatable(solver%Mdir)
    call deallocate_allocatable(solver%Mdiff)
    call deallocate_allocatable(solver%Mth)
    call deallocate_allocatable(solver%Mdir_perm)
    call deallocate_allocatable(solver%Mdiff_perm)
    call deallocate_allocatable(solver%Mth_perm)
  end subroutine

  !> @brief define physical coordinates for DMDA to allow for geometric multigrid
  subroutine set_dmda_cell_coordinates(solver, atm, da, ierr)
    class(t_solver), intent(in) :: solver
    type(t_atmosphere), intent(in) :: atm
    type(tDM), intent(in) :: da
    integer(mpiint), intent(out) :: ierr

    type(tVec) :: coordinates
    real(ireals), pointer, dimension(:, :, :, :) :: xv => null()
    real(ireals), pointer, dimension(:) :: xv1d => null()

    type(tDM) :: coordDA
    integer(iintegers) :: zs, zm, xs, xm, ys, ym, k, i, j

    real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()

    call DMGetCoordinatesLocal(da, coordinates, ierr); call CHKERR(ierr)
    if (coordinates .ne. PETSC_NULL_VEC) return

    if (.not. allocated(solver%atm%hhl)) &
      & call CHKERR(1_mpiint, 'solver%atm%hhl has to be allocated')

    call DMDASetUniformCoordinates(da, &
                                   zero, one, &
                                   zero, one, &
                                   zero, one, ierr); call CHKERR(ierr)
    call DMGetCoordinateDM(da, coordDA, ierr); call CHKERR(ierr)
    call DMDAGetGhostCorners(coordDA, zs, xs, ys, zm, xm, ym, ierr); call CHKERR(ierr)

    call DMGetCoordinatesLocal(da, coordinates, ierr); call CHKERR(ierr)
    call VecGetArrayF90(coordinates, xv1d, ierr); call CHKERR(ierr)
    xv(0:2, zs:zs + zm - 1, xs:xs + xm - 1, ys:ys + ym - 1) => xv1d

    call getVecPointer(solver%C_one_atm1_box%da, atm%hhl, hhl1d, hhl, readonly=.true.)
    do j = ys, ys + ym - 1
      do i = xs, xs + xm - 1
        do k = zs, zs + zm - 1
          xv(i0, k, i, j) = hhl(i0, atmk(solver%atm, k), i, j)
          xv(i1, k, i, j) = (real(i, ireals) + 0.5_ireals) * atm%dx
          xv(i2, k, i, j) = (real(j, ireals) + 0.5_ireals) * atm%dy
        end do
      end do
    end do
    nullify (xv)
    call restoreVecPointer(solver%C_one_atm1_box%da, atm%hhl, hhl1d, hhl, readonly=.true.)
    call VecRestoreArrayF90(coordinates, xv1d, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(coordinates, da, "-pprts_show_coordinates", ierr); call CHKERR(ierr)
  end subroutine

  !> @brief compute horizontal gradient from dz3d
  !> @details build horizontal gradient of height information, i.e. [dz/dx, dz/dy]
  subroutine compute_gradient(comm, atm, C_hhl, vhhl, C_grad, vgrad)
    integer(mpiint), intent(in) :: comm
    type(t_atmosphere), intent(in) :: atm
    type(tVec), intent(in) :: vhhl
    type(t_coord), intent(in) :: C_hhl, C_grad
    type(tVec), allocatable :: vgrad

    real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()
    real(ireals), pointer :: grad(:, :, :, :) => null(), grad_1d(:) => null()

    integer(iintegers) :: i, j, k

    real(ireals) :: zm(4)
    integer(mpiint) :: myid, ierr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    if (.not. allocated(atm%dz)) &
      & call CHKERR(1_mpiint, toStr(myid)//' You called  compute_gradient()'// &
      & ' but the atm struct is not yet up, make sure we have atm%dz before')

    if (.not. allocated(vgrad)) then
      allocate (vgrad)
      call DMCreateLocalVector(C_grad%da, vgrad, ierr); call CHKERR(ierr)
    end if

    call getVecPointer(C_hhl%da, vhhl, hhl1d, hhl, readonly=.true.)
    call getVecPointer(C_grad%da, vgrad, grad_1d, grad)

    do j = C_grad%ys, C_grad%ye
      do i = C_grad%xs, C_grad%xe
        do k = C_grad%zs, C_grad%ze
          ! Mean heights of adjacent columns
          zm(1) = hhl(i0, atmk(atm, k), i - 1, j)
          zm(2) = hhl(i0, atmk(atm, k), i + 1, j)

          zm(3) = hhl(i0, atmk(atm, k), i, j - 1)
          zm(4) = hhl(i0, atmk(atm, k), i, j + 1)

          ! Gradient of height field
          grad(i0, k, i, j) = (zm(2) - zm(1)) / (2._ireals * atm%dx)
          grad(i1, k, i, j) = (zm(4) - zm(3)) / (2._ireals * atm%dy)
        end do
      end do
    end do

    call restoreVecPointer(C_hhl%da, vhhl, hhl1d, hhl, readonly=.true.)

    call restoreVecPointer(C_grad%da, vgrad, grad_1d, grad)

    call PetscObjectViewFromOptions(vgrad, C_hhl%da, "-pprts_show_grad", ierr); call CHKERR(ierr)
  end subroutine

  !> @brief horizontally interpolate/average cell values onto vertices
  !> @details averages the 4 adjacent cell values onto a vertex.
  !>          C_cells should be a box stencil da
  !>          dof on C_cells has to be the same as on C_verts
  !>          and C_verts has to have Nz, Nx+1 and Ny+1 entries
  subroutine interpolate_cell_values_to_vertices(C_cells, cell_vals, C_verts, vert_vals)
    type(tVec), intent(in) :: cell_vals
    type(t_coord), intent(in) :: C_cells
    type(t_coord), intent(in) :: C_verts
    type(tVec), intent(inout) :: vert_vals

    type(tVec) :: vcells
    real(ireals), pointer :: xc(:, :, :, :) => null(), xc1d(:) => null()
    real(ireals), pointer :: xv(:, :, :, :) => null(), xv1d(:) => null()
    integer(iintegers) :: k, i, j
    logical :: is_local
    integer(mpiint) :: ierr

    call CHKERR(int(C_cells%glob_xm - C_verts%glob_xm + i1, mpiint), &
      & 'nonconforming size: cells/verts %xe ', C_cells%glob_xm, C_verts%glob_xm)
    call CHKERR(int(C_cells%glob_ym - C_verts%glob_ym + i1, mpiint), &
      & 'nonconforming size: cells/verts %ye ', C_cells%glob_ym, C_verts%glob_ym)
    call CHKERR(int(C_cells%ze - C_verts%ze, mpiint), &
      & 'nonconforming size: cells/verts %ze ', C_cells%ze, C_verts%ze)

    call CHKERR(int(C_cells%dof - C_verts%dof, mpiint), &
      & 'nonconforming size: cells/verts %dof ', C_cells%dof, C_verts%dof)

    call is_local_vec(C_cells%da, cell_vals, is_local, ierr); call CHKERR(ierr)
    if (is_local) then
      vcells = cell_vals
    else
      call DMGetLocalVector(C_cells%da, vcells, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_cells%da, cell_vals, INSERT_VALUES, vcells, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C_cells%da, cell_vals, INSERT_VALUES, vcells, ierr); call CHKERR(ierr)
    end if

    call getVecPointer(C_cells%da, vcells, xc1d, xc)
    call getVecPointer(C_verts%da, vert_vals, xv1d, xv)

    do j = C_verts%ys, C_verts%ye
      do i = C_verts%xs, C_verts%xe
        do k = C_verts%zs, C_verts%ze
          xv(:, k, i, j) = (xc(:, k, i - 1, j - 1) + xc(:, k, i, j - 1) &
                      & + xc(:, k, i - 1, j) + xc(:, k, i, j))*.25_ireals
        end do
      end do
    end do

    call restoreVecPointer(C_cells%da, vert_vals, xv1d, xv)
    call restoreVecPointer(C_verts%da, vcells, xc1d, xc)

    if (.not. is_local) then
      call DMRestoreLocalVector(C_cells%da, vcells, ierr); call CHKERR(ierr)
    end if
  end subroutine

  pure function atmk(atm, k) ! return vertical index value for DMDA grid on atmosphere grid
    integer(iintegers) :: atmk
    integer(iintegers), intent(in) :: k
    type(t_atmosphere), intent(in) :: atm
    atmk = k + atm%icollapse - 1
  end function

  !> @brief: determine tolerances for solvers
  subroutine determine_ksp_tolerances(C, unconstrained_fraction, rtol, atol, maxit, ksp)
    type(t_coord), intent(in) :: C
    real(ireals), intent(in) :: unconstrained_fraction
    real(ireals), intent(out) :: rtol, atol
    integer(iintegers), intent(out) :: maxit
    type(tKSP), intent(in), allocatable, optional :: ksp
    real(ireals) :: rel_atol = 1e-4_ireals
    integer(mpiint) :: myid, ierr
    real(ireals) :: dtol
    logical, parameter :: ldebug = .false.

    integer(iintegers), parameter :: maxiter = 1000

    if (present(ksp)) then
      if (allocated(ksp)) then
        call KSPGetTolerances(ksp, rtol, atol, dtol, maxit, ierr); call CHKERR(ierr)
        if (ldebug) then
          call mpi_comm_rank(C%comm, myid, ierr); call CHKERR(ierr)
          if (myid .eq. 0) print *, cstr('Read tolerances from ksp', 'red'), rtol, atol, dtol, maxit
        end if
        return
      end if
    end if

    maxit = maxiter
    rtol = 1e-5_ireals
    atol = rel_atol * real(C%glob_xm * C%glob_ym * C%glob_zm, ireals) * unconstrained_fraction
    atol = max(1e-8_ireals, atol)

    if (ldebug) then
      call mpi_comm_rank(C%comm, myid, ierr); call CHKERR(ierr)
      if (myid .eq. 0) &
        print *, 'KSP ', &
        & '-- tolerances:', rtol, atol, &
        & ':: rel_atol', rel_atol, &
        & ':: total dof', C%dof * C%glob_xm * C%glob_ym * C%glob_zm, &
        & ':: unconstrained fraction', unconstrained_fraction
    end if
  end subroutine

  !> @brief set solar incoming radiation at Top_of_Atmosphere
  !> @details todo: in case we do not have periodic boundaries, we should shine light in from the side of the domain...
  subroutine setup_incSolar(solver, edirTOA, incSolar)
    class(t_solver), target, intent(in) :: solver
    real(ireals), intent(in) :: edirTOA
    type(tVec), intent(inout) :: incSolar

    real(ireals), pointer :: x1d(:) => null(), x4d(:, :, :, :) => null()

    integer(mpiint) :: ierr
    real(ireals) :: fac
    integer(iintegers) :: i, j, src
    logical, parameter :: ldebug = .false.

    associate ( &
        & atm => solver%atm, &
        & sun => solver%sun, &
        & C_dir => solver%C_dir)

      fac = edirTOA * atm%dx * atm%dy / real(solver%dirtop%area_divider, ireals) * sun%costheta

      call VecSet(incSolar, 0._ireals, ierr); call CHKERR(ierr)
      call getVecPointer(C_dir%da, incSolar, x1d, x4d)

      do j = C_dir%ys, C_dir%ye
        do i = C_dir%xs, C_dir%xe
          do src = 0, solver%dirtop%dof - 1
            x4d(src, C_dir%zs, i, j) = fac
          end do
        end do
      end do

      call restoreVecPointer(C_dir%da, incSolar, x1d, x4d)

      if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'Setup of IncSolar done', edirTOA, &
        & '(', fac, ')'
    end associate

    call set_open_bc()

  contains

    subroutine set_open_bc()
      logical :: lsun_north, lsun_east
      integer(mpiint) :: ierr
      integer(iintegers) :: i, j, k, src, ioff
      type(tMat) :: A
      type(tVec) :: b, local_incSolar
      type(tKSP) :: ksp
      character(len=default_str_len) :: prefix
      real(ireals), pointer :: xv_b(:)
      real(ireals), pointer :: xg1d(:) => null(), xg4d(:, :, :, :) => null()

      if (solver%lopen_bc) then ! need to update side fluxes somewhere in the domain

        associate ( &
            & atm => solver%atm, &
            & sun => solver%sun, &
            & C_dir => solver%C_dir)

          call DMGetLocalVector(C_dir%da, local_incSolar, ierr); call CHKERR(ierr)
          call VecSet(local_incSolar, 0._ireals, ierr); call CHKERR(ierr)
          !call DMGlobalToLocalBegin(C_dir%da, incSolar, INSERT_VALUES, local_incSolar, ierr); call CHKERR(ierr)
          !call DMGlobalToLocalEnd(C_dir%da, incSolar, INSERT_VALUES, local_incSolar, ierr); call CHKERR(ierr)

          call getVecPointer(C_dir%da, local_incSolar, x1d, x4d)

          ! note here we copy the local parts directly, otherwise, with GlobalToLocal we would create TOA incoming energy
          ! at the domain edges which leads to double counting later when doing the communication with add values
          call getVecPointer(C_dir%da, incSolar, xg1d, xg4d, readonly=.true.)
          x4d(0:solver%dirtop%dof - 1, C_dir%zs, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye) = &
            & xg4d(0:solver%dirtop%dof - 1, C_dir%zs, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye)
          call restoreVecPointer(C_dir%da, incSolar, xg1d, xg4d, readonly=.true.)

          lsun_north = sun%yinc .eq. i0
          lsun_east = sun%xinc .eq. i0

          if (C_dir%xs .eq. i0 .or. &
            & C_dir%ys .eq. i0 .or. &
            & C_dir%xe + 1 .eq. C_dir%glob_xm .or. &
            & C_dir%ye + 1 .eq. C_dir%glob_ym) then ! have an outer domain boundary

            call MatCreate(PETSC_COMM_SELF, A, ierr); call CHKERR(ierr)
            call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE,&
              & C_dir%dof * C_dir%zm, C_dir%dof * C_dir%zm, ierr); call CHKERR(ierr)

            if (len_trim(solver%solvername) .gt. 0) then
              prefix = trim(solver%solvername)//'_open_bc_'
            else
              prefix = 'open_bc_'
            end if
            call MatSetOptionsPrefix(A, prefix, ierr); call CHKERR(ierr)

            call MatSetFromOptions(A, ierr); call CHKERR(ierr)
            call MatSeqAIJSetPreallocation(A, C_dir%dof + i1, PETSC_NULL_INTEGER, ierr); call CHKERR(ierr)

            call MatSetUp(A, ierr); call CHKERR(ierr)

            call MatCreateVecs(A, b, PETSC_NULL_VEC, ierr); call CHKERR(ierr)

            call KSPCreate(PETSC_COMM_SELF, ksp, ierr); call CHKERR(ierr)

            if (C_dir%ys .eq. i0 .or. C_dir%ye + 1 .eq. C_dir%glob_ym) then ! we have open boundary in north/south
              if (lsun_north) then
                j = C_dir%ye
              else
                j = C_dir%ys
              end if
              do i = C_dir%xs, C_dir%xe

                call single_column_solve(solver%OPP, C_dir, i, j, A, ksp, b)

                call VecGetArrayReadF90(b, xv_b, ierr)
                do k = C_dir%zs, C_dir%ze - 1
                  do src = 0, solver%dirside%dof - 1
                    ioff = solver%dirtop%dof + solver%dirside%dof + src
                    x4d(ioff, k, i, j + 1 - sun%yinc) = fac * xv_b(i1 + k * C_dir%dof + ioff)
                  end do
                end do
                call VecRestoreArrayReadF90(b, xv_b, ierr)
              end do
            end if

            if (C_dir%xs .eq. i0 .or. C_dir%xe + 1 .eq. C_dir%glob_xm) then ! we have open boundary in east / west
              if (lsun_east) then
                i = C_dir%xe
              else
                i = C_dir%xs
              end if
              do j = C_dir%ys, C_dir%ye

                call single_column_solve(solver%OPP, C_dir, i, j, A, ksp, b)

                call VecGetArrayReadF90(b, xv_b, ierr)
                do k = C_dir%zs, C_dir%ze - 1
                  do src = 0, solver%dirside%dof - 1
                    ioff = solver%dirtop%dof + src
                    x4d(ioff, k, i + 1 - sun%xinc, j) = fac * xv_b(i1 + k * C_dir%dof + ioff)
                  end do
                end do
                call VecRestoreArrayReadF90(b, xv_b, ierr)
              end do
            end if

            call MatDestroy(A, ierr); call CHKERR(ierr)
            call VecDestroy(b, ierr); call CHKERR(ierr)
            call KSPDestroy(ksp, ierr); call CHKERR(ierr)
          end if ! have an outer domain boundary

          call restoreVecPointer(C_dir%da, local_incSolar, x1d, x4d)

          call VecSet(incSolar, 0._ireals, ierr); call CHKERR(ierr)
          call DMLocalToGlobalBegin(C_dir%da, local_incSolar, ADD_VALUES, incSolar, ierr); call CHKERR(ierr)
          call DMLocalToGlobalEnd(C_dir%da, local_incSolar, ADD_VALUES, incSolar, ierr); call CHKERR(ierr)
          call DMRestoreLocalVector(C_dir%da, local_incSolar, ierr); call CHKERR(ierr)

        end associate

      end if
    end subroutine

    subroutine single_column_solve(OPP, C_dir, i, j, A, ksp, b)
      class(t_optprop_cube), intent(in) :: OPP
      type(t_coord), intent(in) :: C_dir
      integer(iintegers), intent(in) :: i, j
      type(tMat), intent(inout) :: A
      type(tVec), intent(inout) :: b
      type(tKSP), intent(inout) :: ksp

      integer(mpiint) :: ierr
      integer(iintegers) :: k, ak, src, dst, irow, icol, ioff, ioffsrc
      real(ireals) :: dtau, v
      real(irealLUT) :: lutcoeff(C_dir%dof**2)
      real(ireals), target :: coeff(C_dir%dof**2)
      real(ireals), pointer :: cdir2dir(:, :), xv_b(:)

      call MatZeroEntries(A, ierr); call CHKERR(ierr)
      do irow = 0, C_dir%dof * C_dir%zm - 1
        call MatSetValues(A, i1, irow, i1, irow, -1._ireals, ADD_VALUES, ierr); call CHKERR(ierr)
      end do

      associate ( &
          & atm => solver%atm, &
          & sun => solver%sun)

        call VecGetArrayF90(b, xv_b, ierr)
        xv_b(:) = 0
        do src = 1, solver%dirtop%dof
          xv_b(src) = -1._ireals
        end do
        call VecRestoreArrayF90(b, xv_b, ierr)

        do k = C_dir%zs, C_dir%ze - 1
          ak = atmk(atm, k)

          if (atm%l1d(atmk(atm, k))) then

            do src = 0, solver%dirtop%dof - 1
              irow = (k + 1) * C_dir%dof + src
              icol = (k) * C_dir%dof + src

              dtau = atm%kabs(ak, i, j) * atm%dz(ak, i, j) / sun%costheta
              v = exp(-dtau)
              call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
            end do
          else

            call get_coeff( &
              & OPP, &
              & atm%kabs(ak, i, j), &
              & atm%ksca(ak, i, j), &
              & atm%g(ak, i, j), &
              & atm%dz(ak, i, j), &
              & atm%dx, &
              & .true., &
              & lutcoeff, &
              & [real(sun%symmetry_phi, irealLUT), real(sun%theta, irealLUT)], &
              & lswitch_east=sun%xinc .eq. 0, &
              & lswitch_north=sun%yinc .eq. 0 &
              & )
            coeff = real(lutcoeff, ireals)
            cdir2dir(0:C_dir%dof - 1, 0:C_dir%dof - 1) => coeff(:)
            ! cdir2dir(0:C_dir%dof - 1, 0:C_dir%dof - 1) => solver%dir2dir(:, k, i, j)

            do src = 0, solver%dirtop%dof - 1
              icol = k * C_dir%dof + src
              do dst = 0, solver%dirtop%dof - 1 ! top2bot
                v = cdir2dir(src, dst)
                irow = (k + 1) * C_dir%dof + dst
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do

              do dst = 0, solver%dirside%dof - 1 ! top2x
                ioff = solver%dirtop%dof + dst
                v = cdir2dir(src, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
              do dst = 0, solver%dirside%dof - 1 ! top2y
                ioff = solver%dirtop%dof + solver%dirside%dof + dst
                v = cdir2dir(src, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
            end do

            do src = 0, solver%dirside%dof - 1
              ioffsrc = solver%dirtop%dof + src
              icol = k * C_dir%dof + ioffsrc
              do dst = 0, solver%dirtop%dof - 1 ! side2bot
                v = cdir2dir(ioffsrc, dst)
                irow = (k + 1) * C_dir%dof + dst
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do

              do dst = 0, solver%dirside%dof - 1 ! side2x
                ioff = solver%dirtop%dof + dst
                v = cdir2dir(ioffsrc, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
              do dst = 0, solver%dirside%dof - 1 ! side2y
                ioff = solver%dirtop%dof + solver%dirside%dof + dst
                v = cdir2dir(ioffsrc, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
            end do

            do src = 0, solver%dirside%dof - 1
              ioffsrc = solver%dirtop%dof + solver%dirside%dof + src
              icol = k * C_dir%dof + ioffsrc
              do dst = 0, solver%dirtop%dof - 1 ! side2bot
                v = cdir2dir(ioffsrc, dst)
                irow = (k + 1) * C_dir%dof + dst
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do

              do dst = 0, solver%dirside%dof - 1 ! side2x
                ioff = solver%dirtop%dof + dst
                v = cdir2dir(ioffsrc, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
              do dst = 0, solver%dirside%dof - 1 ! side2y
                ioff = solver%dirtop%dof + solver%dirside%dof + dst
                v = cdir2dir(ioffsrc, ioff)
                irow = k * C_dir%dof + ioff
                call MatSetValues(A, i1, irow, i1, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
              end do
            end do

          end if
        end do
        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)

        call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)
        call KSPSolve(ksp, b, b, ierr); call CHKERR(ierr)
      end associate

    end subroutine

  end subroutine

  function get_solution_uid(solutions, opt_solution_uid) result(uid)
    type(t_state_container), allocatable :: solutions(:)
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    integer(iintegers) :: uid
    logical :: lflg1, lflg2
    integer(mpiint) :: ierr

    uid = get_arg(0_iintegers, opt_solution_uid)
    call get_petsc_opt('', '-override_solution_uid', uid, lflg1, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-pprts_override_solution_uid', uid, lflg2, ierr); call CHKERR(ierr)
    !if (lflg1 .or. lflg2) then
    !  print *, 'Override solutions uid, returning '//toStr(uid)//' instead of '//toStr(get_arg(0_iintegers, opt_solution_uid))
    !end if

    if (.not. is_inrange(uid, lbound(solutions, 1, kind=iintegers), ubound(solutions, 1, kind=iintegers))) then
      call CHKWARN(int(uid, mpiint), "uid ("//toStr(uid)//") is not in range of "// &
        & "preallocated solutions container [ "//&
        & toStr(lbound(solutions, 1))//", "//&
        & toStr(ubound(solutions, 1))//" ]."//new_line('')// &
        & "I will set it to 0 but this has some implications"//new_line('')// &
        & "If you are not sure what you are doing, "// &
        & "you can increase the size of the solutions container to fit all spectral bands")
      uid = 0
    end if
  end function

  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(OPP, kabs, ksca, g, dz, dx, ldir, coeff, &
                       angles, lswitch_east, lswitch_north, opt_vertices)
    class(t_optprop_cube), intent(in) :: OPP
    real(ireals), intent(in) :: kabs, ksca, g, dz, dx
    logical, intent(in) :: ldir
    real(irealLUT), intent(out) :: coeff(:)

    real(irealLUT), intent(in), optional :: angles(2)
    logical, intent(in), optional :: lswitch_east, lswitch_north
    real(ireals), intent(in), optional :: opt_vertices(:)

    logical, parameter :: lenable_cache = .false.
    real(irealLUT), save :: coeff_cache(1000)
    logical :: lcurrent

    real(irealLUT) :: aspect_zx, tauz, w0
    integer(mpiint) :: ierr

    if (lenable_cache) then
      call check_cache(lcurrent)
      if (lcurrent) then
        coeff = coeff_cache(1:size(coeff))
        return
      end if
    end if

    aspect_zx = real(dz / dx, irealLUT)
    w0 = real(ksca / max(kabs + ksca, epsilon(kabs)), irealLUT)
    tauz = real((kabs + ksca) * dz, irealLUT)

    if (present(angles)) then
      aspect_zx = max(OPP%dev%dirconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(OPP%dev%dirconfig%dims(1)%vrange(1), &
           & min(OPP%dev%dirconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%dev%dirconfig%dims(2)%vrange(1), &
           & min(OPP%dev%dirconfig%dims(2)%vrange(2), w0))
    else
      aspect_zx = max(OPP%dev%diffconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(OPP%dev%diffconfig%dims(1)%vrange(1), &
           & min(OPP%dev%diffconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%dev%diffconfig%dims(2)%vrange(1), &
           & min(OPP%dev%diffconfig%dims(2)%vrange(2), w0))
    end if
    !print *,'hit cache'//new_line(''),  &
    !  & 'kabs', kabs,new_line(''), &
    !  & 'ksca', ksca,new_line(''), &
    !  & 'g',    g   ,new_line(''), &
    !  & 'dz',   dz  ,new_line(''), &
    !  & 'ldir', ldir

    call OPP%get_coeff(tauz, w0, real(g, irealLUT), aspect_zx, ldir, coeff, ierr, &
                       angles, lswitch_east, lswitch_north, opt_vertices); call CHKERR(ierr)

    if (lenable_cache) then
      coeff_cache(1:size(coeff)) = coeff
    end if

  contains
    subroutine check_cache(lcurrent)
      logical, intent(out) :: lcurrent

      logical, save :: lpresent_angles = .false.
      logical, save :: c_ldir = .false., c_lswitch_east = .false., c_lswitch_north = .false.
      real(ireals), save :: c_kabs = -1, c_ksca = -1, c_g = -1, c_dz = -1, c_vertices(24) = -1
      real(irealLUT), save :: c_angles(2) = -1

      real(ireals), parameter :: cache_limit = 1e-3 ! rel change cache limit
      real(irealLUT), parameter :: cache_limit2 = cache_limit ! rel change cache limit

      if (.not. lenable_cache) then
        lcurrent = .false.
        return
      end if

      if (.not. rel_approx(c_kabs, kabs, cache_limit)) goto 99
      if (.not. rel_approx(c_ksca, ksca, cache_limit)) goto 99
      if (.not. rel_approx(c_g, g, cache_limit)) goto 99
      if (.not. rel_approx(c_dz, dz, cache_limit)) goto 99
      if (lpresent_angles .neqv. present(angles)) goto 99

      if (present(opt_vertices)) then
        if (any(.not. rel_approx(c_vertices, opt_vertices, cache_limit))) goto 99
      end if

      if (present(angles)) then
        if (any(.not. rel_approx(c_angles, angles, cache_limit2))) goto 99
      end if

      if (c_ldir .neqv. ldir) goto 99
      if (present(lswitch_east)) then
        if (c_lswitch_east .neqv. lswitch_east) goto 99
      end if
      if (present(lswitch_north)) then
        if (c_lswitch_north .neqv. lswitch_north) goto 99
      end if
      lcurrent = .true.
      !print *,'hit cache'//new_line(''),  &
      !  & 'kabs', c_kabs, kabs,new_line(''), &
      !  & 'ksca', c_ksca, ksca,new_line(''), &
      !  & 'g',    c_g   , g   ,new_line(''), &
      !  & 'dz',   c_dz  , dz  ,new_line(''), &
      !  & 'ldir', c_ldir, ldir
      !call CHKERR(1_mpiint, 'DEBUG')
      return

99    continue ! update sample pts and go on to compute them
      !print *,'missed cache'//new_line(''),  &
      !  & 'kabs', c_kabs, kabs,new_line(''), &
      !  & 'ksca', c_ksca, ksca,new_line(''), &
      !  & 'g',    c_g   , g   ,new_line(''), &
      !  & 'dz',   c_dz  , dz  ,new_line(''), &
      !  & 'ldir', c_ldir, ldir

      lcurrent = .false.
      c_kabs = kabs
      c_ksca = ksca
      c_g = g
      c_dz = dz
      c_ldir = ldir
      lpresent_angles = present(angles)
      if (present(angles)) c_angles = angles
      if (present(opt_vertices)) c_vertices = opt_vertices
      if (present(lswitch_east)) c_lswitch_east = lswitch_east
      if (present(lswitch_north)) c_lswitch_north = lswitch_north

    end subroutine
  end subroutine
end module
