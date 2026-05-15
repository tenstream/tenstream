module m_pprts_base
#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#endif

  use m_data_parameters, only: &
    & ireals, irealLUT, iintegers, mpiint, &
    & imp_ireals, &
    & zero, one, &
    & i0, i1, i2, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKWARN, CHKERR, &
    & cstr, &
    & deallocate_allocatable, &
    & get_arg, &
    & get_petsc_opt, &
    & is_inrange, &
    & rel_approx, &
    & toStr

#ifdef HAVE_PETSC
  use m_petsc_helpers, only: &
    & getvecpointer, restorevecpointer, is_local_vec
#endif

  use m_optprop, only: t_optprop_cube

  implicit none

  private
  public :: &
    & allocate_pprts_solver_from_commandline, &
    & atmk, &
    & destroy_pprts, &
    & destroy_solution, &
    & determine_ksp_tolerances, &
    & get_coeff, &
    & get_solution_uid, &
    & halo_fill_5pt, &
    & halo_reduce_5pt, &
    & prepare_solution, &
    & print_solution, &
    & setup_coord_native, &
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
#ifdef HAVE_PETSC
  public :: &
    & compute_gradient, &
    & interpolate_cell_values_to_vertices, &
    & set_dmda_cell_coordinates
#endif

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
#ifdef HAVE_PETSC
    type(tDM) :: da                      ! The Domain Decomposition Object
#endif
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
#ifdef HAVE_PETSC
    type(tVec), allocatable :: hhl   ! vertical level heights, local vec on C_one_atm1_box
    type(tVec), allocatable :: vert_heights ! vertical level heights, local vec on Cvert_one_atm1
    type(tVec), allocatable :: hgrad ! horizontal gradient of heights, C_two1
#endif
  end type

  type t_suninfo
    !type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    real(ireals) :: sundir(3)
    real(ireals) :: symmetry_phi, theta, phi, costheta ! note that phi/theta may be rounded or overwritten
    real(ireals) :: mu ! original cos(theta) used for normalization of incoming solar
    integer(iintegers) :: xinc, yinc
    logical :: luse_topography = .false.
  end type

  type t_state_container
    integer(iintegers) :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
#ifdef HAVE_PETSC
    type(tVec), allocatable :: edir   ! global Vec on C_dir DMDA
    type(tVec), allocatable :: ediff  ! global Vec on C_diff DMDA
    type(tVec), allocatable :: abso   ! global Vec on C_one DMDA
#else
    real(ireals), allocatable :: edir(:, :, :, :)   ! (0:dof-1, zs:ze, xs:xe, ys:ye)
    real(ireals), allocatable :: ediff(:, :, :, :)  ! (0:dof-1, zs:ze, xs:xe, ys:ye)
    real(ireals), allocatable :: abso(:, :, :, :)   ! (0:dof-1, zs:ze, xs:xe, ys:ye)
#endif

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

    real(ireals) :: diff_sor_omega = -one  ! warm-start omega for adaptive SOR (-1 = unset)
  end type

  type t_dof
    integer(iintegers) :: dof, area_divider, streams
    logical, allocatable :: is_inward(:)
  end type

  type t_solver_log_events
#ifdef HAVE_PETSC
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
#endif
  end type

  type t_pprts_shell_ctx
    class(t_solver), pointer :: solver
  end type

  type t_mat_permute_info
#ifdef HAVE_PETSC
    type(tIS) :: is ! col/row permutations
#endif
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
#ifdef HAVE_PETSC
    type(tMat), allocatable :: Mdir, Mdiff, Mth
    type(tMat), allocatable :: Mdir_perm, Mdiff_perm, Mth_perm
    type(t_mat_permute_info), allocatable :: perm_dir, perm_diff ! col/row permutations
    type(tKSP), allocatable :: ksp_solar_dir
    type(tKSP), allocatable :: ksp_solar_diff
    type(tKSP), allocatable :: ksp_thermal_diff
#endif
    class(t_optprop_cube), allocatable :: OPP
    class(t_optprop_cube), allocatable :: OPP1d

    type(t_dof) :: difftop, diffside, dirtop, dirside
    real(ireals), allocatable, dimension(:, :, :, :) :: dir2dir, dir2diff, diff2diff

    logical :: lenable_solutions_err_estimates = .true.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
#ifdef HAVE_PETSC
    type(tVec), allocatable :: incSolar       ! global Vec on C_dir DMDA
    type(tVec), allocatable :: b              ! global Vec on C_diff DMDA
    type(tVec), allocatable :: dir_scalevec_Wm2_to_W
    type(tVec), allocatable :: diff_scalevec_Wm2_to_W
    type(tVec), allocatable :: dir_scalevec_W_to_Wm2
    type(tVec), allocatable :: diff_scalevec_W_to_Wm2
    type(tVec), allocatable :: abso_scalevec
#else
    real(ireals), allocatable :: incSolar(:, :, :, :)         ! (0:dof-1, zs:ze, xs:xe, ys:ye) on C_dir
    real(ireals), allocatable :: b(:, :, :, :)                ! (0:dof-1, zs:ze, xs:xe, ys:ye) on C_diff
    real(ireals), allocatable :: dir_scalevec_Wm2_to_W(:, :, :, :)
    real(ireals), allocatable :: diff_scalevec_Wm2_to_W(:, :, :, :)
    real(ireals), allocatable :: dir_scalevec_W_to_Wm2(:, :, :, :)
    real(ireals), allocatable :: diff_scalevec_W_to_Wm2(:, :, :, :)
    real(ireals), allocatable :: abso_scalevec(:, :, :, :)
#endif

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
#ifdef HAVE_PETSC
  subroutine prepare_solution(edir_dm, ediff_dm, abso_dm, lsolar, lthermal, solution, uid)
    type(tDM), intent(in) :: edir_dm, ediff_dm, abso_dm
#else
    subroutine prepare_solution(C_dir, C_diff, C_one, lsolar, lthermal, solution, uid)
      type(t_coord), intent(in) :: C_dir, C_diff, C_one
#endif
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

#ifdef HAVE_PETSC
      if (solution%lsolar_rad) then
        allocate (solution%edir)
        call DMCreateGlobalVector(edir_dm, solution%edir, ierr); call CHKERR(ierr)
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
      end if
      allocate (solution%ediff)
      call DMCreateGlobalVector(ediff_dm, solution%ediff, ierr); call CHKERR(ierr)
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      allocate (solution%abso)
      call DMCreateGlobalVector(abso_dm, solution%abso, ierr); call CHKERR(ierr)
      call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)
#else
      if (solution%lsolar_rad) then
        allocate (solution%edir(0:C_dir%dof - 1, C_dir%zs:C_dir%ze, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye))
        solution%edir = zero
      end if
      allocate (solution%ediff(0:C_diff%dof - 1, C_diff%zs:C_diff%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye))
      solution%ediff = zero
      allocate (solution%abso(0:C_one%dof - 1, C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye))
      solution%abso = zero
#endif

      solution%lset = .true.
    end subroutine
    subroutine destroy_solution(solution)
      type(t_state_container), intent(inout) :: solution
      integer(mpiint) :: ierr
      if (solution%lset) then
#ifdef HAVE_PETSC
        if (allocated(solution%edir)) then
          call VecDestroy(solution%edir, ierr); call CHKERR(ierr)
          deallocate (solution%edir)
        end if
        call VecDestroy(solution%ediff, ierr); call CHKERR(ierr)
        deallocate (solution%ediff)
        call VecDestroy(solution%abso, ierr); call CHKERR(ierr)
        deallocate (solution%abso)
#else
        call deallocate_allocatable(solution%edir)
        call deallocate_allocatable(solution%ediff)
        call deallocate_allocatable(solution%abso)
#endif
        solution%lsolar_rad = .false.
        solution%lthermal_rad = .false.
        solution%lset = .false.
      end if
    end subroutine
    subroutine print_solution(solution)
      type(t_state_container), intent(inout) :: solution
      character(len=30) :: header
      header = cstr('Solution('//toStr(solution%uid)//') ', 'blue')
      print *, trim(header)//'is initialized?', solution%lset
      if (.not. solution%lset) return
      print *, trim(header)//'is a solar solution?', solution%lsolar_rad
      print *, trim(header)//'has changed?', solution%lchanged
      print *, trim(header)//'direct  radiation is in W/m2?', solution%lWm2_dir
      print *, trim(header)//'diffuse radiation is in W/m2?', solution%lWm2_diff

#ifdef HAVE_PETSC
      call investigate_vec(solution%edir, 'Edir ')
      call investigate_vec(solution%ediff, 'Ediff')
      call investigate_vec(solution%abso, 'Abso ')
#else
      call investigate_arr(solution%edir, 'Edir ')
      call investigate_arr(solution%ediff, 'Ediff')
      call investigate_arr(solution%abso, 'Abso ')
#endif
    contains
#ifdef HAVE_PETSC
      subroutine investigate_vec(v, title)
        type(tVec), allocatable, intent(in) :: v
        character(len=*), intent(in) :: title
        real(ireals) :: n2
        integer(mpiint) :: ierr
        n2 = -1
        if (allocated(v)) call VecNorm(v, NORM_2, n2, ierr)
        print *, trim(header)//title// &
          ' (alloc='//toStr(allocated(v))//' 2-norm = '//toStr(n2)
      end subroutine
#else
      subroutine investigate_arr(v, title)
        real(ireals), allocatable, intent(in) :: v(:, :, :, :)
        character(len=*), intent(in) :: title
        real(ireals) :: n2
        n2 = -1
        if (allocated(v)) n2 = norm2(v)
        print *, trim(header)//title// &
          ' (alloc='//toStr(allocated(v))//' 2-norm = '//toStr(n2)
      end subroutine
#endif
    end subroutine

    subroutine setup_log_events(logs, solvername)
      type(t_solver_log_events), intent(inout) :: logs
      character(len=*), optional :: solvername
#ifdef HAVE_PETSC
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
#endif
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
      pref = get_arg('', prefix)
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

    !> @brief Fill a t_coord without PETSc using MPI Cartesian topology.
    !> @details Replaces DMDACreate3d+setup_coords for non-PETSc builds.
    !>   Domain layout mirrors DMDA: z is never decomposed (each rank owns full column).
    !>   Ghost stencil width is 1 in x and y.  neighbors() uses the same index
    !>   formula as DMDAGetNeighbors: idx = (1+diz) + (1+dix)*3 + (1+diy)*9.
    subroutine setup_coord_native(icomm, Nz, Nx, Ny, dof, C, nxproc, nyproc)
      use mpi, only: MPI_PROC_NULL, MPI_Dims_create, MPI_Cart_create, MPI_Cart_get, MPI_Cart_shift, MPI_Comm_free
      integer(mpiint), intent(in) :: icomm
      integer(iintegers), intent(in) :: Nz, Nx, Ny, dof
      type(t_coord), allocatable, intent(inout) :: C
      integer(iintegers), optional, intent(in) :: nxproc(:), nyproc(:)

      integer(mpiint) :: cart_comm, ierr, myid, nproc
      integer(mpiint) :: dims(2), coords(2)
      logical :: periods(2)
      integer(mpiint) :: rank_west, rank_east, rank_south, rank_north
      integer(iintegers) :: nxp, nyp, xi, yi

      if (allocated(C)) deallocate (C)
      allocate (C)

      C%comm = icomm
      C%dof = dof
      C%dim = 3
      C%glob_zm = Nz
      C%glob_xm = Nx
      C%glob_ym = Ny

      call mpi_comm_rank(icomm, myid, ierr); call CHKERR(ierr)

      if (present(nxproc) .and. present(nyproc)) then
        nxp = size(nxproc, kind=iintegers)
        nyp = size(nyproc, kind=iintegers)
        dims = int([nxp, nyp], mpiint)
      else
        call mpi_comm_size(icomm, nproc, ierr); call CHKERR(ierr)
        dims = 0_mpiint
        call MPI_Dims_create(nproc, 2_mpiint, dims, ierr); call CHKERR(ierr)
        nxp = int(dims(1), iintegers)
        nyp = int(dims(2), iintegers)
      end if

      periods = .true.  ! match PETSc DM_BOUNDARY_PERIODIC in x and y
      call MPI_Cart_create(icomm, 2_mpiint, dims, periods, .false., cart_comm, ierr); call CHKERR(ierr)
      call MPI_Cart_get(cart_comm, 2_mpiint, dims, periods, coords, ierr); call CHKERR(ierr)

      xi = int(coords(1), iintegers)
      yi = int(coords(2), iintegers)

      ! z: never decomposed — full column on every rank
      C%zs = i0
      C%zm = Nz
      C%ze = Nz - i1
      C%gzs = i0
      C%gzm = Nz
      C%gze = Nz - i1

      ! x: use provided distribution or even split
      if (present(nxproc)) then
        C%xs = int(sum(nxproc(1:xi)), iintegers)
        C%xm = nxproc(xi + i1)
      else
        C%xs = (xi * Nx) / nxp
        C%xm = ((xi + i1) * Nx) / nxp - C%xs
      end if
      C%xe = C%xs + C%xm - i1

      ! y: use provided distribution or even split
      if (present(nyproc)) then
        C%ys = int(sum(nyproc(1:yi)), iintegers)
        C%ym = nyproc(yi + i1)
      else
        C%ys = (yi * Ny) / nyp
        C%ym = ((yi + i1) * Ny) / nyp - C%ys
      end if
      C%ye = C%ys + C%ym - i1

      ! ghost bounds: 1-cell stencil in x and y
      C%gxs = C%xs - i1
      C%gxe = C%xe + i1
      C%gxm = C%gxe - C%gxs + i1
      C%gys = C%ys - i1
      C%gye = C%ye + i1
      C%gym = C%gye - C%gys + i1

      ! neighbors(0:26): DMDA-compatible layout
      !   idx = (1+diz)*1 + (1+dix)*3 + (1+diy)*9
      !   self=13, west=10, east=16, south=4, north=22
      allocate (C%neighbors(0:3**C%dim - 1))
      C%neighbors = int(MPI_PROC_NULL, iintegers)
      C%neighbors(13) = int(myid, iintegers)

      call MPI_Cart_shift(cart_comm, 0_mpiint, 1_mpiint, rank_west, rank_east, ierr); call CHKERR(ierr)
      call MPI_Cart_shift(cart_comm, 1_mpiint, 1_mpiint, rank_south, rank_north, ierr); call CHKERR(ierr)

      C%neighbors(10) = int(rank_west, iintegers)
      C%neighbors(16) = int(rank_east, iintegers)
      C%neighbors(4) = int(rank_south, iintegers)
      C%neighbors(22) = int(rank_north, iintegers)

      call MPI_Comm_free(cart_comm, ierr); call CHKERR(ierr)
    end subroutine

    subroutine destroy_coord(C)
      type(t_coord), allocatable, intent(inout) :: C
      integer(mpiint) :: ierr
      if (allocated(C)) then
#ifdef HAVE_PETSC
        call DMDestroy(C%da, ierr); call CHKERR(ierr)
#endif
        if (allocated(C%neighbors)) deallocate (C%neighbors)
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
#ifdef HAVE_PETSC
      call deallocate_allocatable(atm%hhl)
      call deallocate_allocatable(atm%vert_heights)
      call deallocate_allocatable(atm%hgrad)
#endif
    end subroutine
    subroutine destroy_pprts(solver, lfinalizepetsc)
      class(t_solver) :: solver
      logical, optional :: lfinalizepetsc
      logical :: lfinalize
#ifdef HAVE_PETSC
      PetscBool :: lpetsc_is_initialized
#endif
      integer(iintegers) :: uid
      integer(mpiint) :: ierr

      if (solver%linitialized) then
#ifdef HAVE_PETSC
        call deallocate_allocatable(solver%ksp_solar_dir)
        call deallocate_allocatable(solver%ksp_solar_diff)
        call deallocate_allocatable(solver%ksp_thermal_diff)
#endif

#ifdef HAVE_PETSC
        if (allocated(solver%incSolar)) then
          call VecDestroy(solver%incSolar, ierr); call CHKERR(ierr)
          deallocate (solver%incSolar)
        end if
        if (allocated(solver%b)) then
          call VecDestroy(solver%b, ierr); call CHKERR(ierr)
          deallocate (solver%b)
        end if
#else
        call deallocate_allocatable(solver%incSolar)
        call deallocate_allocatable(solver%b)
#endif

        call destroy_matrices(solver)

        if (allocated(solver%solutions)) then
          do uid = lbound(solver%solutions, 1), ubound(solver%solutions, 1)
            call destroy_solution(solver%solutions(uid))
          end do
          deallocate (solver%solutions)
        end if

#ifdef HAVE_PETSC
        if (allocated(solver%dir_scalevec_Wm2_to_W)) then
          call VecDestroy(solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
          deallocate (solver%dir_scalevec_Wm2_to_W)
        end if
        if (allocated(solver%diff_scalevec_Wm2_to_W)) then
          call VecDestroy(solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
          deallocate (solver%diff_scalevec_Wm2_to_W)
        end if
        if (allocated(solver%dir_scalevec_W_to_Wm2)) then
          call VecDestroy(solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
          deallocate (solver%dir_scalevec_W_to_Wm2)
        end if
        if (allocated(solver%diff_scalevec_W_to_Wm2)) then
          call VecDestroy(solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
          deallocate (solver%diff_scalevec_W_to_Wm2)
        end if
        if (allocated(solver%abso_scalevec)) then
          call VecDestroy(solver%abso_scalevec, ierr); call CHKERR(ierr)
          deallocate (solver%abso_scalevec)
        end if
#else
        call deallocate_allocatable(solver%dir_scalevec_Wm2_to_W)
        call deallocate_allocatable(solver%diff_scalevec_Wm2_to_W)
        call deallocate_allocatable(solver%dir_scalevec_W_to_Wm2)
        call deallocate_allocatable(solver%diff_scalevec_W_to_Wm2)
        call deallocate_allocatable(solver%abso_scalevec)
#endif

        if (allocated(solver%atm)) then
          call destroy_atm(solver%atm)
          deallocate (solver%atm)
        end if

#ifdef HAVE_PETSC
        if (allocated(solver%perm_dir)) then
          call ISDestroy(solver%perm_dir%is, ierr); call CHKERR(ierr)
          deallocate (solver%perm_dir)
        end if
        if (allocated(solver%perm_diff)) then
          call ISDestroy(solver%perm_diff%is, ierr); call CHKERR(ierr)
          deallocate (solver%perm_diff)
        end if
#endif

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
#ifdef HAVE_PETSC
      if (lfinalize) then
        call PetscInitialized(lpetsc_is_initialized, ierr); call CHKERR(ierr)
        if (lpetsc_is_initialized) call PetscFinalize(ierr); call CHKERR(ierr)
      end if
#endif
    end subroutine

    subroutine destroy_matrices(solver)
      class(t_solver) :: solver
#ifdef HAVE_PETSC
      call deallocate_allocatable(solver%Mdir)
      call deallocate_allocatable(solver%Mdiff)
      call deallocate_allocatable(solver%Mth)
      call deallocate_allocatable(solver%Mdir_perm)
      call deallocate_allocatable(solver%Mdiff_perm)
      call deallocate_allocatable(solver%Mth_perm)
#endif
    end subroutine

#ifdef HAVE_PETSC
    !> @brief define physical coordinates for DMDA to allow for geometric multigrid
    subroutine set_dmda_cell_coordinates(solver, atm, da, ierr)
      class(t_solver), intent(in) :: solver
      type(t_atmosphere), intent(in) :: atm
      type(tDM), intent(in) :: da
      integer(mpiint), intent(out) :: ierr

      type(tVec) :: coordinates
      real(ireals), pointer, dimension(:, :, :, :) :: xv
      real(ireals), pointer, dimension(:) :: xv1d

      type(tDM) :: coordDA
      integer(iintegers) :: zs, zm, xs, xm, ys, ym, k, i, j

      real(ireals), pointer :: hhl(:, :, :, :), hhl1d(:)

      xv => null()
      xv1d => null()
      hhl => null()
      hhl1d => null()

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
      call VecGetArray(coordinates, xv1d, ierr); call CHKERR(ierr)
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
      call VecRestoreArray(coordinates, xv1d, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(PetscObjectCast(coordinates), PetscObjectCast(da), "-pprts_show_coordinates", ierr)
      call CHKERR(ierr)
    end subroutine
#endif

#ifdef HAVE_PETSC
    !> @brief compute horizontal gradient from dz3d
    !> @details build horizontal gradient of height information, i.e. [dz/dx, dz/dy]
    subroutine compute_gradient(comm, atm, C_hhl, vhhl, C_grad, vgrad)
      integer(mpiint), intent(in) :: comm
      type(t_atmosphere), intent(in) :: atm
      type(tVec), intent(in) :: vhhl
      type(t_coord), intent(in) :: C_hhl, C_grad
      type(tVec), allocatable :: vgrad

      real(ireals), pointer :: hhl(:, :, :, :), hhl1d(:)
      real(ireals), pointer :: grad(:, :, :, :), grad_1d(:)

      integer(iintegers) :: i, j, k

      real(ireals) :: zm(4)
      integer(mpiint) :: myid, ierr

      hhl => null()
      hhl1d => null()
      grad => null()
      grad_1d => null()

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

      call PetscObjectViewFromOptions(PetscObjectCast(vgrad), PetscObjectCast(C_hhl%da), "-pprts_show_grad", ierr)
      call CHKERR(ierr)
    end subroutine
#endif

#ifdef HAVE_PETSC
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
      real(ireals), pointer :: xc(:, :, :, :), xc1d(:)
      real(ireals), pointer :: xv(:, :, :, :), xv1d(:)
      integer(iintegers) :: k, i, j
      logical :: is_local
      integer(mpiint) :: ierr

      xc => null()
      xc1d => null()
      xv => null()
      xv1d => null()

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
#endif

    pure function atmk(atm, k) ! return vertical index value for DMDA grid on atmosphere grid
      integer(iintegers) :: atmk
      integer(iintegers), intent(in) :: k
      type(t_atmosphere), intent(in) :: atm
      atmk = k + atm%icollapse - 1
    end function

    !> @brief: determine tolerances for solvers
#ifdef HAVE_PETSC
    subroutine determine_ksp_tolerances(C, unconstrained_fraction, rtol, atol, maxit, ksp)
      type(tKSP), intent(in), allocatable, optional :: ksp
#else
      subroutine determine_ksp_tolerances(C, unconstrained_fraction, rtol, atol, maxit)
#endif
        type(t_coord), intent(in) :: C
        real(ireals), intent(in) :: unconstrained_fraction
        real(ireals), intent(out) :: rtol, atol
        integer(iintegers), intent(out) :: maxit
        real(ireals) :: rel_atol = 1e-4_ireals
        integer(mpiint) :: myid, ierr
#ifdef HAVE_PETSC
        real(ireals) :: dtol
#endif
        logical, parameter :: ldebug = .false.

        integer(iintegers), parameter :: maxiter = 1000

#ifdef HAVE_PETSC
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
#endif

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
#ifdef HAVE_PETSC
        type(tVec), intent(inout) :: incSolar
#else
        real(ireals), target, contiguous, intent(inout) :: incSolar(:, :, :, :)
#endif

        integer(mpiint) :: ierr
        real(ireals) :: fac
        integer(iintegers) :: i, j, src
        logical, parameter :: ldebug = .false.
#ifdef HAVE_PETSC
        real(ireals), pointer :: x1d(:), x4d(:, :, :, :)
#else
        real(ireals), pointer :: x4d(:, :, :, :)
#endif

        associate ( &
            & atm => solver%atm, &
            & sun => solver%sun, &
            & C_dir => solver%C_dir)

          fac = edirTOA * atm%dx * atm%dy / real(solver%dirtop%area_divider, ireals)

#ifdef HAVE_PETSC
          x1d => null()
          x4d => null()
          call VecSet(incSolar, zero, ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, incSolar, x1d, x4d)

          do j = C_dir%ys, C_dir%ye
            do i = C_dir%xs, C_dir%xe
              do src = 0, solver%dirtop%dof - 1
                x4d(src, C_dir%zs, i, j) = fac
              end do
            end do
          end do

          call restoreVecPointer(C_dir%da, incSolar, x1d, x4d)
#else
          incSolar = zero
          x4d => null()
          x4d(0:C_dir%dof - 1, C_dir%zs:C_dir%ze, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye) => incSolar

          do j = C_dir%ys, C_dir%ye
            do i = C_dir%xs, C_dir%xe
              do src = 0, solver%dirtop%dof - 1
                x4d(src, C_dir%zs, i, j) = fac
              end do
            end do
          end do
          nullify (x4d)
#endif

          if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'Setup of IncSolar done', edirTOA, &
            & '(', fac, ')'
        end associate

#ifdef HAVE_PETSC
        call set_open_bc()
#endif

      contains

#ifdef HAVE_PETSC
        subroutine set_open_bc()
          logical :: lsun_north, lsun_east
          integer(mpiint) :: ierr
          integer(iintegers) :: i, j, k, src, ioff
          type(tMat) :: A
          type(tVec) :: b, local_incSolar, vIncSolar
          type(tKSP) :: ksp
          character(len=default_str_len) :: prefix
          real(ireals), pointer :: xv_b(:)
          real(ireals), pointer :: xg1d(:), xg4d(:, :, :, :)
          real(ireals), pointer :: x1d(:), x4d(:, :, :, :)

          xg1d => null()
          xg4d => null()
          x1d => null()
          x4d => null()

          if (solver%lopen_bc) then ! need to update side fluxes somewhere in the domain

            associate ( &
                & atm => solver%atm, &
                & sun => solver%sun, &
                & C_dir => solver%C_dir)

              call DMGetLocalVector(C_dir%da, local_incSolar, ierr); call CHKERR(ierr)
              call VecSet(local_incSolar, 0._ireals, ierr); call CHKERR(ierr)

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
                call MatSeqAIJSetPreallocation(A, C_dir%dof + i1, PETSC_NULL_INTEGER_ARRAY, ierr); call CHKERR(ierr)

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

                    call VecGetArrayRead(b, xv_b, ierr)
                    do k = C_dir%zs, C_dir%ze - 1
                      do src = 0, solver%dirside%dof - 1
                        ioff = solver%dirtop%dof + solver%dirside%dof + src
                        x4d(ioff, k, i, j + 1 - sun%yinc) = fac * xv_b(i1 + k * C_dir%dof + ioff)
                      end do
                    end do
                    call VecRestoreArrayRead(b, xv_b, ierr)
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

                    call VecGetArrayRead(b, xv_b, ierr)
                    do k = C_dir%zs, C_dir%ze - 1
                      do src = 0, solver%dirside%dof - 1
                        ioff = solver%dirtop%dof + src
                        x4d(ioff, k, i + 1 - sun%xinc, j) = fac * xv_b(i1 + k * C_dir%dof + ioff)
                      end do
                    end do
                    call VecRestoreArrayRead(b, xv_b, ierr)
                  end do
                end if

                call MatDestroy(A, ierr); call CHKERR(ierr)
                call VecDestroy(b, ierr); call CHKERR(ierr)
                call KSPDestroy(ksp, ierr); call CHKERR(ierr)
              end if ! have an outer domain boundary

              call restoreVecPointer(C_dir%da, local_incSolar, x1d, x4d)

              ! scatter local (with side BCs) back to plain incSolar array via temp tVec
              call DMCreateGlobalVector(C_dir%da, vIncSolar, ierr); call CHKERR(ierr)
              call VecSet(vIncSolar, 0._ireals, ierr); call CHKERR(ierr)
              call DMLocalToGlobalBegin(C_dir%da, local_incSolar, ADD_VALUES, vIncSolar, ierr); call CHKERR(ierr)
              call DMLocalToGlobalEnd(C_dir%da, local_incSolar, ADD_VALUES, vIncSolar, ierr); call CHKERR(ierr)
              call DMRestoreLocalVector(C_dir%da, local_incSolar, ierr); call CHKERR(ierr)

              call VecCopy(vIncSolar, incSolar, ierr); call CHKERR(ierr)
              call VecDestroy(vIncSolar, ierr); call CHKERR(ierr)

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
            call MatSetValue(A, irow, irow, -1._ireals, ADD_VALUES, ierr); call CHKERR(ierr)
          end do

          associate ( &
              & atm => solver%atm, &
              & sun => solver%sun)

            call VecGetArray(b, xv_b, ierr)
            xv_b(:) = 0
            do src = 1, solver%dirtop%dof
              xv_b(src) = -1._ireals
            end do
            call VecRestoreArray(b, xv_b, ierr)

            do k = C_dir%zs, C_dir%ze - 1
              ak = atmk(atm, k)

              if (atm%l1d(atmk(atm, k))) then

                do src = 0, solver%dirtop%dof - 1
                  irow = (k + 1) * C_dir%dof + src
                  icol = (k) * C_dir%dof + src

                  dtau = atm%kabs(ak, i, j) * atm%dz(ak, i, j) / sun%costheta
                  v = exp(-dtau)
                  call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
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
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do

                  do dst = 0, solver%dirside%dof - 1 ! top2x
                    ioff = solver%dirtop%dof + dst
                    v = cdir2dir(src, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do
                  do dst = 0, solver%dirside%dof - 1 ! top2y
                    ioff = solver%dirtop%dof + solver%dirside%dof + dst
                    v = cdir2dir(src, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do
                end do

                do src = 0, solver%dirside%dof - 1
                  ioffsrc = solver%dirtop%dof + src
                  icol = k * C_dir%dof + ioffsrc
                  do dst = 0, solver%dirtop%dof - 1 ! side2bot
                    v = cdir2dir(ioffsrc, dst)
                    irow = (k + 1) * C_dir%dof + dst
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do

                  do dst = 0, solver%dirside%dof - 1 ! side2x
                    ioff = solver%dirtop%dof + dst
                    v = cdir2dir(ioffsrc, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do
                  do dst = 0, solver%dirside%dof - 1 ! side2y
                    ioff = solver%dirtop%dof + solver%dirside%dof + dst
                    v = cdir2dir(ioffsrc, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do
                end do

                do src = 0, solver%dirside%dof - 1
                  ioffsrc = solver%dirtop%dof + solver%dirside%dof + src
                  icol = k * C_dir%dof + ioffsrc
                  do dst = 0, solver%dirtop%dof - 1 ! side2bot
                    v = cdir2dir(ioffsrc, dst)
                    irow = (k + 1) * C_dir%dof + dst
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do

                  do dst = 0, solver%dirside%dof - 1 ! side2x
                    ioff = solver%dirtop%dof + dst
                    v = cdir2dir(ioffsrc, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
                  end do
                  do dst = 0, solver%dirside%dof - 1 ! side2y
                    ioff = solver%dirtop%dof + solver%dirside%dof + dst
                    v = cdir2dir(ioffsrc, ioff)
                    irow = k * C_dir%dof + ioff
                    call MatSetValue(A, irow, icol, v, ADD_VALUES, ierr); call CHKERR(ierr)
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
#endif

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

99        continue ! update sample pts and go on to compute them
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

      !> Forward halo fill (interior → ghost). Mirrors DMGlobalToLocal.
      !> v is (1:dof, 1:zm, 1:gxm, 1:gym) in assumed-shape (1-based).
      !> Interior: 2..gxm-1 (x), 2..gym-1 (y). Ghosts: 1 and gxm (x), 1 and gym (y).
      subroutine halo_fill_5pt(comm, C, v, ierr)
        use mpi
        integer(mpiint), intent(in) :: comm
        type(t_coord), intent(in) :: C
        real(ireals), intent(inout) :: v(:, :, :, :)
        integer(mpiint), intent(out) :: ierr

        real(ireals), allocatable :: se(:, :, :), sw(:, :, :), sn(:, :, :), ss(:, :, :)
        real(ireals), allocatable :: re(:, :, :), rw(:, :, :), rn(:, :, :), rs(:, :, :)
        integer(mpiint) :: rqs(8), sts(MPI_STATUS_SIZE, 8)
        integer(mpiint) :: nw, ne, ns, nn
        integer(iintegers) :: gxm, gym

        nw = int(C%neighbors(10), mpiint)
        ne = int(C%neighbors(16), mpiint)
        ns = int(C%neighbors(4), mpiint)
        nn = int(C%neighbors(22), mpiint)

        gxm = size(v, 3, kind=iintegers)
        gym = size(v, 4, kind=iintegers)

        allocate (se(size(v, 1), size(v, 2), size(v, 4))); se = zero
        allocate (sw(size(v, 1), size(v, 2), size(v, 4))); sw = zero
        allocate (re(size(v, 1), size(v, 2), size(v, 4))); re = zero
        allocate (rw(size(v, 1), size(v, 2), size(v, 4))); rw = zero
        allocate (sn(size(v, 1), size(v, 2), size(v, 3) - 2)); sn = zero
        allocate (ss(size(v, 1), size(v, 2), size(v, 3) - 2)); ss = zero
        allocate (rn(size(v, 1), size(v, 2), size(v, 3) - 2)); rn = zero
        allocate (rs(size(v, 1), size(v, 2), size(v, 3) - 2)); rs = zero

        se = v(:, :, gxm - 1, :)
        sw = v(:, :, 2, :)
        sn = v(:, :, 2:gxm - 1, gym - 1)
        ss = v(:, :, 2:gxm - 1, 2)

        call MPI_Irecv(rw, size(rw, kind=mpiint), imp_ireals, nw, 16, comm, rqs(1), ierr)
        call MPI_Irecv(re, size(re, kind=mpiint), imp_ireals, ne, 10, comm, rqs(2), ierr)
        call MPI_Irecv(rs, size(rs, kind=mpiint), imp_ireals, ns, 22, comm, rqs(3), ierr)
        call MPI_Irecv(rn, size(rn, kind=mpiint), imp_ireals, nn, 4, comm, rqs(4), ierr)
        call MPI_Isend(se, size(se, kind=mpiint), imp_ireals, ne, 16, comm, rqs(5), ierr)
        call MPI_Isend(sw, size(sw, kind=mpiint), imp_ireals, nw, 10, comm, rqs(6), ierr)
        call MPI_Isend(sn, size(sn, kind=mpiint), imp_ireals, nn, 22, comm, rqs(7), ierr)
        call MPI_Isend(ss, size(ss, kind=mpiint), imp_ireals, ns, 4, comm, rqs(8), ierr)
        call MPI_Waitall(8_mpiint, rqs, sts, ierr)

        v(:, :, 1, :) = rw
        v(:, :, gxm, :) = re
        v(:, :, 2:gxm - 1, 1) = rs
        v(:, :, 2:gxm - 1, gym) = rn
      end subroutine

      !> Reverse halo reduce (ghost → interior, ADD_VALUES). Mirrors DMLocalToGlobal(ADD_VALUES).
      !> Sends ghost cell values to neighbors who add them onto their interior boundary cells,
      !> then zeros own ghost cells.
      subroutine halo_reduce_5pt(comm, C, v, ierr)
        use mpi
        integer(mpiint), intent(in) :: comm
        type(t_coord), intent(in) :: C
        real(ireals), intent(inout) :: v(:, :, :, :)
        integer(mpiint), intent(out) :: ierr

        real(ireals), allocatable :: se(:, :, :), sw(:, :, :), sn(:, :, :), ss(:, :, :)
        real(ireals), allocatable :: re(:, :, :), rw(:, :, :), rn(:, :, :), rs(:, :, :)
        integer(mpiint) :: rqs(8), sts(MPI_STATUS_SIZE, 8)
        integer(mpiint) :: nw, ne, ns, nn
        integer(iintegers) :: gxm, gym

        nw = int(C%neighbors(10), mpiint)
        ne = int(C%neighbors(16), mpiint)
        ns = int(C%neighbors(4), mpiint)
        nn = int(C%neighbors(22), mpiint)

        gxm = size(v, 3, kind=iintegers)
        gym = size(v, 4, kind=iintegers)

        allocate (se(size(v, 1), size(v, 2), size(v, 4))); se = zero
        allocate (sw(size(v, 1), size(v, 2), size(v, 4))); sw = zero
        allocate (re(size(v, 1), size(v, 2), size(v, 4))); re = zero
        allocate (rw(size(v, 1), size(v, 2), size(v, 4))); rw = zero
        allocate (sn(size(v, 1), size(v, 2), size(v, 3) - 2)); sn = zero
        allocate (ss(size(v, 1), size(v, 2), size(v, 3) - 2)); ss = zero
        allocate (rn(size(v, 1), size(v, 2), size(v, 3) - 2)); rn = zero
        allocate (rs(size(v, 1), size(v, 2), size(v, 3) - 2)); rs = zero

        ! Send ghost cells to owners; use same tag scheme as halo_fill_5pt.
        se = v(:, :, gxm, :)
        sw = v(:, :, 1, :)
        sn = v(:, :, 2:gxm - 1, gym)
        ss = v(:, :, 2:gxm - 1, 1)

        call MPI_Irecv(rw, size(rw, kind=mpiint), imp_ireals, nw, 16, comm, rqs(1), ierr)
        call MPI_Irecv(re, size(re, kind=mpiint), imp_ireals, ne, 10, comm, rqs(2), ierr)
        call MPI_Irecv(rs, size(rs, kind=mpiint), imp_ireals, ns, 22, comm, rqs(3), ierr)
        call MPI_Irecv(rn, size(rn, kind=mpiint), imp_ireals, nn, 4, comm, rqs(4), ierr)
        call MPI_Isend(se, size(se, kind=mpiint), imp_ireals, ne, 16, comm, rqs(5), ierr)
        call MPI_Isend(sw, size(sw, kind=mpiint), imp_ireals, nw, 10, comm, rqs(6), ierr)
        call MPI_Isend(sn, size(sn, kind=mpiint), imp_ireals, nn, 22, comm, rqs(7), ierr)
        call MPI_Isend(ss, size(ss, kind=mpiint), imp_ireals, ns, 4, comm, rqs(8), ierr)
        call MPI_Waitall(8_mpiint, rqs, sts, ierr)

        v(:, :, 2, :) = v(:, :, 2, :) + rw
        v(:, :, gxm - 1, :) = v(:, :, gxm - 1, :) + re
        v(:, :, 2:gxm - 1, 2) = v(:, :, 2:gxm - 1, 2) + rs
        v(:, :, 2:gxm - 1, gym - 1) = v(:, :, 2:gxm - 1, gym - 1) + rn

        v(:, :, 1, :) = zero
        v(:, :, gxm, :) = zero
        v(:, :, 2:gxm - 1, 1) = zero
        v(:, :, 2:gxm - 1, gym) = zero
      end subroutine

      end module
