module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options, ltwostr_only, lschwarzschild

  use m_helper_functions, only: CHKERR, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, strF2C, get_arg, &
    vec_proj_on_plane, cross_3d, norm, rotation_matrix_world_to_local_basis, &
    approx, swap, delta_scale, delta_scale_optprop, itoa, ftoa, &
    imp_allreduce_min, rotation_matrix_around_axis_vec

  use m_data_parameters, only : ireals, iintegers, mpiint, irealLUT, &
    i0, i1, i2, i3, i4, i5, i6, i7, i8, default_str_len, &
    zero, one, pi, EXP_MINVAL, EXP_MAXVAL

  use m_plex_grid, only: t_plexgrid, compute_face_geometry, &
    setup_edir_dmplex, setup_ediff_dmplex, setup_abso_dmplex, &
    orient_face_normals_along_sundir, compute_wedge_orientation, is_solar_src, get_inward_face_normal, &
    facevec2cellvec, icell_icon_2_plex, iface_top_icon_2_plex, get_consecutive_vertical_cell_idx, &
    get_top_bot_face_of_cell, destroy_plexgrid, determine_diff_incoming_outgoing_offsets, &
    TOAFACE, BOTFACE, SIDEFACE

  use m_optprop, only : t_optprop, t_optprop_wedge_5_8, OPP_1D_RETCODE
  use m_optprop_parameters, only : ldebug_optprop

  use m_schwarzschild, only: schwarzschild
  use m_twostream, only: delta_eddington_twostream

  use m_pprts_base, only : t_state_container, prepare_solution, destroy_solution, &
    t_solver_log_events, setup_log_events

  implicit none

  private

  public :: t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, &
    compute_face_geometry, set_plex_rt_optprop, destroy_plexrt_solver, &
    plexrt_get_result, scale_flx

  type t_plex_solver
    type(t_plexgrid), allocatable :: plex
    class(t_optprop), allocatable :: OPP

    type(tVec), allocatable :: kabs, ksca, g       ! in each cell [pStart..pEnd-1]
    type(tVec), allocatable :: albedo              ! on each surface face [defined on plex%srfc_boundary_dm]

    type(tVec), allocatable :: plck ! Planck Radiation in cell in each cell [W] [pStart..pEnd-1]
    ! srfc_emission is needed if plck is allocated:
    ! defined as planck(T_srfc) [plex%srfc_boundary_dm]
    type(tVec), allocatable :: srfc_emission

    type(t_state_container) :: solutions(-1000:1000)

    type(tVec), allocatable :: dirsrc, diffsrc
    type(tMat), allocatable :: Mdir
    type(tMat), allocatable :: Mdiff
    type(tKSP), allocatable :: kspdir
    type(tKSP), allocatable :: kspdiff

    type(tVec), allocatable :: dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2
    type(tVec), allocatable :: diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2
    type(tIS),  allocatable :: IS_diff_in_out_dof

    logical :: lenable_solutions_err_estimates=.True.

    type(t_solver_log_events) :: logs
  end type

  logical, parameter :: ldebug=.False.
  contains
    subroutine init_plex_rt_solver(plex, solver)
      type(t_plexgrid), intent(in) :: plex
      class(t_plex_solver), allocatable, intent(inout) :: solver

      call read_commandline_options(plex%comm)

      if(allocated(solver)) call CHKERR(1_mpiint, 'Should not call init_plex_rt_solver with already allocated solver object')
      allocate(solver)

      call setup_log_events(solver%logs, 'plexrt')

      allocate(solver%plex)
      solver%plex = plex

      allocate(t_optprop_wedge_5_8::solver%OPP)
      call solver%OPP%init(plex%comm)
    end subroutine

    subroutine destroy_plexrt_solver(solver, lfinalizepetsc)
      class(t_plex_solver), allocatable,  intent(inout) :: solver
      logical, intent(in) :: lfinalizepetsc

      integer(iintegers) :: uid
      integer(mpiint) :: ierr

      if(allocated(solver)) then
        if(allocated(solver%plex)) then
          call destroy_plexgrid(solver%plex)
          deallocate(solver%plex)
        endif

        if(allocated(solver%OPP)) then
          call solver%OPP%destroy()
          deallocate(solver%OPP)
        endif

        call dealloc_vec(solver%kabs)
        call dealloc_vec(solver%ksca)
        call dealloc_vec(solver%g   )
        call dealloc_vec(solver%albedo)
        call dealloc_vec(solver%plck)
        call dealloc_vec(solver%srfc_emission)
        call dealloc_vec(solver%dirsrc)
        call dealloc_vec(solver%diffsrc)

        call dealloc_vec(solver%dir_scalevec_Wm2_to_W)
        call dealloc_vec(solver%dir_scalevec_W_to_Wm2)
        call dealloc_vec(solver%diff_scalevec_Wm2_to_W)
        call dealloc_vec(solver%diff_scalevec_W_to_Wm2)

        if(allocated(solver%IS_diff_in_out_dof)) then
          call ISDestroy(solver%IS_diff_in_out_dof, ierr); call CHKERR(ierr)
          deallocate(solver%IS_diff_in_out_dof)
        endif

        if(allocated(solver%Mdir)) then
          call MatDestroy(solver%Mdir, ierr); call CHKERR(ierr)
          deallocate(solver%Mdir)
        endif

        if(allocated(solver%Mdiff)) then
          call MatDestroy(solver%Mdiff, ierr); call CHKERR(ierr)
          deallocate(solver%Mdiff)
        endif
        if(allocated(solver%kspdir)) then
          call KSPDestroy(solver%kspdir, ierr); call CHKERR(ierr)
          deallocate(solver%kspdir)
        endif
        if(allocated(solver%kspdiff)) then
          call KSPDestroy(solver%kspdiff, ierr); call CHKERR(ierr)
          deallocate(solver%kspdiff)
        endif

        do uid=lbound(solver%solutions,1),ubound(solver%solutions,1)
          call destroy_solution(solver%solutions(uid))
        enddo
        deallocate(solver)
      endif

      if(lfinalizepetsc) then
        call PetscFinalize(ierr) ;call CHKERR(ierr)
      endif
    contains
      subroutine dealloc_vec(v)
        type(tVec), allocatable, intent(inout) :: v
        if(allocated(v)) then
          call VecDestroy(v, ierr); call CHKERR(ierr)
          deallocate(v)
        endif
      end subroutine
    end subroutine

    subroutine set_plex_rt_optprop(solver, vlwc, viwc, vert_integrated_kabs, vert_integrated_ksca)
      use m_helper_functions, only : delta_scale
      class(t_plex_solver), allocatable, intent(inout) :: solver
      type(tVec),intent(in), optional :: vlwc, viwc
      real(ireals), optional :: vert_integrated_kabs, vert_integrated_ksca
      real(ireals), pointer :: xlwc(:), xiwc(:)
      real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:)
      integer(iintegers) :: geom_offset
      real(ireals) :: dz

      real(ireals),parameter :: w0 = .99, reff_w=10, reff_i=20, rayleigh=1e-4
      real(ireals) :: ksca_tot, kabs_tot, g_tot
      real(ireals) :: ksca_cld, kabs_cld, g_cld

      integer(iintegers) :: i, cStart, cEnd
      integer(mpiint) :: myid, ierr

      if(.not.allocated(solver)) call CHKERR(1_mpiint, 'set_plex_rt_optprop::solver has to be allocated')
      if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'set_plex_rt_optprop::solver%plex has to be allocated')
      if(.not.allocated(solver%plex%cell1_dm)) call CHKERR(1_mpiint, 'cell1_dm should be allocated already')
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(solver%plex%dm, i3, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

      if(.not.allocated(solver%kabs)) then
        allocate(solver%kabs)
        call DMCreateGlobalVector(solver%plex%cell1_dm, solver%kabs, ierr); call CHKERR(ierr)
      endif
      if(.not.allocated(solver%ksca)) then
        allocate(solver%ksca)
        call DMCreateGlobalVector(solver%plex%cell1_dm, solver%ksca, ierr); call CHKERR(ierr)
      endif
      if(.not.allocated(solver%g   )) then
        allocate(solver%g)
        call DMCreateGlobalVector(solver%plex%cell1_dm, solver%g, ierr); call CHKERR(ierr)
      endif

      if(present(vlwc)) then
        call VecGetArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      endif
      if(present(viwc)) then
        call VecGetArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)
      endif

      call DMGetSection(solver%plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call VecGetArrayF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solver%g   , xg   , ierr); call CHKERR(ierr)
      do i = cStart, cEnd-1
        call PetscSectionGetFieldOffset(geomSection, i, i3, geom_offset, ierr); call CHKERR(ierr)
        dz = geoms(i1+geom_offset)

        kabs_tot = get_arg(rayleigh*(one-w0)*dz, vert_integrated_kabs/solver%plex%Nlay/dz)
        ksca_tot = get_arg(rayleigh*     w0 *dz, vert_integrated_ksca/solver%plex%Nlay/dz)

        g_tot    = zero

        if(present(vlwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xlwc(i+1) * 1e3) / reff_w * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xlwc(i+1) * 1e3) / reff_w * (w0)
          g_cld    = .5_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = g_cld * ksca_cld / (ksca_cld + ksca_tot)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        if(present(viwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xiwc(i+1) * 1e3) / reff_i * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xiwc(i+1) * 1e3) / reff_i * (w0)
          g_cld    = .3_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = (g_tot * ksca_tot + g_cld * ksca_cld) / (ksca_tot + ksca_cld)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        xkabs(i1+i) = kabs_tot
        xksca(i1+i) = ksca_tot
        xg(i1+i)    = g_tot
        call delta_scale(xkabs(i1+i), xksca(i1+i), xg(i1+i))
      enddo

      if(present(vlwc)) then
        call VecRestoreArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      endif
      if(present(viwc)) then
        call VecRestoreArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)
      endif

      print *,'Min/Max of kabs', minval(xkabs), maxval(xkabs)
      print *,'Min/Max of ksca', minval(xksca), maxval(xksca)
      print *,'Min/Max of g   ', minval(xg   ), maxval(xg   )

      call VecRestoreArrayF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(solver%g   , xg   , ierr); call CHKERR(ierr)
    end subroutine

    subroutine get_in_out_dof_offsets(IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, opt_xv)
      type(tIS), intent(in) :: IS_diff_in_out_dof
      integer(iintegers), intent(in) :: icell
      integer(iintegers), allocatable, intent(inout) :: incoming_offsets(:), outgoing_offsets(:)
      integer(iintegers), pointer, intent(in), optional :: opt_xv(:)
      integer(iintegers), pointer :: xv(:)
      integer(iintegers) :: bs, idx_offset
      integer(mpiint) :: ierr

      if(present(opt_xv)) then
        bs = opt_xv(size(opt_xv))
        if(.not.allocated(incoming_offsets)) allocate(incoming_offsets(bs))
        if(.not.allocated(outgoing_offsets)) allocate(outgoing_offsets(bs))

        idx_offset = icell*bs
        incoming_offsets = opt_xv(i1+idx_offset:idx_offset+bs/2)
        outgoing_offsets = opt_xv(i1+idx_offset+bs/2:idx_offset+bs)
      else
        call ISGetIndicesF90(IS_diff_in_out_dof, xv, ierr); call CHKERR(ierr)
        bs = xv(size(xv))
        if(.not.allocated(incoming_offsets)) allocate(incoming_offsets(bs))
        if(.not.allocated(outgoing_offsets)) allocate(outgoing_offsets(bs))

        idx_offset = icell*bs
        incoming_offsets = xv(i1+idx_offset:idx_offset+bs/2)
        outgoing_offsets = xv(i1+idx_offset+bs/2:idx_offset+bs)
        call ISRestoreIndicesF90(IS_diff_in_out_dof, xv, ierr); call CHKERR(ierr)
      endif
    end subroutine

    subroutine setup_IS_diff_in_out_dof(plex, ediffdm, IS_diff_in_out_dof)
      type(t_plexgrid), allocatable, intent(in) :: plex
      type(tDM), allocatable, intent(in) :: ediffdm
      type(tIS), allocatable, intent(inout) :: IS_diff_in_out_dof

      integer(iintegers) :: icell, bs, N, idx_offset
      integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:), idx(:)
      integer(mpiint) :: comm, myid, ierr

      integer(iintegers),pointer :: faces_of_cell(:)
      type(tPetscSection) :: ediffSection
      integer(iintegers) :: ndof, num_dof, max_num_dof, iface

      if(allocated(IS_diff_in_out_dof)) call CHKERR(1_mpiint, 'IS_diff_in_out_dof already allocated')
      if(.not.allocated(ediffdm)) call CHKERR(1_mpiint, 'ediffdm has to be allocated')
      if(.not.allocated(plex)) call CHKERR(1_mpiint, 'plex has to be allocated')

      call PetscObjectGetComm(ediffdm, comm, ierr); call CHKERR(ierr)
      call DMGetSection(ediffdm, ediffSection, ierr); call CHKERR(ierr)

      ! Compute maximum number of dof per cell
      max_num_dof = 0
      do icell = plex%cStart, plex%cEnd-1
        call DMPlexGetCone(ediffdm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
        num_dof = 0
        do iface = 1, size(faces_of_cell)
          call PetscSectionGetDof(ediffSection, faces_of_cell(iface), ndof, ierr); call CHKERR(ierr)
          num_dof = num_dof + ndof
        enddo
        call DMPlexRestoreCone(ediffdm, icell, faces_of_cell,ierr); call CHKERR(ierr)
        max_num_dof = max(max_num_dof, num_dof)
      enddo
      bs = max_num_dof

      N = (plex%cEnd - plex%cStart) * bs +1 ! plus one for the blocksize, we will write that to the end of the IS
      allocate(idx(N))

      idx(N) = bs

      do icell = plex%cStart, plex%cEnd-1
        call determine_diff_incoming_outgoing_offsets(plex, ediffdm, icell, incoming_offsets, outgoing_offsets)
        idx_offset = icell*bs
        do ndof = 1, size(incoming_offsets)
          idx(idx_offset+ndof)      = incoming_offsets(ndof)
          idx(idx_offset+bs/2+ndof) = outgoing_offsets(ndof)
        enddo
        idx(idx_offset+ndof     :idx_offset+bs/2) = -i1
        idx(idx_offset+bs/2+ndof:idx_offset+bs  ) = -i1
        !print *,bs, 'icell', icell, 'incoming', incoming_offsets, 'out', outgoing_offsets, &
        !            'idx_offset',idx(i1+idx_offset:idx_offset+bs)
      enddo

      allocate(IS_diff_in_out_dof)
      call ISCreateGeneral(comm, N, idx, PETSC_COPY_VALUES, IS_diff_in_out_dof, ierr); call CHKERR(ierr)
      deallocate(idx)
      if(ldebug) then
        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        if(myid.eq.0) print *,'IS_diff_in_out_dof blocksize = '//itoa(bs)
      endif
    end subroutine

    subroutine run_plex_rt_solver(solver, lthermal, lsolar, sundir, opt_solution_uid, opt_solution_time)
      class(t_plex_solver), allocatable, intent(inout) :: solver
      logical, intent(in) :: lthermal, lsolar
      real(ireals), intent(in) :: sundir(3) ! cartesian direction of sun rays, norm of vector is the energy in W/m2
      integer(iintegers), intent(in), optional :: opt_solution_uid
      real(ireals), intent(in), optional :: opt_solution_time

      integer(iintegers) :: suid

      integer(mpiint) :: myid, ierr

      real(ireals), save :: last_sundir(3) = [zero,zero,zero]

      if(.not.allocated(solver)) call CHKERR(1_mpiint, 'run_plex_rt_solver::solver has to be allocated')

      if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'run_plex_rt_solver::plex has to be allocated first')
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(solver%plex%geom_dm)) call CHKERR(1_mpiint, 'run_plex_rt_solver::geom_dm has to be allocated first')
      if(.not.allocated(solver%plex%srfc_boundary_dm)) call CHKERR(1_mpiint, 'run_plex_rt_solver::srfc_boundary_dm has to be allocated first')
      if(.not.allocated(solver%kabs  )) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, kabs, ksca, g have to be allocated first')
      if(.not.allocated(solver%ksca  )) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, kabs, ksca, g have to be allocated first')
      if(.not.allocated(solver%g     )) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, kabs, ksca, g have to be allocated first')
      if(.not.allocated(solver%albedo)) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, albedo has to be allocated first')
      if(lthermal) then
        if(.not.allocated(solver%plck)) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, planck radiation vec has to be allocated first')
        if(.not.allocated(solver%srfc_emission)) call CHKERR(1_mpiint, 'run_plex_rt_solver::optprop, srfc emission vec has to be allocated first')
      endif

      call PetscLogStagePush(solver%logs%stage_solve, ierr); call CHKERR(ierr)

      if(.not.allocated(solver%plex%geom_dm))  call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      if(.not.allocated(solver%plex%edir_dm))  call setup_edir_dmplex(solver%plex, solver%plex%dm, solver%plex%edir_dm)
      if(.not.allocated(solver%plex%ediff_dm)) call setup_ediff_dmplex(solver%plex, solver%plex%dm, solver%plex%ediff_dm)
      if(.not.allocated(solver%plex%abso_dm))  call setup_abso_dmplex(solver%plex%dm, solver%plex%abso_dm)

      if(.not.allocated(solver%IS_diff_in_out_dof)) &
        call setup_IS_diff_in_out_dof(solver%plex, solver%plex%ediff_dm, solver%IS_diff_in_out_dof)

      ! Prepare the space for the solution
      suid = get_arg(i0, opt_solution_uid)

      if(.not.solver%solutions(suid)%lset) then
        call prepare_solution(solver%plex%edir_dm, solver%plex%ediff_dm, solver%plex%abso_dm, &
          lsolar=lsolar, solution=solver%solutions(suid), uid=suid)
      endif

      call dump_optical_properties(solver%kabs, solver%ksca, solver%g, solver%albedo, &
        plck=solver%plck, srfc_emission=solver%srfc_emission, postfix='_'//itoa(suid))

      ! Wedge Orientation is used in solar and thermal case alike
      if(any(last_sundir.ne.sundir).or.&                       ! update wedge orientations if sundir has changed
        .not.allocated(solver%plex%wedge_orientation_dm)) then ! or if we lost the info somehow... e.g. happens after destroy_solver
        call PetscLogEventBegin(solver%logs%orient_face_normals, ierr)
        call orient_face_normals_along_sundir(solver%plex, sundir)
        call PetscLogEventEnd(solver%logs%orient_face_normals, ierr)

        call PetscLogEventBegin(solver%logs%compute_orientation, ierr)
        call compute_wedge_orientation(solver%plex, sundir, solver%plex%wedge_orientation_dm, &
          solver%plex%wedge_orientation)
        call PetscLogEventEnd(solver%logs%compute_orientation, ierr)

        last_sundir = sundir
      endif

      associate( solution => solver%solutions(suid) )

        ! --------- Calculate 1D Radiative Transfer ------------
        if(ltwostr_only) then
          call PetscLogEventBegin(solver%logs%solve_twostream, ierr)
          if( (solution%lsolar_rad.eqv..False.) .and. lschwarzschild ) then
            call schwarz(solver, solution)
          else
            call twostream(solver%plex, solver%kabs, solver%ksca, solver%g, &
              solver%albedo, sundir, solution, plck=solver%plck, srfc_emission=solver%srfc_emission)
          endif
          call PetscLogEventEnd(solver%logs%solve_twostream, ierr)

          ! if(ldebug) print *,'1D calculation done', suid
          goto 99
        endif

        if(solution%lsolar_rad) then
          ! Output of wedge_orient vec
          call PetscLogEventBegin(solver%logs%setup_dir_src, ierr)
          call create_edir_src_vec(solver%plex, solver%plex%edir_dm, norm(sundir), &
            solver%kabs, solver%ksca, solver%g, &
            sundir/norm(sundir), solver%dirsrc)
          call PetscLogEventEnd(solver%logs%setup_dir_src, ierr)

          ! Output of srcVec
          if(ldebug) then
            call PetscLogEventBegin(solver%logs%debug_output, ierr)
            call scale_facevec(solver%plex, solver%plex%edir_dm, solver%dirsrc, lW_to_Wm2=.True.)
            call facevec2cellvec(solver%plex%dm, solver%plex%edir_dm, solver%dirsrc)
            call scale_facevec(solver%plex, solver%plex%edir_dm, solver%dirsrc, lW_to_Wm2=.False.)
            call PetscLogEventEnd(solver%logs%debug_output, ierr)
          endif

          ! Create Direct Matrix
          call PetscLogEventBegin(solver%logs%setup_Mdir, ierr)
          call create_edir_mat(solver, solver%plex, solver%OPP, solver%kabs, solver%ksca, solver%g, sundir, solver%Mdir)
          call PetscLogEventEnd(solver%logs%setup_Mdir, ierr)

          call scale_flx(solver%plex, &
            solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
            solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
            solution, lWm2=.False., logevent=solver%logs%scale_flx)

          ! Solve Direct Matrix
          call PetscLogEventBegin(solver%logs%solve_Mdir, ierr)
          call solve_plex_rt(solver%plex%edir_dm, solver%dirsrc, solver%Mdir, solver%kspdir, solution%edir, &
            ksp_residual_history=solution%dir_ksp_residual_history, prefix='dir_')
          call PetscLogEventEnd(solver%logs%solve_Mdir, ierr)
          call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, &
            '-show_edir_vec_global', ierr); call CHKERR(ierr)
          solution%lWm2_dir = .False.
          solution%lchanged = .True.

          ! Output of Edir Vec
          if(ldebug) then
            call PetscLogEventBegin(solver%logs%debug_output, ierr)
            call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.True.)
            call facevec2cellvec(solver%plex%dm, solver%plex%edir_dm, solution%edir)
            call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.False.)
            call PetscLogEventEnd(solver%logs%debug_output, ierr)
          endif
        endif

        ! Create Diffuse Src
        call PetscLogEventBegin(solver%logs%setup_diff_src, ierr)
        call create_ediff_src_vec(solver, solver%plex, solver%OPP, solver%plex%ediff_dm, &
          solver%kabs, solver%ksca, solver%g, solver%albedo, &
          lthermal, lsolar, solver%diffsrc, &
          solver%plck, solver%srfc_emission, &
          solver%plex%edir_dm, solution%edir)
        call PetscLogEventEnd(solver%logs%setup_diff_src, ierr)

        ! Output of Diffuse Src Vec
        if(ldebug) then
          call PetscLogEventBegin(solver%logs%debug_output, ierr)
          call scale_facevec(solver%plex, solver%plex%ediff_dm, solver%diffsrc, lW_to_Wm2=.True.)
          call facevec2cellvec(solver%plex%dm, solver%plex%ediff_dm, solver%diffsrc)
          call scale_facevec(solver%plex, solver%plex%ediff_dm, solver%diffsrc, lW_to_Wm2=.False.)
          call PetscLogEventEnd(solver%logs%debug_output, ierr)
        endif

        ! Create Diffuse Matrix
        call PetscLogEventBegin(solver%logs%setup_Mdiff, ierr)
        call create_ediff_mat(solver, solver%plex, solver%OPP, &
          solver%kabs, solver%ksca, solver%g, solver%albedo, solver%Mdiff)
        call PetscObjectViewFromOptions(solver%Mdiff, PETSC_NULL_MAT, &
          '-show_Mediff_'//itoa(solution%uid), ierr); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%setup_Mdiff, ierr)

        call scale_flx(solver%plex, &
          solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
          solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
          solution, lWm2=.False., logevent=solver%logs%scale_flx)

        ! Solve Diffuse Matrix
        call PetscLogEventBegin(solver%logs%solve_Mdiff, ierr)
        call solve_plex_rt(solver%plex%ediff_dm, solver%diffsrc, solver%Mdiff, solver%kspdiff, solution%ediff, &
          ksp_residual_history=solution%diff_ksp_residual_history, prefix='diff_')
        call PetscLogEventEnd(solver%logs%solve_Mdiff, ierr)
        call PetscObjectSetName(solution%ediff, 'ediff', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, &
          '-show_ediff_vec_global', ierr); call CHKERR(ierr)
        solution%lWm2_diff = .False.
        solution%lchanged = .True.

        99 continue ! this is the quick exit final call where we clean up before the end of the routine

        ! Bring solution into a coherent state, i.e. update absorption etc.
        call restore_solution(solver, solution, opt_solution_time)

      end associate

      call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop solver%logs%stage_solve
    end subroutine

    subroutine dump_optical_properties(kabs, ksca, g, albedo, plck, srfc_emission, postfix)
      character(len=*), intent(in) :: postfix
      type(tVec), allocatable, intent(in) :: kabs, ksca, g, albedo, plck, srfc_emission
      call dump_var(kabs, 'kabs')
      call dump_var(ksca, 'ksca')
      call dump_var(g, 'g')
      call dump_var(plck, 'plck')
      call dump_var(srfc_emission, 'srfc_emission')
      call dump_var(albedo, 'albedo')
    contains
      subroutine dump_var(var, varname)
        type(tVec), allocatable, intent(in) :: var
        character(len=*), intent(in) :: varname
        character(len=default_str_len) :: oldname
        integer(mpiint) :: ierr
        if(.not.allocated(var)) return
        call PetscObjectGetName(var, oldname, ierr); call CHKERR(ierr)
        call PetscObjectSetName(var, trim(varname)//trim(postfix), ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(var, PETSC_NULL_VEC, '-dump_optprop_'//trim(varname), ierr); call CHKERR(ierr)
        call PetscObjectSetName(var, trim(oldname), ierr); call CHKERR(ierr)
      end subroutine
    end subroutine

    subroutine create_edir_src_vec(plex, edirdm, E0, kabs, ksca, g, sundir, srcVec)
      type(t_plexgrid), allocatable, intent(in) :: plex
      type(tDM),allocatable, intent(in) :: edirdm
      real(ireals), intent(in) :: E0, sundir(3)
      type(tVec), allocatable, intent(in) :: kabs, ksca, g
      type(tVec), allocatable, intent(inout) :: srcVec

      integer(iintegers) :: i, voff
      type(tPetscSection) :: edirsection

      type(tVec) :: localVec
      real(ireals), pointer :: xv(:)

      type(tIS) :: boundary_ids
      integer(iintegers) :: iface, icell

      integer(iintegers), pointer :: xx_v(:)
      integer(iintegers), pointer :: cell_support(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: area, mu, face_normal(3)

      integer(mpiint) :: myid, ierr

      if(.not.allocated(plex)) stop 'called create_src_vec but plex is not allocated'
      call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(edirdm)) call CHKERR(myid+1, 'called create_src_vec but edirdm is not allocated')

      if(.not.allocated(plex%geom_dm)) call CHKERR(myid+1, 'get_normal_of_first_TOA_face::needs allocated geom_dm first')

      call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetSection(edirdm, edirsection, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(edirsection, PETSC_NULL_SECTION, '-show_src_section', ierr); call CHKERR(ierr)

      ! Now lets get vectors!
      if(.not.allocated(srcVec)) then
        allocate(srcVec)
        call DMCreateGlobalVector(edirdm, srcVec,ierr); call CHKERR(ierr)
        call PetscObjectSetName(srcVec, 'srcVecGlobal', ierr);call CHKERR(ierr)
      endif
      call VecSet(srcVec, zero, ierr); call CHKERR(ierr)


      call DMGetLocalVector(edirdm, localVec,ierr); call CHKERR(ierr)
      call VecSet(localVec, zero, ierr); call CHKERR(ierr)

      call VecGetArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      call DMGetStratumIS(edirdm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
      if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
      else
        call ISGetIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)

        ! First set the TOA boundary fluxes on faces
        do i = 1, size(xx_v)
          iface = xx_v(i)
          call DMPlexGetSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

          call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
          area = geoms(i1+geom_offset)

          call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)

          if(is_solar_src(face_normal, sundir)) then
            mu = dot_product(sundir, face_normal)
            call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
            xv(voff+1) = E0 * area * mu
          endif
        enddo
        call ISRestoreIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)
      endif ! TOA boundary ids

      call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      !call set_sidewards_direct_fluxes()

      call DMLocalToGlobalBegin(edirdm, localVec, INSERT_VALUES, srcVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (edirdm, localVec, INSERT_VALUES, srcVec, ierr); call CHKERR(ierr)
      call PetscObjectSetName(srcVec, 'srcVec', ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(srcVec, PETSC_NULL_VEC, '-show_src_vec_global', ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(edirdm, localVec, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      contains
        subroutine set_sidewards_direct_fluxes()
          type(tIS) :: IS_side_faces
          integer(iintegers), pointer :: iside_faces(:), cell_support(:), faces_of_cell(:)
          integer(iintegers), allocatable :: vert_cell_idx(:)
          integer(iintegers) :: topface, botface, iface_side, k, o, topoffset, botoffset
          real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
          real(ireals) :: side_face_normal(3)
          real(ireals) :: tau, dz, transport_this_cell, mu_side

          type(tVec) :: vedir
          real(ireals), pointer :: xedir(:)

          call DMGetStratumIS(plex%edir_dm, 'DomainBoundary', SIDEFACE, IS_side_faces, ierr); call CHKERR(ierr)
          if (IS_side_faces.eq.PETSC_NULL_IS) then ! dont have side boundary faces
          else

            call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
            call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
            call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)


            ! First we compute the edir vec on all columns that are on the domain side and the sun is shining on it
            call DMGetGlobalVector(edirdm, vedir, ierr); call CHKERR(ierr)
            call VecSet(vedir, -one, ierr); call CHKERR(ierr)
            call VecGetArrayF90(vedir, xedir, ierr); call CHKERR(ierr)

            call ISGetIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)

            do o = 1, size(iside_faces)
              iface_side = iside_faces(o)

              call DMPlexGetSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)
              icell = cell_support(1)
              call DMPlexRestoreSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)

              k = plex%zindex(icell)
              if(k.ne.1) cycle ! only start if side face, i.e. cell is close to TOA

              call get_inward_face_normal(iface_side, icell, geomSection, geoms, side_face_normal)
              if(.not.is_solar_src(side_face_normal, sundir)) cycle ! if sun is not shining on this side face, skip this column

                call get_consecutive_vertical_cell_idx(plex, icell, vert_cell_idx)

                call DMPlexGetCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
                topface = faces_of_cell(1)
                botface = faces_of_cell(2)
                call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)

                call get_inward_face_normal(topface, icell, geomSection, geoms, face_normal)
                mu = dot_product(sundir, face_normal)
                if(.not.is_solar_src(face_normal, sundir)) cycle ! if sun is not shining on top of TOA face, skip it

                ! compute direct radiation in 1D fashion
                call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)
                xedir(i1+topoffset) = E0
                print *,icell,'edir', 1, xedir(i1+topoffset)

                do k = 1, plex%Nlay
                  icell = vert_cell_idx(k)

                  call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
                  dz = geoms(i1+geom_offset)

                  tau = (xkabs((i1+icell)) + xksca((i1+icell)))

                  call DMPlexGetCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
                  topface = faces_of_cell(1)
                  botface = faces_of_cell(2)
                  call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)

                  call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)
                  call PetscSectionGetOffset(edirsection, botface, botoffset, ierr); call CHKERR(ierr)

                  xedir(i1+botoffset) = xedir(i1+topoffset) * exp(-tau*dz/mu)

                  print *,icell,'edir', k, xedir(i1+botoffset)
                enddo
            enddo

            ! Now we go ahead and try to estimate the fluxes along the side boundary

            call VecGetArrayF90(localVec, xv, ierr); call CHKERR(ierr)


            do o = 1, size(iside_faces)
              iface_side = iside_faces(o)

              call DMPlexGetSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)
              icell = cell_support(1)
              call DMPlexRestoreSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)

              call get_inward_face_normal(iface_side, icell, geomSection, geoms, side_face_normal)
              if(.not.is_solar_src(side_face_normal, sundir)) cycle ! if sun is not shining on this side face, skip this

              call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
              dz = geoms(i1+geom_offset)

              tau = (xkabs((i1+icell)) + xksca((i1+icell)))

              call DMPlexGetCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
              topface = faces_of_cell(1)
              botface = faces_of_cell(2)
              call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)

              call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)

              call PetscSectionGetFieldOffset(geomSection, iface_side, i2, geom_offset, ierr); call CHKERR(ierr)
              area = geoms(i1+geom_offset)

              call get_inward_face_normal(topface, icell, geomSection, geoms, face_normal)
              mu = dot_product(sundir, face_normal)
              transport_this_cell = exp(-tau*dz/mu)

              mu_side = dot_product(sundir, side_face_normal)

              call PetscSectionGetOffset(edirsection, iface_side, voff, ierr); call CHKERR(ierr)
              xv(i1+voff) = xedir(i1+topoffset) * area * mu * transport_this_cell
              print *,icell,'edir top', xedir(i1+topoffset), 'sidesrc', xv(i1+voff)
            enddo

            call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

            call ISRestoreIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)

            call VecRestoreArrayF90(vedir, xedir, ierr); call CHKERR(ierr)

            call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
            call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
            call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
            !call CHKERR(1_mpiint, 'DEBUG')
          endif ! TOA boundary ids
      end subroutine
    end subroutine

    !> @brief setup source term for diffuse radiation
    !> @details this is either direct radiation scattered into one of the diffuse coeffs:
    !> \n direct source term is
    !> \n   direct radiation times the dir2diff coeffs
    !> \n or it may be that we have a source term due to thermal emission --
    !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
    !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
    subroutine create_ediff_src_vec(solver, plex, OPP, ediffdm, kabs, ksca, g, albedo, &
        lthermal, lsolar, srcVec, &
        plckVec, srfc_emission, edirdm, edirVec)
      class(t_plex_solver), allocatable, intent(in) :: solver
      type(t_plexgrid), intent(in) :: plex
      class(t_optprop), intent(in) :: OPP
      type(tDM), allocatable, intent(in) :: ediffdm
      type(tVec), allocatable, intent(in) :: kabs, ksca, g ! cell1_dm
      type(tVec), allocatable, intent(in) :: albedo ! srfc_boundary_dm
      logical, intent(in) :: lthermal, lsolar
      type(tVec), allocatable, intent(inout) :: srcVec

      type(tVec), allocatable, intent(in), optional :: plckVec ! cell1_dm
      type(tVec), allocatable, intent(in), optional :: srfc_emission  ! srfc_boundary_dm

      type(tDM), allocatable, intent(in), optional :: edirdm
      type(tVec), allocatable, intent(in), optional :: edirVec

      type(tVec) :: lsrcVec

      type(tPetscSection) :: geomSection, wedgeSection, ediffSection, srfcSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec

      real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

      real(ireals), pointer :: xb(:)

      integer(mpiint) :: ierr

      if(.not.allocated(ediffdm)) call CHKERR(1_mpiint, 'ediff dm has to be allocated before')
      if(.not.allocated(kabs  )) call CHKERR(1_mpiint, 'kabs   has to be allocated')
      if(.not.allocated(ksca  )) call CHKERR(1_mpiint, 'ksca   has to be allocated')
      if(.not.allocated(g     )) call CHKERR(1_mpiint, 'g      has to be allocated')
      if(.not.allocated(albedo)) call CHKERR(1_mpiint, 'albedo has to be allocated')

      if(.not.allocated(srcVec)) then
        allocate(srcVec)
        call DMCreateGlobalVector(ediffdm, srcVec, ierr); call CHKERR(ierr)
        call PetscObjectSetName(srcVec, 'DiffSrcVec', ierr);call CHKERR(ierr)
      endif
      call DMGetSection(ediffdm, ediffSection, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
      call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
      call DMGetSection(plex%srfc_boundary_dm, srfcSection, ierr); CHKERRQ(ierr)

      call VecGetArrayReadF90(kabs  , xkabs  , ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(ksca  , xksca  , ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(g     , xg     , ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

      call DMGetLocalVector(ediffdm, lsrcVec, ierr); call CHKERR(ierr)
      call VecSet(lsrcVec, zero, ierr); call CHKERR(ierr)
      call VecGetArrayF90(lsrcVec, xb, ierr); call CHKERR(ierr)

      if(lsolar) call set_solar_source()
      if(lthermal) call set_thermal_source()

      call VecRestoreArrayF90(lsrcVec, xb, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call VecRestoreArrayReadF90(kabs  , xkabs  , ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(ksca  , xksca  , ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(g     , xg     , ierr); call CHKERR(ierr)

      call VecSet(srcVec, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(ediffdm, lsrcVec, ADD_VALUES, srcVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (ediffdm, lsrcVec, ADD_VALUES, srcVec, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(ediffdm, lsrcVec, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(srcVec, PETSC_NULL_VEC, '-show_diff_src_vec', ierr); call CHKERR(ierr)

      contains
        subroutine set_solar_source()
          type(tVec) :: ledirVec
          type(tPetscSection) :: edirSection

          integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

          integer(iintegers), pointer :: faces_of_cell(:)
          integer(iintegers) :: i, icell, iface, isrc, idst, numDof

          integer(iintegers) :: face_plex2bmc(5), wedge_offset, geom_offset, diff_plex2bmc(8)
          real(ireals), pointer :: xedir(:)
          real(ireals) :: zenith, azimuth, area_bot, area_top
          real(ireals) :: dir2diff(8)
          logical :: lsrc(5)

          real(ireals) :: dz, coeff(5*8)
          integer(iintegers), pointer :: xinoutdof(:)

          integer(mpiint) :: myid

          if(.not.allocated(edirVec)) return

          call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
          if(any([present(edirdm),present(edirVec)]).and. &
            .not.all([present(edirdm),present(edirVec)])) &
            call CHKERR(1_mpiint, 'either provide all vars for direct radiation or none')

          call DMGetSection(edirdm, edirSection, ierr); call CHKERR(ierr)
          call DMGetLocalVector(edirdm, ledirVec, ierr); call CHKERR(ierr)

          call DMGlobalToLocalBegin(edirdm, edirVec, INSERT_VALUES, ledirVec, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd  (edirdm, edirVec, INSERT_VALUES, ledirVec, ierr); call CHKERR(ierr)


          call ISGetIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
          call VecGetArrayReadF90(ledirVec, xedir, ierr); call CHKERR(ierr)

          do icell = plex%cStart, plex%cEnd-1
            call DMPlexGetCone(ediffdm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

            call get_in_out_dof_offsets(solver%IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)
            !print *,myid,'icell', icell, 'faces_of_cell', faces_of_cell
            !print *,myid,'icell', icell, 'incoming_offsets', incoming_offsets
            !print *,myid,'icell', icell, 'outgoing_offsets', outgoing_offsets

            call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)
            zenith  = wedgeorient(wedge_offset+i1)
            azimuth = wedgeorient(wedge_offset+i2)

            do iface = 1, size(faces_of_cell)
              face_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
              lsrc(iface) = nint(wedgeorient(wedge_offset+i8+iface), iintegers) .eq. i1
            enddo
            call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            call PetscLogEventBegin(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff, ierr, &
              angles=[real(rad2deg(azimuth), irealLUT), real(rad2deg(zenith), irealLUT)])

            if(ierr.eq.OPP_1D_RETCODE) then
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
              area_top = geoms(i1+geom_offset)
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
              area_bot = geoms(i1+geom_offset)
              coeff(4*8+1) = coeff(4*8+1) * area_bot / area_top
            endif

            call PetscLogEventEnd(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)

            do i = 1, size(faces_of_cell)
              iface = faces_of_cell(i)
              !print *,'icell '//itoa(icell)//' i '//itoa(i)//' face ',iface, ':plex2bmc', face_plex2bmc
              call PetscSectionGetFieldDof(edirSection, iface, i0, numDof, ierr); call CHKERR(ierr)
              if(numDof.eq.i0) cycle ! dont have dofs on this incoming face

              if(lsrc(i)) then
                dir2diff = coeff(face_plex2bmc(i):size(coeff):i5) ! dir2diff in src ordering
                call PetscSectionGetOffset(edirSection, iface, isrc, ierr); call CHKERR(ierr)

                do idst = 1, size(outgoing_offsets)
                  if(dir2diff(diff_plex2bmc(idst)).gt.zero) then
                    if(outgoing_offsets(idst).lt.0) cycle

                    !print *,'Setting diffsrc for dst', outgoing_offsets(idst),'= src, edir', isrc, xedir(i1+isrc),&
                    !  '* idst, p2bmc coeff', idst, diff_plex2bmc(idst), dir2diff(diff_plex2bmc(idst))
                    xb(i1+outgoing_offsets(idst)) = xb(i1+outgoing_offsets(idst)) + &
                      xedir(i1+isrc) * dir2diff(diff_plex2bmc(idst))
                  endif
                enddo
              endif
            enddo ! enddo iface

            call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
          enddo

          call set_Edir_srfc_reflection(edirSection, xedir)

          call VecRestoreArrayReadF90(ledirVec, xedir, ierr); call CHKERR(ierr)
          call ISRestoreIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
          call DMRestoreLocalVector(edirdm, ledirVec, ierr); call CHKERR(ierr)
        end subroutine

        subroutine set_Edir_srfc_reflection(edirSection, xedir)
          type(tPetscSection), intent(in) :: edirSection
          real(ireals), pointer, intent(in) :: xedir(:)
          type(tIS) :: srfc_ids
          integer(iintegers), pointer :: xi(:)
          integer(iintegers) :: i, iface, offset_Eup, offset_Edir, offset_srfc
          real(ireals), pointer :: xalbedo(:)

          call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)

          call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', BOTFACE, srfc_ids, ierr); call CHKERR(ierr)
          if (srfc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
          else
            call ISGetIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
            do i = 1, size(xi)
              iface = xi(i)
              ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
              ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
              call PetscSectionGetFieldOffset(ediffSection, iface, i0, offset_Eup, ierr); call CHKERR(ierr)
              call PetscSectionGetOffset(edirSection, iface, offset_Edir, ierr); call CHKERR(ierr)

              call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

              if(ldebug) then
                if(xb(i1+offset_Eup).gt.zero) then
                  print *,'Srfc Reflection of Edir',offset_Eup, offset_Edir, offset_srfc, &
                    '::', xedir(i1+offset_Edir), '*', xalbedo(i1+offset_srfc), '::', xb(i1+offset_Eup)
                  call CHKERR(1_mpiint, 'hah, I am setting '// &
                    'edir -> eup but found a value in src vec that is already larger zero. '// &
                    ' I did not expect that?!')
                endif
              endif

              xb(i1+offset_Eup) = xb(i1+offset_Eup) + xedir(i1+offset_Edir) * xalbedo(i1+offset_srfc)
            enddo
            call ISRestoreIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
          endif
          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)

        end subroutine

        subroutine set_thermal_source()
          type(tVec) :: thermal_src_vec
          real(ireals), pointer :: xsrc(:)

          integer(iintegers) :: geom_offset, face_offset

          integer(iintegers) :: fStart, fEnd, idof, iface, num_dof
          real(ireals) :: area

          if(.not.present(plckVec)) return
          if(.not.present(srfc_emission)) return
          if(.not.allocated(plckVec)) return

          call DMGetLocalVector(ediffdm, thermal_src_vec, ierr); call CHKERR(ierr)
          call VecSet(thermal_src_vec, zero, ierr); call CHKERR(ierr)
          call VecGetArrayF90(thermal_src_vec, xsrc, ierr); call CHKERR(ierr)

          call set_cell_emissions(xsrc)
          call thermal_srfc_emission(xsrc)

          ! Scaling from [W/m2] to Energy [W]
          call DMPlexGetDepthStratum(ediffdm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
          do iface = fStart, fEnd-1
            call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
            area = geoms(i1+geom_offset)

            call PetscSectionGetOffset(ediffSection, iface, face_offset, ierr); call CHKERR(ierr)
            call PetscSectionGetDof(ediffSection, iface, num_dof, ierr); call CHKERR(ierr)
            do idof = 1, num_dof
              xsrc(face_offset+idof) = xsrc(face_offset+idof) * area
            enddo
          enddo

          xb = xb + xsrc

          call VecRestoreArrayF90(thermal_src_vec, xsrc, ierr); call CHKERR(ierr)
          call DMRestoreLocalVector(ediffdm, thermal_src_vec, ierr); call CHKERR(ierr)
        end subroutine
        subroutine set_cell_emissions(xsrc)
          real(ireals), pointer :: xsrc(:)

          integer(iintegers), pointer :: faces_of_cell(:)
          integer(iintegers) :: face_plex2bmc(5)
          integer(iintegers) :: diff_plex2bmc(8)
          integer(iintegers) :: wedge_offset, geom_offset

          integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)
          integer(iintegers) :: i, j, icell, icol, iface, idof, numDof
          real(ireals) :: dz, coeff(8**2) ! coefficients for each src=[1..8] and dst[1..8]

          integer(iintegers), pointer :: xinoutdof(:)
          real(ireals), pointer :: xplanck(:)
          real(ireals) :: diff2diff(8), emissivity, area_bot, area_top

          call ISGetIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
          call VecGetArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)

          do icell = plex%cStart, plex%cEnd-1

            call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
            call get_in_out_dof_offsets(solver%IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)

            call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)

            do i = 1, size(faces_of_cell)
              face_plex2bmc(int(wedgeorient(wedge_offset+i3+i), iintegers)) = i
            enddo
            call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff, ierr)
            if(ierr.eq.OPP_1D_RETCODE) then
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
              area_top = geoms(i1+geom_offset)
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
              area_bot = geoms(i1+geom_offset)
              if(area_top.gt.area_bot) then
                ! send overhanging energy to side faces
                coeff(7*8+[2,4,6]) = coeff(7*8+[2,4,6]) + coeff(7*8+1) * (one - area_bot / area_top) / i3
                ! reduce transmission bc receiver is smaller than top face
                coeff(7*8+1) = coeff(7*8+1) * area_bot / area_top
              else
                ! send overhanging energy to side faces
                coeff([3,5,7]) = coeff([3,5,7]) + coeff(8) * (one - area_top / area_bot) / i3
                ! transmission from bot to top face
                coeff(8) = coeff(8) * area_top / area_bot
              endif
            endif
            call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

            i = 1
            do j = 1, size(faces_of_cell)
              iface = faces_of_cell(j)
              call PetscSectionGetDof(ediffSection, iface, numDof, ierr); call CHKERR(ierr)

              do idof = 1, numDof/2
                icol = outgoing_offsets(i)

                diff2diff = coeff(diff_plex2bmc(i): size(coeff): i8)

                emissivity = max(zero, one - sum(diff2diff))

                xsrc(i1+icol) = xplanck(i1+icell) * pi * emissivity / (numDof/2)
                i = i+1
              enddo
            enddo
            call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
          enddo ! icell
          call VecRestoreArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)
          call ISRestoreIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
        end subroutine

        subroutine thermal_srfc_emission(xsrc)
          real(ireals), pointer :: xsrc(:)
          real(ireals), pointer :: xalbedo(:), xsrfc_emission(:)
          type(tIS) :: srfc_ids
          integer(iintegers), pointer :: xi(:)
          integer(iintegers) :: i, iface, offset_Eup, offset_srfc

          if(.not.allocated(srfc_emission)) &
            call CHKERR(1_mpiint, 'srfc emission has to be allocated to compute thermal radiation')
          call VecGetArrayReadF90(srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)

          call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
          call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', BOTFACE, srfc_ids, ierr); call CHKERR(ierr)
          if (srfc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
          else
            call ISGetIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
            do i = 1, size(xi)
              iface = xi(i)
              ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
              ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
              call PetscSectionGetFieldOffset(ediffSection, iface, i0, offset_Eup, ierr); call CHKERR(ierr)

              call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

              xsrc(i1+offset_Eup) = xsrfc_emission(i1+offset_srfc) * (one-xalbedo(i1+offset_srfc)) * pi
            enddo
            call ISRestoreIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
          endif

          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
          call VecRestoreArrayReadF90(srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
        end subroutine
      end subroutine

      subroutine face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)
        integer(iintegers), intent(in) :: face_plex2bmc(:)  ! mapping from faces_of_cell to bmc_faces ordering (dim=5)
        integer(iintegers), intent(out) :: diff_plex2bmc(:) ! mapping from 8 diff streams in plex ordering to bmc dof (dim=8)
        integer(iintegers) :: i, j

        ! basic mapping would be e.g.:
        ! top face  -> 1
        ! base_face -> 2,3
        ! left_face -> 4,5
        ! right_face-> 6,7
        ! bot face  -> 8

        ! however, we also need to consider the side face permuations

        j = 1
        do i = 1, size(face_plex2bmc)
          select case(face_plex2bmc(i)) ! case switch on bmc indices
          case(1)
            diff_plex2bmc(j) = 1
            j = j+1
          case(5)
            diff_plex2bmc(j) = 8
            j = j+1
          case(2)
            diff_plex2bmc(j:j+1) = [2,3]
            j = j+2
          case(3)
            diff_plex2bmc(j:j+1) = [4,5]
            j = j+2
          case(4)
            diff_plex2bmc(j:j+1) = [6,7]
            j = j+2
          end select
        enddo
      end subroutine

    subroutine solve_plex_rt(dm, b, A, ksp, x, ksp_residual_history, prefix)
      type(tDM), intent(inout) :: dm
      type(tVec), allocatable, intent(in) :: b
      type(tMat), allocatable, intent(in) :: A
      type(tKSP), allocatable, intent(inout) :: ksp
      type(tVec), allocatable, intent(inout) :: x
      real(ireals), allocatable, intent(inout), optional :: ksp_residual_history(:)
      character(len=*),optional :: prefix

      real(ireals),parameter :: rtol=1e-5_ireals, rel_atol=1e-5_ireals
      integer(iintegers),parameter  :: maxiter=1000
      real(ireals) :: atol
      type(tPC) :: prec

      integer(iintegers) :: Nrows_global

      integer(iintegers), parameter :: Nmaxhistory=1000
      integer(mpiint) :: comm, myid, numnodes, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)

      if(.not.allocated(b)) call CHKERR(1_mpiint, 'Src Vector has to be allocated before running solve')
      if(.not.allocated(A)) call CHKERR(1_mpiint, 'System Matrix has to be allocated before running solve')

      if(.not.allocated(ksp)) then
        allocate(ksp)
        call KSPCreate(comm, ksp, ierr); call CHKERR(ierr)
        if(present(prefix)) then
          call KSPAppendOptionsPrefix(ksp, trim(prefix), ierr); call CHKERR(ierr)
        endif

        call MatGetSize(A, Nrows_global, PETSC_NULL_INTEGER, ierr); call CHKERR(ierr)
        call imp_allreduce_min(comm, rel_atol*Nrows_global, atol)
        atol = max(1e-8_ireals, atol)

        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

        if(myid.eq.0.and.ldebug) &
          print *,'Setup KSP -- tolerances:',rtol,atol,'::',rel_atol, Nrows_global

        call KSPSetType(ksp,KSPFGMRES,ierr); call CHKERR(ierr)
        call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr); call CHKERR(ierr)
        call KSPGetPC(ksp,prec,ierr); call CHKERR(ierr)
        if(numnodes.eq.0) then
          call PCSetType(prec, PCILU, ierr); call CHKERR(ierr)
        else
          call PCSetType(prec, PCBJACOBI, ierr); call CHKERR(ierr)
        endif

        call KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT_REAL, maxiter, ierr); call CHKERR(ierr)

        call KSPSetDM(ksp, dm, ierr); call CHKERR(ierr)
        call KSPSetDMActive(ksp, PETSC_FALSE, ierr); call CHKERR(ierr)
      endif

      if(present(ksp_residual_history)) then
        if(.not.allocated(ksp_residual_history)) &
          allocate(ksp_residual_history(Nmaxhistory), source=-one)
        call KSPSetResidualHistory(ksp, ksp_residual_history, Nmaxhistory, PETSC_TRUE, ierr); call CHKERR(ierr)
      endif

      if(.not.allocated(x)) then
        allocate(x)
        call VecDuplicate(b, x, ierr); call CHKERR(ierr)
      endif

      call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)
      call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)
      call KSPSetUp(ksp, ierr); call CHKERR(ierr)

      call KSPSolve(ksp, b, x, ierr); call CHKERR(ierr)

      call handle_diverged_solve()

      call handle_reuse_solver()

      contains
        subroutine handle_reuse_solver()
          logical :: ldestroy_solver, lflg
          type(tMat) :: Amat, Pmat

          ldestroy_solver = .False.
          call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-destroy_solver",&
            ldestroy_solver, lflg, ierr) ;call CHKERR(ierr)
          if(ldestroy_solver) then
            if(myid.eq.0.and.ldebug) print *,'Zeroing KSP Operators...'
            call KSPGetOperators(ksp, Amat, Pmat, ierr); call CHKERR(ierr)
            call MatZeroEntries(Amat, ierr); call CHKERR(ierr)
            call MatZeroEntries(Pmat, ierr); call CHKERR(ierr)
            if(myid.eq.0.and.ldebug) print *,'Zeroing KSP Operators... done'
          endif
        end subroutine
        subroutine handle_diverged_solve()
          KSPConvergedReason :: reason
          KSPType :: old_ksp_type
          integer(iintegers) :: iter
          integer(mpiint) :: myid
          character(len=*), parameter :: alternate_prefix='diverged_alternate_'
          character(len=default_str_len) :: old_prefix

          type(tMat) :: A2

          call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)
          if(reason.le.0) then
            call PetscObjectViewFromOptions(b  , PETSC_NULL_VEC, '-show_diverged_b', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(x  , PETSC_NULL_VEC, '-show_diverged_x', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(A  , PETSC_NULL_MAT, '-show_diverged_A', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(ksp, PETSC_NULL_KSP, '-show_diverged_ksp', ierr); call CHKERR(ierr)
          endif

          if(reason.le.0) then
            call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
            if(myid.eq.0.and.ldebug) &
              print *,myid,'Resetted initial guess to zero and try again with gmres:'
            call VecSet(x,zero,ierr) ;call CHKERR(ierr)
            call KSPGetType(ksp,old_ksp_type,ierr); call CHKERR(ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr) ;call CHKERR(ierr)

            call KSPGetOptionsPrefix(ksp, old_prefix, ierr); call CHKERR(ierr)
            call KSPAppendOptionsPrefix(ksp, trim(alternate_prefix), ierr); call CHKERR(ierr)
            call DMCreateMatrix(dm, A2, ierr); call CHKERR(ierr)
            call MatCopy(A, A2, DIFFERENT_NONZERO_PATTERN, ierr); call CHKERR(ierr)
            call KSPSetOperators(ksp, A2, A2, ierr); call CHKERR(ierr)
            call KSPSetFromOptions(ksp, ierr) ;call CHKERR(ierr)

            call KSPSetUp(ksp,ierr) ;call CHKERR(ierr)
            call KSPSolve(ksp,b,x,ierr) ;call CHKERR(ierr)
            call KSPGetIterationNumber(ksp,iter,ierr) ;call CHKERR(ierr)
            call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)

            ! And return to normal solver...
            call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)
            call KSPSetOptionsPrefix(ksp, trim(old_prefix), ierr); call CHKERR(ierr)
            call KSPSetType(ksp, old_ksp_type,ierr) ;call CHKERR(ierr)
            call KSPSetFromOptions(ksp, ierr) ;call CHKERR(ierr)
            call KSPSetUp(ksp, ierr) ;call CHKERR(ierr)
            if(myid.eq.0.and.ldebug) &
              print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
          endif

          if(reason.le.0) then
            call CHKERR(1_mpiint, '***** SOLVER did NOT converge :( ********'//itoa(reason))
          endif
        end subroutine

    end subroutine

    subroutine scale_facevec(plex, face_dm, globalfaceVec, lW_to_Wm2)
    type(t_plexgrid), intent(in) :: plex
    type(tDM), intent(in) :: face_dm
    type(tVec), intent(inout) :: globalfaceVec
    logical, intent(in) :: lW_to_Wm2    ! convert from W to W/m2 or vice versa

    type(tVec) :: faceVec
    real(ireals), pointer :: xv(:)

    type(tPetscSection) :: geomSection, faceSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    integer(iintegers) :: fStart, fEnd, geom_offset, face_offset, iface, idof, num_dof

    real(ireals) :: area

    integer(mpiint) :: myid, ierr

    !if(ldebug) print *,'plex_rt::scale_facevec...'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(face_dm, faceSection, ierr); call CHKERR(ierr)

    call DMGetLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)

    call DMGlobalToLocalBegin(face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

    call VecGetArrayF90(faceVec, xv, ierr); call CHKERR(ierr)

    call DMPlexGetDepthStratum(face_dm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
    do iface = fStart, fEnd-1
      call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
      area = geoms(i1+geom_offset)

      call PetscSectionGetOffset(faceSection, iface, face_offset, ierr); call CHKERR(ierr)
      call PetscSectionGetDof(faceSection, iface, num_dof, ierr); call CHKERR(ierr)
      if(lW_to_Wm2) then
        do idof = 1, num_dof
          xv(face_offset+idof) = xv(face_offset+idof) / area
        enddo
      else
        do idof = 1, num_dof
          xv(face_offset+idof) = xv(face_offset+idof) * area
        enddo
      endif
    enddo

    call VecRestoreArrayF90(faceVec, xv, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(face_dm, faceVec, INSERT_VALUES, globalfaceVec, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (face_dm, faceVec, INSERT_VALUES, globalfaceVec, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::scale_facevec... finished'
    end subroutine

  subroutine create_edir_mat(solver, plex, OPP, kabs, ksca, g, sundir, A)
    class(t_plex_solver), intent(in) :: solver
    type(t_plexgrid), intent(inout) :: plex
    class(t_optprop), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g
    real(ireals), intent(in) :: sundir(3)
    type(tMat), allocatable, intent(inout) :: A

    type(tPetscSection) :: sec

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: icell, iface, irow, icol, idst, fStart, fEnd
    real(ireals) :: c

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: geom_offset, wedge_offset

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)

    real(ireals) :: zenith, azimuth, area_top, area_bot

    real(ireals) :: dir2dir(5)
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)
    integer(iintegers) :: numSrc, numDst

    real(ireals) :: dz, coeff(5**2) ! coefficients for each src=[1..5] and dst[1..5]

    logical, parameter :: lonline=.False.
    logical :: lflg, ldestroy_mat

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    !if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...'

    if(.not.allocated(plex%geom_dm)) call CHKERR(1_mpiint, 'geom_dm has to allocated in order to create an Edir Matrix')
    if(.not.allocated(plex%wedge_orientation_dm)) call CHKERR(1_mpiint, 'wedge_orientation_dm has to allocated in order to create an Edir Matrix')
    if(.not.allocated(plex%edir_dm)) call CHKERR(1_mpiint, 'edir_dm has to allocated in order to create an Edir Matrix')
    if(.not.allocated(kabs  )) call CHKERR(1_mpiint, 'kabs   has to be allocated')
    if(.not.allocated(ksca  )) call CHKERR(1_mpiint, 'ksca   has to be allocated')
    if(.not.allocated(g     )) call CHKERR(1_mpiint, 'g      has to be allocated')

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call DMGetSection(plex%edir_dm, sec, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(sec, PETSC_NULL_SECTION, '-show_edir_loc_section', ierr); call CHKERR(ierr)

    if(allocated(A)) then
      ldestroy_mat = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-destroy_mat",&
        ldestroy_mat, lflg, ierr) ;call CHKERR(ierr)
      if(ldestroy_mat) then
        if(myid.eq.0.and.ldebug) print *,'Destroying Old Direct Matrix'
        call MatDestroy(A, ierr); call CHKERR(ierr)
        deallocate(A)
      endif
    endif

    if(.not.allocated(A)) then
      allocate(A)
      if(myid.eq.0.and.ldebug) print *,'Creating Direct Matrix...'
      call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)
      call MatSetBlockSize(A,i1,ierr); call CHKERR(ierr)
      if(myid.eq.0.and.ldebug) print *,'Creating Direct Matrix... done'
    endif

    call MatZeroEntries(A, ierr); call CHKERR(ierr)

    call set_side_incoming_boundary_condition()

    ! Set Diagonal Entries
    call DMPlexGetDepthStratum(plex%edir_dm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
    do iface = fStart, fEnd-1
      call PetscSectionGetDof(sec, iface, numDst, ierr); call CHKERR(ierr)
      do idst = 0, numDst-1
        call PetscSectionGetOffset(sec, iface+idst, irow, ierr); call CHKERR(ierr)
        call MatSetValuesLocal(A, i1, irow, i1, irow, [one], INSERT_VALUES, ierr); call CHKERR(ierr)
      enddo
    enddo

    do icell = plex%cStart, plex%cEnd-1
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)
      zenith  = wedgeorient(wedge_offset+i1)
      azimuth = wedgeorient(wedge_offset+i2)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
        lsrc(iface) = nint(wedgeorient(wedge_offset+i8+iface), iintegers) .eq. i1
      enddo

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      call PetscLogEventBegin(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
        dz, wedgeorient(wedge_offset+14:wedge_offset+19), .True., coeff, ierr, &
        angles=[real(rad2deg(azimuth), irealLUT), real(rad2deg(zenith), irealLUT)])

      if(ldebug) then
        do iface=1,i5
          dir2dir = coeff(face_plex2bmc(iface):size(coeff):i5)
          if(sum(dir2dir).gt.one+10*sqrt(epsilon(one))) then
            print *,iface,': bmcface', face_plex2bmc(iface), 'dir2dir gt one', dir2dir,':',sum(dir2dir)
            call CHKERR(1_mpiint, 'energy conservation violated! '//ftoa(sum(dir2dir)))
          endif
        enddo
      endif

      if(ierr.eq.OPP_1D_RETCODE) then
        ! for the case of 1D spherical radiative transfer,
        ! need to consider the change in area between upper and lower face
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
        area_top = geoms(i1+geom_offset)
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
        area_bot = geoms(i1+geom_offset)

        ! super compensate direct radiation
        coeff(21) = coeff(21) * (area_top / area_bot)
      endif
      call PetscLogEventEnd(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)

      do iface = 1, size(faces_of_cell)

        if(lsrc(iface)) then
          ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering
          dir2dir = coeff(face_plex2bmc(iface):size(coeff):i5)

          call PetscSectionGetDof(sec, faces_of_cell(iface), numSrc, ierr); call CHKERR(ierr)
          if(numSrc.eq.i0) cycle
          if(ldebug.and.numSrc.gt.i1) call CHKERR(1_mpiint, 'direct dof more than 1 is not implemented')
          call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr)

          do idst = 1, size(faces_of_cell)
            if(.not.lsrc(idst)) then
              c = -dir2dir(face_plex2bmc(idst))

              call PetscSectionGetDof(sec, faces_of_cell(idst), numDst, ierr); call CHKERR(ierr)
              if(numDst.eq.i0) cycle

              call PetscSectionGetOffset(sec, faces_of_cell(idst), irow, ierr); call CHKERR(ierr)
              !print *,'isrc', face_plex2bmc(iface), 'idst', face_plex2bmc(idst), &
              !  'if', iface, 'id', idst, &
              !  'srcface->dstface', faces_of_cell(iface),faces_of_cell(idst), &
              !  'col -> row', icol, irow, c
              call MatSetValuesLocal(A, i1, irow, i1, icol, c, INSERT_VALUES, ierr); call CHKERR(ierr)
              if(ldebug.and.irow.eq.icol) call CHKERR(1_mpiint, &
                'src and dst are the same :( ... should not happen here row '//itoa(irow)//' col '//itoa(icol))
            endif
          enddo
        endif

      enddo ! enddo iface

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Mdir', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    contains
      subroutine set_side_incoming_boundary_condition()
        use m_helper_functions, only: rotation_matrix_local_basis_to_world, rotation_matrix_world_to_local_basis
        use m_plex_grid, only: compute_local_wedge_ordering, compute_local_vertex_coordinates
        type(tIS) :: IS_side_faces
        integer(iintegers), pointer :: iside_faces(:), cell_support(:)
        real(ireals) :: dz
        integer(iintegers) :: o, i, j, icell, iface_top, iface_side, irow, icol, idst, isrc
        real(ireals) :: dir2dir(5)
        real(ireals) :: inv_sundir(3), mean_dx
        real(ireals) :: top_face_normal(3), side_face_normal(3), U(3), Mrot(3,3), coords_2d(6)

        integer(iintegers) :: upper_face, base_face, left_face, right_face, bottom_face, forient(5)

        call DMGetStratumIS(plex%edir_dm, 'DomainBoundary', SIDEFACE, IS_side_faces, ierr); call CHKERR(ierr)
        if (IS_side_faces.eq.PETSC_NULL_IS) return ! dont have side boundary faces
        call ISGetIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)

        do o = 1, size(iside_faces)
          iface_side = iside_faces(o)

          call PetscSectionGetDof(sec, iface_side, numDst, ierr); call CHKERR(ierr)
          if(numDst.eq.i0) cycle

          call DMPlexGetSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)
          icell = cell_support(1)
          call DMPlexRestoreSupport(plex%edir_dm, iface_side, cell_support, ierr); call CHKERR(ierr)

          call get_inward_face_normal(iface_side, icell, geomSection, geoms, side_face_normal)
          if(is_solar_src(side_face_normal, sundir)) then

            call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
            iface_top = faces_of_cell(1)

            call get_inward_face_normal(iface_top, icell, geomSection, geoms, top_face_normal)
            U = cross_3d(top_face_normal, side_face_normal)

            Mrot = rotation_matrix_world_to_local_basis(top_face_normal, side_face_normal, U)
            inv_sundir = matmul(Mrot, sundir)
            inv_sundir = [inv_sundir(1), -inv_sundir(2), -inv_sundir(3)] ! rotate by 180
            Mrot = rotation_matrix_local_basis_to_world(top_face_normal, side_face_normal, U)
            inv_sundir = matmul(Mrot, inv_sundir)

            !print *,sundir,'inv_sundir',inv_sundir

            !print *,'faces_of_cell', faces_of_cell, 'sideface', iface_side

            call compute_local_wedge_ordering(plex, icell, &
              geomSection, geoms, inv_sundir, &
              lsrc, zenith, azimuth, mean_dx, &
              upper_face, bottom_face, &
              base_face, left_face, right_face)

            forient = [upper_face, base_face, left_face, right_face, bottom_face]
            do i=1, size(faces_of_cell)
              face_plex2bmc(forient(i)) = i
            enddo

            call compute_local_vertex_coordinates(plex, &
              faces_of_cell(upper_face), &
              faces_of_cell(base_face), &
              faces_of_cell(left_face), &
              faces_of_cell(right_face), &
              coords_2d)

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            !print *,'az', rad2deg(azimuth), rad2deg(zenith), dz
            !print *,'upper_face, bottom_face', upper_face, bottom_face
            !print *,'base_face, left_face, right_face', base_face, left_face, right_face
            !print *,'coords2d', coords_2d
            !print *,'face_plex2bmc', face_plex2bmc

            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              dz, coords_2d, .True., coeff, ierr, &
              angles=[real(rad2deg(azimuth),irealLUT), real(rad2deg(zenith),irealLUT)])

            !do i=1,5
            !  print *,'coeff for bmc dst',i, coeff((i-1)*5+1:i*5), ':', sum(coeff((i-1)*5+1:i*5))
            !enddo

            do i = 1, 5
              if(faces_of_cell(i).eq.iface_side) then
                call PetscSectionGetOffset(sec, faces_of_cell(i), irow, ierr); call CHKERR(ierr)
                idst = face_plex2bmc(i)
                !print *,i, 'idst', idst
                dir2dir = coeff((idst-1)*5+1:idst*5)
                !print *,'dir2dir', dir2dir
                do j = 1, 5
                  isrc = face_plex2bmc(j)

                  call PetscSectionGetOffset(sec, faces_of_cell(j), icol, ierr); call CHKERR(ierr)

                  !call get_inward_face_normal(faces_of_cell(j), icell, geomSection, geoms, side_face_normal)
                  !print *,i,j,':', idst, isrc, ': faces', faces_of_cell(j),'->', faces_of_cell(i), &
                  !  'icol', icol, '->', irow, '(',dir2dir(isrc),')', &
                  !  is_solar_src(side_face_normal, sundir), is_solar_src(side_face_normal, inv_sundir)

                  !if(.not.is_solar_src(side_face_normal, inv_sundir)) cycle

                  !if(dir2dir(isrc).le.zero) then
                  !  call CHKERR(1_mpiint, 'encountered a small value... weird')
                  !endif
                  if(dir2dir(isrc).gt.zero) then
                    call MatSetValuesLocal(A, i1, irow, i1, icol, -dir2dir(isrc), INSERT_VALUES, ierr); call CHKERR(ierr)
                    if(ldebug.and.irow.eq.icol) call CHKERR(1_mpiint, &
                      'src and dst are the same :( ... should not happen here row '//itoa(irow)//' col '//itoa(icol))
                  endif
                enddo

              endif
            enddo ! i
          endif ! boundary side face is sunlit
        enddo
      end subroutine
  end subroutine

  subroutine create_ediff_mat(solver, plex, OPP, kabs, ksca, g, albedo, A)
    class(t_plex_solver), intent(in) :: solver
    type(t_plexgrid), intent(in) :: plex
    class(t_optprop), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g ! cell1_dm
    type(tVec), allocatable, intent(in) :: albedo        ! srfc_boundary_dm
    type(tMat), allocatable, intent(inout) :: A

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: j, icell, iface, irow, icol, i_inoff
    real(ireals) :: c

    type(tPetscSection) :: ediffSection, geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: wedge_offset, geom_offset

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: diff_plex2bmc(8)

    real(ireals) :: diff2diff(8)

    real(ireals) :: dz, coeff(8**2) ! coefficients for each src=[1..8] and dst[1..8]
    integer(iintegers), pointer :: xinoutdof(:)
    real(ireals) :: area_top, area_bot
    real(ireals), parameter :: coeff_norm_err_tolerance=one+100*sqrt(epsilon(one))

    logical :: lflg, ldestroy_mat, lreset_mat

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    !if(ldebug.and.myid.eq.0) print *,'plex_rt::create_ediff_mat...'

    if(.not.allocated(plex%ediff_dm)) call CHKERR(1_mpiint, 'ediff_dm has to allocated in order to create an Ediff Matrix')
    if(.not.allocated(kabs  )) call CHKERR(1_mpiint, 'kabs   has to be allocated')
    if(.not.allocated(ksca  )) call CHKERR(1_mpiint, 'ksca   has to be allocated')
    if(.not.allocated(g     )) call CHKERR(1_mpiint, 'g      has to be allocated')
    if(.not.allocated(albedo)) call CHKERR(1_mpiint, 'albedo has to be allocated')

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call DMGetSection(plex%ediff_dm, ediffSection, ierr); call CHKERR(ierr)

    if(allocated(A)) then
      ldestroy_mat = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-destroy_mat",&
        ldestroy_mat, lflg, ierr) ;call CHKERR(ierr)
      if(ldestroy_mat) then
        if(myid.eq.0.and.ldebug) print *,'Destroying Old Diffuse Matrix'
        call MatDestroy(A, ierr); call CHKERR(ierr)
        deallocate(A)
      endif
    endif

    if(.not.allocated(A)) then
      allocate(A)
      if(myid.eq.0.and.ldebug) print *,'Creating Diffuse Matrix...'
      call DMCreateMatrix(plex%ediff_dm, A, ierr); call CHKERR(ierr)
      call MatSetBlockSize(A,i2,ierr); call CHKERR(ierr)
      call MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, ierr); call CHKERR(ierr)
      call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr); call CHKERR(ierr)
      if(myid.eq.0.and.ldebug) print *,'Creating Diffuse Matrix... done'
    endif

    lreset_mat = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-reset_mat",&
      lreset_mat, lflg, ierr) ;call CHKERR(ierr)
    if(lreset_mat) then
      if(myid.eq.0.and.ldebug) print *,'Resetting Diffuse Matrix'
      call MatResetPreallocation(A, ierr); call CHKERR(ierr)
    endif
    call MatZeroEntries(A, ierr); call CHKERR(ierr)

    call set_boundary_conditions()
    call set_diagonal_entries()

    call ISGetIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
    do icell = plex%cStart, plex%cEnd-1
      call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      call get_in_out_dof_offsets(solver%IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)
      !print *,'icell', icell, 'faces_of_cell', faces_of_cell
      !print *,'icell', icell, 'incoming_offsets', incoming_offsets
      !print *,'icell', icell, 'outgoing_offsets', outgoing_offsets

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
      enddo
      call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      !print *,'icell',icell,': foc',faces_of_cell
      !print *,'icell',icell,':',face_plex2bmc
      call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
        dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff, ierr)

      if(ldebug_optprop) then
        do i_inoff = 1, size(incoming_offsets)
          diff2diff = coeff(diff_plex2bmc(i_inoff):size(coeff):i8)
          if(sum(diff2diff).gt.coeff_norm_err_tolerance) then
            print *,i_inoff,': bmcface', diff_plex2bmc(i_inoff), 'diff2diff gt one', diff2diff, &
              ':', sum(diff2diff), 'l1d', ierr.eq.OPP_1D_RETCODE, 'tol', coeff_norm_err_tolerance
            call CHKERR(1_mpiint, '1 energy conservation violated! '//ftoa(sum(diff2diff)))
          endif
        enddo
      endif

      if(ierr.eq.OPP_1D_RETCODE) then
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
        area_top = geoms(i1+geom_offset)
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
        area_bot = geoms(i1+geom_offset)

        !if(area_top.gt.area_bot) then
        !  ! send overhanging energy to side faces
        !  coeff(7*8+[2,4,6]) = coeff(7*8+[2,4,6]) + coeff(7*8+1) * (one - area_bot / area_top) / i3
        !  ! reduce transmission bc receiver is smaller than top face
        !  coeff(7*8+1) = coeff(7*8+1) * area_bot / area_top
        !else
        !  ! send overhanging energy to side faces
        !  coeff([3,5,7]) = coeff([3,5,7]) + coeff(8) * (one - area_top / area_bot) / i3
        !  ! transmission from bot to top face
        !  coeff(8) = coeff(8) * area_top / area_bot
        !endif
      endif
      call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

      do i_inoff = 1, size(incoming_offsets)
        icol = incoming_offsets(i_inoff)
        if(icol.lt.0) cycle

        ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering
        diff2diff = coeff(diff_plex2bmc(i_inoff):size(coeff):i8)
        if(ldebug_optprop) then
          if(sum(diff2diff).gt.coeff_norm_err_tolerance) then
            print *,i_inoff,': bmcface', diff_plex2bmc(i_inoff), 'diff2diff gt one', diff2diff,':',sum(diff2diff)
            call CHKERR(1_mpiint, '2 energy conservation violated! '//ftoa(sum(diff2diff)))
          endif
        endif
        if(sum(diff2diff).gt.one) &
          diff2diff = diff2diff/(sum(diff2diff)+10*epsilon(diff2diff))

        do j = 1, size(outgoing_offsets)
          c = -diff2diff(diff_plex2bmc(j))
          irow = outgoing_offsets(j)
          if(c.ge.-epsilon(zero) .or. irow.lt.0) cycle

          !print *,'icell',icell,'isrc,jdst',i_inoff,j,'icol', icol, 'irow', irow, '=>', c
          call MatSetValuesLocal(A, i1, irow, i1, icol, c, INSERT_VALUES, ierr); call CHKERR(ierr)
          if(ldebug.and.irow.eq.icol) call CHKERR(1_mpiint, &
            'src and dst are the same :( ... should not happen here row '//itoa(irow)//' col '//itoa(icol))
        enddo

      enddo ! enddo i_inoff

      call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo
    call ISRestoreIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Mediff', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    !if(ldebug.and.myid.eq.0) print *,'plex_rt::create_ediff_mat...finished'
    contains
      subroutine set_boundary_conditions()
        type(tIS) :: bc_ids
        integer(iintegers), pointer :: xi(:)
        integer(iintegers) :: i, iface, istream, offset_Ein, offset_Eout, offset_srfc
        integer(iintegers) :: num_fields, numDof
        type(tPetscSection) :: srfcSection

        real(ireals), pointer :: xalbedo(:)
        real(ireals) :: sideward_bc_coeff
        logical :: lflg

        call PetscSectionGetNumFields(ediffSection, num_fields, ierr); call CHKERR(ierr)

        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', BOTFACE, bc_ids, ierr); call CHKERR(ierr)
        if (bc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
        else
          call DMGetSection(plex%srfc_boundary_dm, srfcSection, ierr); CHKERRQ(ierr)
          call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
          call ISGetIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          do i = 1, size(xi)
            iface = xi(i)
            call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

            do istream = 1, num_fields
              call PetscSectionGetFieldDof(ediffSection, iface, istream-i1, numDof, ierr); call CHKERR(ierr)
              if(numDof.eq.i2) then
                ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
                ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
                call PetscSectionGetFieldOffset(ediffSection, iface, istream-i1, offset_Ein, ierr); call CHKERR(ierr)
                offset_Eout = offset_Ein+i1

                call MatSetValuesLocal(A, i1, offset_Ein, i1, offset_Eout, -xalbedo(i1+offset_srfc), INSERT_VALUES, ierr); call CHKERR(ierr)
              endif
            enddo
          enddo
          call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
        endif


        sideward_bc_coeff = one
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-sideward_bc_coeff", sideward_bc_coeff, lflg, ierr)  ; call CHKERR(ierr)

        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', SIDEFACE, bc_ids, ierr); call CHKERR(ierr)
        if (bc_ids.eq.PETSC_NULL_IS.or.sideward_bc_coeff.le.zero) then ! dont have surface points
        else
          call ISGetIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          do i = 1, size(xi)
            iface = xi(i)
            call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)
            do istream = 1, num_fields
              call PetscSectionGetFieldDof(ediffSection, iface, istream-i1, numDof, ierr); call CHKERR(ierr)
              if(numDof.eq.i2) then
                ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
                ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
                call PetscSectionGetFieldOffset(ediffSection, iface, istream-i1, offset_Ein, ierr); call CHKERR(ierr)
                offset_Eout = offset_Ein+i1

                call MatSetValuesLocal(A, i1, offset_Ein, i1, offset_Eout, -sideward_bc_coeff, INSERT_VALUES, ierr); call CHKERR(ierr)
                if(ldebug.and.offset_Ein.eq.offset_Eout) call CHKERR(1_mpiint, &
                  'src and dst are the same :( ... should not happen here'// &
                  ' row '//itoa(offset_Ein)// &
                  ' col '//itoa(offset_Eout))
              endif
            enddo
          enddo
          call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
        endif
      end subroutine
      subroutine set_diagonal_entries()
        integer(iintegers) :: irow, is, ie
        call PetscSectionGetOffsetRange(ediffSection, is, ie, ierr); call CHKERR(ierr)
        do irow = is, ie-1
          call MatSetValuesLocal(A, i1, irow, i1, irow, one, INSERT_VALUES, ierr); call CHKERR(ierr)
        enddo
      end subroutine
  end subroutine
  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(OPP, kabs, ksca, g, dz, wedge_coords, ldir, coeff, ierr, angles)
    class(t_optprop), intent(in) :: OPP
    real(ireals), intent(in)     :: kabs, ksca, g
    real(ireals),intent(in)      :: dz
    real(ireals), intent(in)     :: wedge_coords(:) ! coordinates of upper triangle pts A,B,C in in (x,y)
    logical,intent(in)           :: ldir
    real(ireals),intent(out)     :: coeff(:)
    integer(mpiint), intent(out) :: ierr

    real(irealLUT),intent(in),optional  :: angles(2)

    real(ireals) :: dkabs, dksca, dg
    real(irealLUT) :: aspect, tauz, w0, dx, relcoords(6), CC(size(coeff))
    integer, parameter :: iC1=4, iC2=5
    real(ireals), parameter :: dither=.0_ireals ! add dither on relcoords
    real(ireals) :: Rdither


    dx = real(wedge_coords(3), irealLUT)

    dkabs = kabs
    dksca = ksca
    dg    = g
    call delta_scale( dkabs, dksca, dg, opt_f=dg)

    aspect = real(dz / dx, irealLUT)
    tauz = real((dkabs+dksca) * dz, irealLUT)
    if(approx(tauz,0._irealLUT)) then
      w0 = 0
    else
      w0 = real(dksca / (dkabs+dksca), irealLUT)
    endif

    relcoords = real(wedge_coords / dx, irealLUT)

    if(present(angles)) then
      tauz = max(OPP%OPP_LUT%dirconfig%dims(1)%vrange(1), &
        min(OPP%OPP_LUT%dirconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%OPP_LUT%dirconfig%dims(2)%vrange(1), &
        min(OPP%OPP_LUT%dirconfig%dims(2)%vrange(2), w0))
    else
      tauz = max(OPP%OPP_LUT%diffconfig%dims(1)%vrange(1), &
        min(OPP%OPP_LUT%diffconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%OPP_LUT%diffconfig%dims(2)%vrange(1), &
        min(OPP%OPP_LUT%diffconfig%dims(2)%vrange(2), w0))
    endif

    if(dither.gt.zero) then
      call random_number(Rdither)
      relcoords(5) = relcoords(5) * real((one+(dither*(-one+2*Rdither))), irealLUT)
      relcoords(6) = relcoords(6) * real((one+(dither*(-one+2*Rdither))), irealLUT)
    endif

    relcoords(5) = max(OPP%OPP_LUT%diffconfig%dims(iC1)%vrange(1), &
                   min(OPP%OPP_LUT%diffconfig%dims(iC1)%vrange(2), relcoords(5)))
    relcoords(6) = max(OPP%OPP_LUT%diffconfig%dims(iC2)%vrange(1), &
                   min(OPP%OPP_LUT%diffconfig%dims(iC2)%vrange(2), relcoords(6)))

    !print *,'DEBUG Coeffs', tauz, w0, g, aspect, angles, ':', norm(wedge_coords(3:4)-wedge_coords(1:2)), ':', relcoords
    call OPP%get_coeff(tauz, w0, real(dg, irealLUT), aspect, ldir, CC, ierr, angles=angles, wedge_coords=relcoords)
    !print *,'DEBUG Lookup Coeffs for', tauz, w0, g, aspect, angles, ':', norm(wedge_coords(3:4)-wedge_coords(1:2)), ':', relcoords, '::', coeff
    coeff = CC
    if(ldebug) then
      if(any(coeff.lt.zero).or.any(coeff.gt.one)) then
        print *,'Lookup Coeffs for', aspect, tauz, w0, dg, angles,'::', coeff
        call CHKERR(1_mpiint, 'Found corrupted coefficients!')
      endif
    endif
  end subroutine

  subroutine restore_solution(solver, solution, time)
    class(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution
    real(ireals),intent(in),optional :: time

    integer(mpiint) :: ierr
    type(tVec) :: abso_old

    if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Create working vec to determine difference between old and new absorption vec
      call DMGetGlobalVector(solver%plex%abso_dm, abso_old, ierr)   ; call CHKERR(ierr)
      call VecCopy( solution%abso, abso_old, ierr)     ; call CHKERR(ierr)
    endif

    call compute_absorption(solver, solution)
    solution%lchanged = .False.

    call update_absorption_norms_for_adaptive_spectral_integration()

    contains
      subroutine update_absorption_norms_for_adaptive_spectral_integration()
        real(ireals) :: norm1, norm2, norm3
        integer(mpiint) :: myid, comm, ierr
        if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
          call VecAXPY(abso_old, -one, solution%abso, ierr); call CHKERR(ierr) ! overwrite abso_old with difference to new one
          call VecNorm(abso_old, NORM_1, norm1, ierr)       ; call CHKERR(ierr)
          call VecNorm(abso_old, NORM_2, norm2, ierr)       ; call CHKERR(ierr)
          call VecNorm(abso_old, NORM_INFINITY, norm3, ierr); call CHKERR(ierr)

          call DMRestoreGlobalVector(solver%plex%abso_dm, abso_old, ierr)   ; call CHKERR(ierr)

          ! Save norm for later analysis
          solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
          solution%twonorm = eoshift ( solution%twonorm, shift = -1) !shift all values by 1 to the right
          solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

          solution%maxnorm( 1 ) = norm3
          solution%twonorm( 1 ) = norm2
          solution%time( 1 )    = time

          if(ldebug) then
            call PetscObjectGetComm(solver%plex%abso_dm, comm, ierr); call CHKERR(ierr)
            call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
            if(myid.eq.0) &
              print *,'Updating error statistics for solution',solution%uid,'at time', &
                      solution%time(1),':: norm 1,2,inf',norm1,norm2,norm3,'[W] :: ', &
                      'hr_norm approx:',norm3*86.1,'[K/d]'
          endif
        endif !present(time) .and. solver%lenable_solutions_err_estimates
      end subroutine
  end subroutine

  subroutine compute_absorption(solver, solution)
    class(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution

    integer(mpiint) :: ierr

    if(.not.solution%lset) call CHKERR(1_mpiint, 'compute_absorption needs to get an initialized solution obj')

    if(.not.solution%lchanged) return
    if(solution%lsolar_rad .and. .not.allocated(solution%edir)) &
      call CHKERR(1_mpiint, 'solution%lsolar_rad true but edir not allocated')
    if(.not.allocated(solution%ediff)) call CHKERR(1_mpiint, 'ediff vec not allocated')

    ! Make sure we have the radiation vecs in plain energy units
    call scale_flx(solver%plex, &
      solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
      solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
      solution, lWm2=.False., logevent=solver%logs%scale_flx)

    if(.not.allocated(solution%abso)) then
      allocate(solution%abso)
      call DMCreateGlobalVector(solver%plex%abso_dm, solution%abso, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad) then
      call PetscLogEventBegin(solver%logs%compute_absorption, ierr)
      call compute_edir_absorption(solver%plex, solution%edir, solution%abso)
      call PetscLogEventEnd(solver%logs%compute_absorption, ierr)

      call PetscLogEventBegin(solver%logs%debug_output, ierr)
      call PetscObjectSetName(solution%abso, 'abso_direct_'//itoa(solution%uid), ierr);call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_abso_direct', ierr); call CHKERR(ierr)
      call PetscLogEventEnd(solver%logs%debug_output, ierr)
    endif

    call PetscLogEventBegin(solver%logs%compute_absorption, ierr)
    call compute_ediff_absorption(solver%plex, solver%IS_diff_in_out_dof, solution%ediff, solution%abso)
    call PetscLogEventEnd(solver%logs%compute_absorption, ierr)

    call PetscLogEventBegin(solver%logs%debug_output, ierr)
    call PetscObjectSetName(solution%abso, 'abso_'//itoa(solution%uid), ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
    call PetscLogEventEnd(solver%logs%debug_output, ierr)
  end subroutine

  subroutine compute_edir_absorption(plex, edir, abso)
    type(t_plexgrid), allocatable, intent(inout) :: plex
    type(tVec), intent(in) :: edir
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_edir

    real(ireals), pointer :: xedir(:), xabso(:)

    type(tPetscSection) :: abso_section, edir_section

    integer(iintegers) :: cStart, cEnd, idof, numDof
    integer(iintegers) :: i, icell, iface
    integer(iintegers),pointer :: faces_of_cell(:)
    integer(mpiint) :: myid, ierr

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:), wedgeorient(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset, edir_offset, wedge_offset

    real(ireals) :: inv_volume

    logical :: lsrc ! is src or destination of solar beam

    if(.not.allocated(plex)) stop 'called compute_edir_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%edir_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_edir_absorption with an unallocated abso vec')

    !if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      do i = 1, size(faces_of_cell)
        iface = faces_of_cell(i)
        lsrc = int(wedgeorient(wedge_offset+i8+i), iintegers) .eq. i1

        call PetscSectionGetDof(edir_section, iface, numDof, ierr); call CHKERR(ierr)
        do idof = 0, numDof-1
          call PetscSectionGetOffset(edir_section, iface, edir_offset, ierr); call CHKERR(ierr)

          if(lsrc) then ! sun is shining into this face
            xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1) * inv_volume
            !print *,icell,'Dir Abso',abso_offset,'=',xabso(abso_offset+i1),'( did +',xedir(edir_offset+i1)*inv_volume,')'
          else
            xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1) * inv_volume
            !print *,icell,'Dir Abso',abso_offset,'=',xabso(abso_offset+i1),'( did -',xedir(edir_offset+i1)*inv_volume,')'
          endif
        enddo
      enddo
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end subroutine

  subroutine compute_ediff_absorption(plex, IS_diff_in_out_dof, ediff, abso)
    type(t_plexgrid), allocatable, intent(in) :: plex
    type(tIS), intent(in) :: IS_diff_in_out_dof
    type(tVec), intent(in) :: ediff
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_ediff

    real(ireals), pointer :: xediff(:), xabso(:)

    type(tPetscSection) :: abso_section, ediff_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface, i, j, idof, numDof
    integer(mpiint) :: myid, ierr

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset

    integer(iintegers), pointer :: faces_of_cell(:), xinoutdof(:)

    real(ireals) :: inv_volume

    if(.not.allocated(plex)) stop 'called compute_ediff_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%ediff_dm)) call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_ediff_absorption with an unallocated abso vec')

    !if(ldebug) print *,'plex_rt::compute_ediff_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call ISGetIndicesF90(IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
    do icell = cStart, cEnd-1

      call get_in_out_dof_offsets(IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)

      call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      i = 1
      do j = 1, size(faces_of_cell)
        iface = faces_of_cell(j)
        call PetscSectionGetDof(ediff_section, iface, numDof, ierr); call CHKERR(ierr)
        do idof = 1, numDof/2
          !print *,'Current Abso', xabso(i1+abso_offset)
          xabso(i1+abso_offset) = xabso(i1+abso_offset) + (xediff(i1+incoming_offsets(i)) * inv_volume)
          xabso(i1+abso_offset) = xabso(i1+abso_offset) - (xediff(i1+outgoing_offsets(i)) * inv_volume)
          !print *,'iface', iface, numDof, incoming_offsets(i), outgoing_offsets(i), 'abso', &
          !  '+', xediff(i1+incoming_offsets(i)) * inv_volume, &
          !  '-', xediff(i1+outgoing_offsets(i)) * inv_volume, &
          !  '=', xabso(i1+abso_offset)
          i = i+1
        enddo
      enddo
      call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo
    call ISRestoreIndicesF90(IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::compute_ediff_absorption....finished'
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution)
    class(t_plex_solver)    :: solver
    type(t_state_container) :: solution

    real(ireals),allocatable :: dtau(:), Blev(:), Edn(:),Eup(:)

    type(tIS) :: boundary_ids
    integer(iintegers), pointer :: xitoa(:), cell_support(:)
    integer(iintegers), allocatable :: cell_idx(:)
    integer(iintegers) :: i, k, icell, iface, voff, ke1, geom_offset
    real(ireals) :: dz
    real(ireals), pointer :: xkabs(:), xalbedo(:), xplck(:), xsrfc_emission(:), xediff(:), xgeoms(:)
    type(tPetscSection) :: ediff_section, geom_section

    integer(mpiint) :: ierr

    if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
    if( .not. allocated(solver%plck) ) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

    associate( plex => solver%plex )

    call DMGetStratumIS(plex%ediff_dm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
    if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
    else
      allocate(dtau(plex%Nlay))
      allocate(Blev(plex%Nlay+1))
      allocate(Edn (plex%Nlay+1))
      allocate(Eup (plex%Nlay+1))

      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geom_section, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)

      call ISGetIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)
      do i = 1, size(xitoa)
        iface = xitoa(i)
        call DMPlexGetSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        icell = cell_support(1)
        call DMPlexRestoreSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

        call get_consecutive_vertical_cell_idx(plex, icell, cell_idx)
        ke1 = size(cell_idx)+1
        do k=1,ke1-1
          icell = cell_idx(k)

          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          dtau(k) = xkabs(i1+icell) * dz
          Blev(k) = xplck(i1+icell)
        enddo
        Blev(k) = xsrfc_emission(i)

        call schwarzschild(dtau, xalbedo(i), Edn, Eup, Blev)

        do k = 0, ke1-2
          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          xediff(i1+voff) = Edn(i1+k)
          xediff(i2+voff) = Eup(i1+k)
        enddo
        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
        call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
        xediff(i1+voff) = Eup(i1+k)
        xediff(i2+voff) = Edn(i1+k)
      enddo
      call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)


      !Schwarzschild solver returns fluxes as [W/m^2]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged  = .True.
    endif ! TOA boundary ids

    end associate
  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(plex, kabs, ksca, g, albedo, sundir, solution, plck, srfc_emission)
    type(t_plexgrid), intent(in) :: plex
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), intent(in), optional :: plck, srfc_emission
    real(ireals), intent(in) :: sundir(:)
    type(t_state_container) :: solution

    real(ireals),allocatable :: vdtau(:), vw0(:), vg(:), Blev(:), Edir(:), Edn(:),Eup(:)

    type(tIS) :: boundary_ids
    integer(iintegers), pointer :: xitoa(:), cell_support(:)
    integer(iintegers), allocatable :: cell_idx(:)
    integer(iintegers) :: i, k, icell, iface, voff, ke1, geom_offset
    real(ireals) :: dz, theta0
    real(ireals), pointer :: xksca(:), xkabs(:), xg(:), xalbedo(:), xplck(:), xsrfc_emission(:)
    real(ireals), pointer :: xedir(:), xediff(:), xgeoms(:)
    type(tPetscSection) :: edir_section, ediff_section, geom_section

    real(ireals) :: face_normal(3)
    logical :: lthermal, lsolar

    integer(mpiint) :: ierr

    if(solution%lsolar_rad) then
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad .and. norm(sundir).le.zero) then
      return
    endif

    lsolar = solution%lsolar_rad
    lthermal = .not.solution%lsolar_rad

    call DMGetStratumIS(plex%edir_dm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
    if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
    else
      allocate(vdtau(plex%Nlay))
      allocate(vw0  (plex%Nlay))
      allocate(vg   (plex%Nlay))
      allocate(Edir(plex%Nlay+1))
      allocate(Edn (plex%Nlay+1))
      allocate(Eup (plex%Nlay+1))

      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geom_section, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)


      if(lthermal) then
        if(.not.present(plck)) call CHKERR(1_mpiint, 'have to provide planck vec to compute thermal rad with twostr')
        if(.not.present(srfc_emission)) call CHKERR(1_mpiint, 'have to provide srfc_emission to compute thermal rad with twostr')
        allocate(Blev(plex%Nlay+1))
        call VecGetArrayReadF90(plck, xplck, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
      endif

      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
      if(lsolar) then
        call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
        call VecGetArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
      endif

      call ISGetIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)
      do i = 1, size(xitoa)
        iface = xitoa(i)

        call DMPlexGetSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        icell = cell_support(1)
        call DMPlexRestoreSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        call get_inward_face_normal(iface, icell, geom_section, xgeoms, face_normal)
        theta0 = angle_between_two_vec(face_normal, sundir)

        call get_consecutive_vertical_cell_idx(plex, icell, cell_idx)
        ke1 = size(cell_idx)+1
        do k=1,ke1-1
          icell = cell_idx(k)

          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          vdtau(k) = (xkabs(i1+icell) + xksca(i1+icell)) * dz
          vw0(k)   = xksca(i1+icell) / max(epsilon(vw0), xkabs(i1+icell) + xksca(i1+icell))
          vg(k)    = xg(i1+icell)

          call delta_scale_optprop(vdtau(k), vw0(k), vg(k), vg(k))

          if(lthermal) Blev(k) = xplck(i1+icell)
        enddo
        if(lthermal) Blev(k) = xsrfc_emission(i)

        if(solution%lsolar_rad) then
          call delta_eddington_twostream(vdtau,vw0,vg,cos(theta0),norm(sundir)*cos(theta0),xalbedo(i), &
            Edir, Edn, Eup )
        else
          call delta_eddington_twostream(vdtau,vw0,vg,cos(theta0),norm(sundir)*cos(theta0),xalbedo(i), &
            Edir, Edn, Eup, Blev)
        endif

        if(lsolar) then
          do k = 0, ke1-1
            call PetscSectionGetOffset(edir_section, iface+k, voff, ierr); call CHKERR(ierr)
            xedir(i1+voff) = edir(i1+k)
          enddo
        endif

        do k = 0, ke1-2
          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          xediff(i1+voff) = Edn(i1+k)
          xediff(i1+voff+i1) = Eup(i1+k)
        enddo
        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
        call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
        xediff(i1+voff) = Eup(i1+k)
        xediff(i1+voff+i1) = Edn(i1+k)
      enddo
      call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)

      if(lsolar) then
        call VecRestoreArrayF90(solution%edir, xedir, ierr); call CHKERR(ierr)
      endif
      if(lthermal) then
        call VecRestoreArrayReadF90(plck, xplck, ierr); call CHKERR(ierr)
        call VecRestoreArrayReadF90(srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
      endif

      !Twostream solver returns fluxes as [W/m^2]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged  = .True.
    endif ! TOA boundary ids
  end subroutine

  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(plex, &
      dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, &
      diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, &
      solution, lWm2, logevent)
    type(t_plexgrid), intent(in) :: plex
    type(tVec), allocatable, intent(inout) :: dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2
    type(tVec), allocatable, intent(inout) :: diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2
    type(t_state_container),intent(inout)  :: solution  !< @param solution container with computed fluxes
    logical,intent(in)                     :: lWm2      !< @param determines direction of scaling, if true, scale to W/m**2
    PetscLogEvent, intent(in), optional    :: logevent
    integer(mpiint) :: ierr

    call gen_scale_vecs()

    if(present(logevent)) then
      call PetscLogEventBegin(logevent, ierr); call CHKERR(ierr)
    endif

    if(solution%lsolar_rad) then
      if(solution%lWm2_dir .neqv. lWm2) then
        if(lWm2) then
          call VecPointwiseMult(solution%edir, solution%edir, dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        else
          call VecPointwiseMult(solution%edir, solution%edir, dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        endif
        solution%lWm2_dir = lWm2
      endif
    endif

    if(solution%lWm2_diff .neqv. lWm2) then
      if(lWm2) then
        call VecPointwiseMult(solution%ediff, solution%ediff, diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      else
        call VecPointwiseMult(solution%ediff, solution%ediff, diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      endif
      solution%lWm2_diff = lWm2
    endif
    if(present(logevent)) then
      call PetscLogEventEnd(logevent, ierr); call CHKERR(ierr)
    endif

    contains
      subroutine gen_scale_vecs()
        if(solution%lsolar_rad) then
          if(.not.allocated(dir_scalevec_Wm2_to_W)) then
            allocate(dir_scalevec_Wm2_to_W)
            call VecDuplicate(solution%edir, dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
            call VecSet(dir_scalevec_Wm2_to_W, one, ierr); call CHKERR(ierr)
            call scale_facevec(plex, plex%edir_dm, dir_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
          endif
          if(.not.allocated(dir_scalevec_W_to_Wm2)) then
            allocate(dir_scalevec_W_to_Wm2)
            call VecDuplicate(dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
            call VecSet(dir_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
            call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
              dir_scalevec_W_to_Wm2, &
              dir_scalevec_W_to_Wm2, &
              dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
          endif
        endif

        if(.not.allocated(diff_scalevec_Wm2_to_W)) then
          allocate(diff_scalevec_Wm2_to_W)
          call VecDuplicate(solution%ediff, diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
          call VecSet(diff_scalevec_Wm2_to_W, one, ierr); call CHKERR(ierr)
          call scale_facevec(plex, plex%ediff_dm, diff_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
        endif
        if(.not.allocated(diff_scalevec_W_to_Wm2)) then
          allocate(diff_scalevec_W_to_Wm2)
          call VecDuplicate(diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
          call VecSet(diff_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
          call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
            diff_scalevec_W_to_Wm2, &
            diff_scalevec_W_to_Wm2, &
            diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        endif
      end subroutine
  end subroutine

  subroutine plexrt_get_result(solver, redn, reup, rabso, redir, opt_solution_uid)
    class(t_plex_solver), intent(inout) :: solver
    real(ireals), allocatable, dimension(:,:), intent(inout) :: redn,reup,rabso
    real(ireals), allocatable, dimension(:,:), intent(inout), optional :: redir
    integer(iintegers),optional,intent(in) :: opt_solution_uid

    type(tIS) :: toa_ids
    integer(iintegers) :: Ncol, ke1, uid, i, k, iface, icell, voff
    integer(iintegers), pointer :: xtoa_faces(:), cell_support(:)

    type(tPetscSection) :: abso_section, ediff_section, edir_section
    real(ireals), pointer :: xediff(:), xedir(:), xabso(:)

    integer(mpiint) :: myid, ierr

    call PetscLogEventBegin(solver%logs%get_result, ierr)

    uid = get_arg(0_iintegers, opt_solution_uid)
    if(.not.solver%solutions(uid)%lset) &
      call CHKERR(1_mpiint, 'You tried to retrieve results from a solution uid which has not yet been calculated')

    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

    ke1 = solver%plex%Nlay+1

    call check_size_or_allocate(redn , [ke1, Ncol])
    call check_size_or_allocate(reup , [ke1, Ncol])
    call check_size_or_allocate(rabso, [ke1-1, Ncol])
    if(present(redir)) call check_size_or_allocate(redir, [ke1, Ncol])

    associate(solution => solver%solutions(uid))

    if(.not.allocated(solution%ediff)) call CHKERR(1_mpiint, 'ediff vec not allocated')
    if(.not.allocated(solution%abso)) call CHKERR(1_mpiint, 'abso vec not allocated')
    if(solution%lchanged) call CHKERR(1_mpiint, 'tried to get results from unrestored solution -- call restore_solution first')

    if(present(redir) .and. .not.solution%lsolar_rad) &
      call CHKERR(1_mpiint, 'you asked for direct radiation but solution does not have it')

    call scale_flx(solver%plex, &
      solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
      solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
      solution, lWm2=.True., logevent=solver%logs%scale_flx)

    if(present(redir)) then
      call DMGetSection(solver%plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    endif
    call DMGetSection(solver%plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
    call DMGetSection(solver%plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    if(present(redir)) then
      call VecGetArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
    endif
    call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
    call VecGetArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

    call ISGetIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
    do i = 1, size(xtoa_faces)
      iface = xtoa_faces(i)
      do k = 0, ke1-2
        call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
        redn(i1+k, i) = xediff(i1+voff)
        reup(i1+k, i) = xediff(i2+voff)
      enddo
      ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
      call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
      redn(i1+k, i) = xediff(i2+voff)
      reup(i1+k, i) = xediff(i1+voff)

      ! Fill Absorption Vec
      call DMPlexGetSupport(solver%plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr)
      icell = cell_support(1)
      call DMPlexRestoreSupport(solver%plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr)

      do k = 0, ke1-2
        call PetscSectionGetOffset(abso_section, icell+k, voff, ierr); call CHKERR(ierr)
        rabso(i1+k, i) = xabso(i1+voff)
      enddo

      if(present(redir)) then
        do k = 0, ke1-1
          call PetscSectionGetFieldOffset(edir_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          redir(i1+k, i) = xedir(i1+voff)
        enddo

        if(ldebug) then
          if(reup(ke1,i)-sqrt(epsilon(reup)) .gt. abs(redir(ke1,i))+abs(redn(ke1,i))) then
            call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)
            do k = 1, size(redir,dim=1)
              print *,'sgr ::', myid, 'i', i, 'k', k, redir(k,i), redn(k,i), reup(k,i)
            enddo
            call CHKERR(1_mpiint, 'solar get_result Eup cannot be bigger than Edir+Edn!')
          endif
        endif
      endif
    enddo

    if(present(redir)) then
      call VecRestoreArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
    endif
    call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

    if(ldebug) then
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)
      if(myid.eq.0) then
        if(present(redir)) then
          print *,'Get Result, k     Edir                   Edn                      Eup                  abso'
          do k = 1, ke1-1
            print *,k, redir(k,1), redn(k,1), reup(k,1), rabso(k,1)!, &
              !redir(k,1)-redir(k+1,1)+redn(k,1)-redn(k+1,1)-reup(k,1)+reup(k+1,1)
          enddo
          print *,k, redir(k,1), redn(k,1), reup(k,1)
        else
          print *,'Get Result, k     Edn                    Eup                      abso'
          do k = 1, ke1-1
            print *,k, redn(k,1), reup(k,1), rabso(k,1)!, &
              !redn(k,1)-redn(k+1,1)-reup(k,1)+reup(k+1,1)
          enddo
          print *,k, redn(k,1), reup(k,1)
        endif
      endif
    endif

    end associate
    call PetscLogEventEnd(solver%logs%get_result, ierr)

    contains
      subroutine check_size_or_allocate(arr, expected_shape)
        real(ireals), allocatable :: arr(:,:)
        integer(iintegers), intent(in) :: expected_shape(:)

        if(allocated(arr)) then
          if(.not.all(shape(arr).eq.expected_shape)) then
            print *,'Expected results array to be of size', expected_shape, 'but you provided mem with', shape(redn)
            call CHKERR(1_mpiint, 'Wrong shape for results array Edn')
          endif
        else
          if(size(expected_shape).eq.2) then
            allocate(arr(expected_shape(1),expected_shape(2)))
          else
            call CHKERR(1_mpiint, 'here is a place where the code is not yet implemented.. but its easy.. go ahead')
          endif
        endif
      end subroutine
  end subroutine
end module
