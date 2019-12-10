module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options, ltwostr_only, lschwarzschild, lcalc_nca

  use m_helper_functions, only: CHKERR, CHKWARN, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, strF2C, get_arg, meanval, &
    vec_proj_on_plane, cross_3d, rotation_matrix_world_to_local_basis, &
    approx, swap, delta_scale, delta_scale_optprop, itoa, ftoa, cstr, &
    imp_reduce_mean, imp_min_mean_max, rotation_matrix_around_axis_vec

  use m_data_parameters, only : ireals, iintegers, mpiint, irealLUT, imp_ireals, &
    i0, i1, i2, i3, i4, i5, i6, i7, i8, default_str_len, &
    zero, one, pi, EXP_MINVAL, EXP_MAXVAL, nan

  use m_plex_grid, only: t_plexgrid, compute_face_geometry, &
    setup_cell1_dmplex, setup_edir_dmplex, setup_ediff_dmplex, setup_abso_dmplex, &
    orient_face_normals_along_sundir, compute_wedge_orientation, is_solar_src, get_inward_face_normal, &
    facevec2cellvec, icell_icon_2_plex, iface_top_icon_2_plex, get_consecutive_vertical_cell_idx, &
    get_top_bot_face_of_cell, destroy_plexgrid, determine_diff_incoming_outgoing_offsets, &
    TOAFACE, BOTFACE, SIDEFACE

  use m_optprop, only : t_optprop_wedge, OPP_1D_RETCODE, OPP_TINYASPECT_RETCODE, &
    t_optprop_wedge_5_8, &
    t_optprop_rectilinear_wedge_5_8, &
    t_optprop_wedge_18_8
  use m_optprop_LUT, only : find_lut_dim_by_name
  use m_optprop_parameters, only : ldebug_optprop

  use m_pprts_base, only : t_state_container, prepare_solution, destroy_solution, &
    t_solver_log_events, setup_log_events

  use m_plex_rt_base, only: &
    t_plex_solver, &
    t_plex_solver_5_8, &
    t_plex_solver_rectilinear_5_8, &
    t_plex_solver_18_8, &
    t_dof

  use m_plexrt_external_solvers, only: plexrt_schwarz, plexrt_twostream, plexrt_disort, plexrt_NCA_wrapper
  use m_plex2rayli, only: rayli_wrapper
  use m_schwarzschild, only: B_eff
  use m_netcdfio, only: ncwrite
  use m_petsc_helpers, only: hegedus_trick

  implicit none

  private

  public :: init_plex_rt_solver, run_plex_rt_solver, &
    compute_face_geometry, set_plex_rt_optprop, destroy_plexrt_solver, &
    plexrt_get_result, scale_flx

  !logical, parameter :: ldebug=.False.
  logical, parameter :: ldebug=.True.

  contains
    subroutine init_plex_rt_solver(plex, solver)
      type(t_plexgrid), intent(in) :: plex
      class(t_plex_solver), allocatable, intent(inout) :: solver
      integer(mpiint) :: myid, ierr
      logical :: lplexrt_skip_loadLUT, lflg

      call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
      if(ldebug.and.myid.eq.0) print *,'Init_plex_rt_solver ... '
      call read_commandline_options(plex%comm)

      if(.not.allocated(solver)) &
        call CHKERR(1_mpiint, 'Should not call init_plex_rt_solver with an unallocated solver object')

      select type(solver)
      class is (t_plex_solver_5_8)
        allocate(t_optprop_wedge_5_8::solver%OPP)
        solver%dirtop%dof = 1
        solver%dirtop%area_divider = 1
        solver%dirside%dof = 1
        solver%dirside%area_divider = 1

        solver%difftop%dof = 2
        solver%difftop%area_divider = 1
        solver%diffside%dof = 4
        solver%diffside%area_divider = 1

      class is (t_plex_solver_rectilinear_5_8)
        allocate(t_optprop_rectilinear_wedge_5_8::solver%OPP)
        solver%dirtop%dof = 1
        solver%dirtop%area_divider = 1
        solver%dirside%dof = 1
        solver%dirside%area_divider = 1

        solver%difftop%dof = 2
        solver%difftop%area_divider = 1
        solver%diffside%dof = 4
        solver%diffside%area_divider = 1

      class is (t_plex_solver_18_8)
        allocate(t_optprop_wedge_18_8::solver%OPP)
        solver%dirtop%dof = 3
        solver%dirtop%area_divider = 3
        solver%dirside%dof = 4
        solver%dirside%area_divider = 4

        solver%difftop%dof = 2
        solver%difftop%area_divider = 1
        solver%diffside%dof = 4
        solver%diffside%area_divider = 1

      class default
        call CHKERR(1_mpiint, 'unexpected type for solver')
      end select
      solver%dirdof = solver%dirtop%dof*2 + solver%dirside%dof*3
      solver%diffdof = solver%difftop%dof*2 + solver%diffside%dof*3

      allocate(solver%plex)
      solver%plex = plex

      call setup_log_events(solver%logs, 'plexrt')

      lplexrt_skip_loadLUT=.False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_skip_loadLUT",&
        lplexrt_skip_loadLUT, lflg, ierr) ;call CHKERR(ierr)
      if(.not.lplexrt_skip_loadLUT) call solver%OPP%init(plex%comm)

      call setup_abso_dmplex (solver%plex%dm, solver%plex%abso_dm)
      call setup_edir_dmplex (solver%plex, solver%plex%dm, &
        i1, i0, i1, &
        solver%plex%horizface1_dm)
      call setup_edir_dmplex (solver%plex, solver%plex%dm, &
        solver%dirtop%dof, solver%dirside%dof, i1, &
        solver%plex%edir_dm)
      call setup_ediff_dmplex(solver%plex, solver%plex%dm, &
        solver%difftop%dof/2, solver%diffside%dof/2, i2, &
        solver%plex%ediff_dm)

      call plexrt_view_geometry(plex%comm)
      if(ldebug.and.myid.eq.0) print *,'Init_plex_rt_solver ... done'
      contains
        subroutine plexrt_view_geometry(comm)
          integer(mpiint) :: comm
          logical :: lview, lflg
          type(tPetscSection) :: geom_section
          real(ireals), pointer :: xgeoms(:) ! pointer to coordinates vec
          integer(iintegers) :: geom_offset
          integer(iintegers), pointer :: xitoa(:), cell_support(:)
          type(tIS) :: boundary_ids

          real(ireals), allocatable, dimension(:,:) :: top_area, bot_area, dz, vol
          real(ireals), dimension(3) :: mtop_area, mbot_area, mdz, mvol
          integer(iintegers), allocatable :: cell_idx(:)
          integer(iintegers) :: i, k, icell, iface

          integer(mpiint) :: myid, ierr

          if(.not.allocated(solver%plex%geom_dm)) &
            call CHKERR(1_mpiint, 'run_plex_rt_solver::geom_dm has to be allocated first')

          lview=.False.
          call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_view_geometry",&
            lview, lflg, ierr) ;call CHKERR(ierr)
          if(.not.lview) return

          associate( plex => solver%plex, geom_dm => solver%plex%geom_dm )

            call DMGetStratumIS(plex%edir_dm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
            if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
              allocate(dz(plex%Nlay,0), vol(plex%Nlay,0), top_area(plex%Nlay,0), bot_area(plex%Nlay,0))
            else
              call DMGetSection(geom_dm, geom_section, ierr); CHKERRQ(ierr)
              call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)
              call ISGetIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

              allocate(&
                dz (plex%Nlay,size(xitoa)), &
                vol(plex%Nlay,size(xitoa)), &
                top_area(plex%Nlay,size(xitoa)), &
                bot_area(plex%Nlay,size(xitoa)))

              do i = 1, size(xitoa)
                iface = xitoa(i)

                call DMPlexGetSupport(geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
                icell = cell_support(1)
                call DMPlexRestoreSupport(geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
                call get_consecutive_vertical_cell_idx(plex, icell, cell_idx)
                do k = 0, size(cell_idx)-1
                  icell = cell_idx(i1+k)

                  call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
                  dz(i1+k,i) = xgeoms(i1+geom_offset)

                  call PetscSectionGetFieldOffset(geom_section, icell, i2, geom_offset, ierr); call CHKERR(ierr)
                  vol(i1+k,i) = xgeoms(i1+geom_offset)

                  call PetscSectionGetFieldOffset(geom_section, iface+k, i2, geom_offset, ierr); call CHKERR(ierr)
                  top_area(i1+k,i) = xgeoms(i1+geom_offset)
                  call PetscSectionGetFieldOffset(geom_section, iface+k+1, i2, geom_offset, ierr); call CHKERR(ierr)
                  bot_area(i1+k,i) = xgeoms(i1+geom_offset)
                enddo
              enddo

              call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)
              call VecRestoreArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)
            endif
          end associate

          call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

          if(myid.eq.0) print *,'*        k  '// &
            '                       '//cstr(' dz (min/mean/max)', 'blue')//'                      '// &
            '                       '//cstr(' volume           ', 'red' )//'                      '// &
            '                       '//cstr(' cell_top_area    ', 'blue')//'                      '// &
            '                       '//cstr(' cell_bot_area    ', 'red')

          do k = 1, size(dz,dim=1)
            call imp_min_mean_max(comm, dz      (k,:), mdz      )
            call imp_min_mean_max(comm, vol     (k,:), mvol     )
            call imp_min_mean_max(comm, top_area(k,:), mtop_area)
            call imp_min_mean_max(comm, bot_area(k,:), mbot_area)

            if(myid.eq.0) then
              print *,k, cstr(ftoa(mdz), 'blue'), cstr(ftoa(mvol), 'red'), &
                cstr(ftoa(mtop_area), 'blue'), cstr(ftoa(mbot_area), 'red')
            endif
          enddo
        end subroutine
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
        if(allocated(solver%ksp_solar_dir)) then
          call KSPDestroy(solver%ksp_solar_dir, ierr); call CHKERR(ierr)
          deallocate(solver%ksp_solar_dir)
        endif
        if(allocated(solver%ksp_solar_dir)) then
          call KSPDestroy(solver%ksp_solar_dir, ierr); call CHKERR(ierr)
          deallocate(solver%ksp_solar_dir)
        endif
        if(allocated(solver%ksp_solar_diff)) then
          call KSPDestroy(solver%ksp_solar_diff, ierr); call CHKERR(ierr)
          deallocate(solver%ksp_solar_diff)
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

        if(present(vert_integrated_kabs)) then
          kabs_tot = vert_integrated_kabs / real(solver%plex%Nlay, ireals) / dz
        else
          kabs_tot = rayleigh * (one - w0) * dz
        endif
        if(present(vert_integrated_ksca)) then
          ksca_tot = vert_integrated_ksca / real(solver%plex%Nlay, ireals) / dz
        else
          ksca_tot = rayleigh * w0 *dz
        endif

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
      integer(mpiint) :: ierr

      if(present(opt_xv)) then
        call with_is_vec(opt_xv)
      else
        call ISGetIndicesF90(IS_diff_in_out_dof, xv, ierr); call CHKERR(ierr)
        call with_is_vec(xv)
        call ISRestoreIndicesF90(IS_diff_in_out_dof, xv, ierr); call CHKERR(ierr)
      endif
    contains
      subroutine with_is_vec(xis)
        integer(iintegers), intent(in), pointer :: xis(:)
        integer(iintegers) :: bs, idx_offset

        bs = xis(size(xis))
        if(.not.allocated(incoming_offsets)) allocate(incoming_offsets(bs))
        if(.not.allocated(outgoing_offsets)) allocate(outgoing_offsets(bs))

        idx_offset = icell*bs
        incoming_offsets = xis(i1+idx_offset:idx_offset+bs/2)
        outgoing_offsets = xis(i1+idx_offset+bs/2:idx_offset+bs)
      end subroutine
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
      logical :: luse_rayli, lvacuum_domain_boundary, luse_disort, lflg

      if(.not.allocated(solver)) call CHKERR(1_mpiint, 'run_plex_rt_solver::solver has to be allocated')

      if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'run_plex_rt_solver::plex has to be allocated first')
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(solver%plex%geom_dm)) call CHKERR(1_mpiint, 'geom_dm has to be allocated first')
      if(.not.allocated(solver%plex%srfc_boundary_dm)) call CHKERR(1_mpiint, 'srfc_boundary_dm has to be allocated first')
      if(.not.allocated(solver%kabs  )) call CHKERR(1_mpiint, 'kabs   has to be allocated first')
      if(.not.allocated(solver%ksca  )) call CHKERR(1_mpiint, 'ksca   has to be allocated first')
      if(.not.allocated(solver%g     )) call CHKERR(1_mpiint, 'g      has to be allocated first')
      if(.not.allocated(solver%albedo)) call CHKERR(1_mpiint, 'albedo has to be allocated first')
      if(lthermal) then
        if(.not.allocated(solver%plck)) call CHKERR(1_mpiint, 'planck radiation vec has to be allocated first')
      endif

      if(.not.allocated(solver%plex%geom_dm))  call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      if(.not.allocated(solver%plex%edir_dm))  &
        call CHKERR(1_mpiint, 'solver%plex%edir_dm not allocated, should have happened in init_solver?')
      if(.not.allocated(solver%plex%ediff_dm)) &
        call CHKERR(1_mpiint, 'solver%plex%ediff_dm not allocated, should have happened in init_solver?')
      if(.not.allocated(solver%plex%abso_dm))  &
        call CHKERR(1_mpiint, 'solver%plex%abso_dm not allocated, should have happened in init_solver?')

      if(.not.allocated(solver%IS_diff_in_out_dof)) &
        call setup_IS_diff_in_out_dof(solver%plex, solver%plex%ediff_dm, solver%IS_diff_in_out_dof)

      ! Prepare the space for the solution
      suid = get_arg(i0, opt_solution_uid)

      if(.not.solver%solutions(suid)%lset) then
        call prepare_solution(solver%plex%edir_dm, solver%plex%ediff_dm, solver%plex%abso_dm, &
          lsolar=lsolar, solution=solver%solutions(suid), uid=suid)
      endif

      lvacuum_domain_boundary=.False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        "-plexrt_vacuum_domain_boundary", lvacuum_domain_boundary, lflg, ierr) ;call CHKERR(ierr)
      if(lvacuum_domain_boundary) then
        call vacuum_domain_boundary(solver%plex%cell1_dm, &
          solver%kabs, solver%ksca, solver%g, &
          fill_kabs=zero, fill_ksca=zero, fill_g=zero)
      endif

      call dump_optical_properties(solver%kabs, solver%ksca, solver%g, solver%albedo, &
        plck=solver%plck, postfix='_'//itoa(suid))

      call print_optical_properties_summary(solver%plex%comm, solver%plex%Nlay, &
        solver%kabs, solver%ksca, solver%g, solver%albedo, solver%plck)

      associate( solution => solver%solutions(suid) )

        if(ldebug) print *,'sundir/norm2(sundir)',sundir/norm2(sundir), 'vs', last_sundir, &
          ':', all(approx(last_sundir, sundir/norm2(sundir), sqrt(epsilon(sundir))))
        ! Wedge Orientation is used in solar and thermal case alike
        if(.not.all(approx(last_sundir, sundir/norm2(sundir), sqrt(epsilon(sundir)))).or.&  ! update wedge orientations if sundir has changed
          .not.allocated(solver%plex%wedge_orientation_dm)) then ! or if we lost the info somehow... e.g. happens after destroy_solver
          if(ldebug.and.myid.eq.0) print *,'Updating sun info for sundir:', sundir
          call PetscLogEventBegin(solver%logs%orient_face_normals, ierr)
          call orient_face_normals_along_sundir(solver%plex, sundir)
          call PetscLogEventEnd(solver%logs%orient_face_normals, ierr)

          call PetscLogEventBegin(solver%logs%compute_orientation, ierr)
          call compute_wedge_orientation(solver%plex, sundir, solver%plex%wedge_orientation_dm, &
            solver%plex%wedge_orientation)
          call PetscLogEventEnd(solver%logs%compute_orientation, ierr)

          last_sundir = sundir/norm2(sundir)
        endif

        ! --------- Calculate 1D Radiative Transfer ------------
        if(ltwostr_only) then
          if( (lsolar.eqv..False.) .and. lschwarzschild ) then
            call PetscLogEventBegin(solver%logs%solve_schwarzschild, ierr)
            call plexrt_schwarz(solver, solution)
            call PetscLogEventEnd(solver%logs%solve_schwarzschild, ierr)
          else
            call PetscLogEventBegin(solver%logs%solve_twostream, ierr)
            call plexrt_twostream(solver, solver%plex, &
              lthermal, lsolar, &
              solver%kabs, solver%ksca, solver%g, &
              solver%albedo, sundir, solution, plck=solver%plck)
            call PetscLogEventEnd(solver%logs%solve_twostream, ierr)
          endif

          if(ldebug) print *,'1D calculation done', suid, ':', solution%lsolar_rad, lschwarzschild
          goto 99
        endif

        ! DISORT interface
        luse_disort=.False.
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_use_disort",&
          luse_disort, lflg, ierr) ;call CHKERR(ierr)
        if(luse_disort) then
          call PetscLogEventBegin(solver%logs%solve_twostream, ierr)
          call plexrt_disort(solver, solver%plex, solver%kabs, solver%ksca, solver%g, &
            solver%albedo, sundir, solution, plck=solver%plck)
          call PetscLogEventEnd(solver%logs%solve_twostream, ierr)

          if(ldebug) print *,'1D disort calculation done', suid
          goto 99
        endif


        ! RayLi Raytracer interface
        luse_rayli=.False.
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_use_rayli",&
          luse_rayli, lflg, ierr) ;call CHKERR(ierr)
        if(luse_rayli) then
          call PetscLogEventBegin(solver%logs%solve_rayli, ierr)
          call rayli_wrapper(solver%plex, solver%kabs, solver%ksca, solver%g, &
              solver%albedo, sundir, solution, plck=solver%plck)
          call PetscLogEventEnd(solver%logs%solve_rayli, ierr)
          goto 99
        endif

        if(solution%lsolar_rad) then
          call PetscLogEventBegin(solver%logs%compute_Edir, ierr)
          ! Setup direct source term
          call PetscLogEventBegin(solver%logs%setup_dir_src, ierr)
          call create_edir_src_vec(solver, solver%plex, solver%plex%edir_dm, norm2(sundir), &
            solver%kabs, solver%ksca, solver%g, &
            sundir/norm2(sundir), solver%dirsrc)
          call PetscLogEventEnd(solver%logs%setup_dir_src, ierr)

          ! Output of srcVec
          if(ldebug) &
            call debug_dump_vec(solver%plex%edir_dm, solver%dirsrc, solver%dir_scalevec_W_to_Wm2)

          ! Create Direct Matrix
          call PetscLogEventBegin(solver%logs%setup_Mdir, ierr)
          call create_edir_mat(solver, solver%plex, solver%OPP, solver%kabs, solver%ksca, solver%g, &
            sundir/norm2(sundir), solver%Mdir)
          call PetscLogEventEnd(solver%logs%setup_Mdir, ierr)


          ! Solve Direct Matrix
          call PetscLogEventBegin(solver%logs%solve_Mdir, ierr)
          call solve_plex_rt(solver%plex%edir_dm, solver%dirsrc, solver%Mdir, solver%ksp_solar_dir, solution%edir, &
            ksp_residual_history=solution%dir_ksp_residual_history, prefix='solar_dir_')
          call PetscLogEventEnd(solver%logs%solve_Mdir, ierr)
          call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, &
            '-show_edir_vec_global', ierr); call CHKERR(ierr)
          solution%lWm2_dir = .False.
          solution%lchanged = .True.

          ! Output of Edir Vec
          if(ldebug) &
            call debug_dump_vec(solver%plex%edir_dm, solution%edir, solver%dir_scalevec_W_to_Wm2)

          call PetscLogEventEnd(solver%logs%compute_Edir, ierr)
        endif

        ! Create Diffuse Src
        call PetscLogEventBegin(solver%logs%compute_Ediff, ierr)
        call PetscLogEventBegin(solver%logs%setup_diff_src, ierr)
        call create_ediff_src_vec(solver, solver%plex, solver%OPP, solver%plex%ediff_dm, &
          solver%kabs, solver%ksca, solver%g, solver%albedo, &
          lthermal, lsolar, solver%diffsrc, &
          solver%plex%horizface1_dm, solver%plck, &
          solver%plex%edir_dm, solution%edir)
        call PetscLogEventEnd(solver%logs%setup_diff_src, ierr)

        ! Output of Diffuse Src Vec
        if(ldebug) &
          call debug_dump_vec(solver%plex%ediff_dm, solver%diffsrc, solver%diff_scalevec_W_to_Wm2)

        ! Create Diffuse Matrix
        call PetscLogEventBegin(solver%logs%setup_Mdiff, ierr)
        call create_ediff_mat(solver, solver%plex, solver%OPP, &
          solver%kabs, solver%ksca, solver%g, solver%albedo, solver%Mdiff)
        call PetscObjectViewFromOptions(solver%Mdiff, PETSC_NULL_MAT, &
          '-show_Mediff_'//itoa(solution%uid), ierr); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%setup_Mdiff, ierr)


        ! Solve Diffuse Matrix
        call PetscLogEventBegin(solver%logs%solve_Mdiff, ierr)
        if(solution%lsolar_rad) then
          call solve_plex_rt(solver%plex%ediff_dm, solver%diffsrc, solver%Mdiff, solver%ksp_solar_diff, solution%ediff, &
            ksp_residual_history=solution%diff_ksp_residual_history, prefix='solar_diff_')
        else
          call solve_plex_rt(solver%plex%ediff_dm, solver%diffsrc, solver%Mdiff, solver%ksp_thermal_diff, solution%ediff, &
            ksp_residual_history=solution%diff_ksp_residual_history, prefix='thermal_diff_')
        endif
        call PetscLogEventEnd(solver%logs%solve_Mdiff, ierr)
        call PetscObjectSetName(solution%ediff, 'ediff', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, &
          '-show_ediff_vec_global', ierr); call CHKERR(ierr)
        solution%lWm2_diff = .False.
        solution%lchanged = .True.

        call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)

        99 continue ! this is the quick exit final call where we clean up before the end of the routine

        ! Bring solution into a coherent state, i.e. update absorption etc.
        call restore_solution(solver, solution, opt_solution_time)

      end associate

      contains

        subroutine debug_dump_vec(dm, vec, scalevec)
          type(tDM), intent(in) :: dm
          type(tVec), intent(in) :: vec
          type(tVec), intent(in), optional :: scalevec
          type(tVec) :: cpy
          character(len=default_str_len) :: vecname
          call PetscLogEventBegin(solver%logs%debug_output, ierr)
          call VecDuplicate(vec, cpy, ierr); call CHKERR(ierr)
          call PetscObjectGetName(vec, vecname, ierr); call CHKERR(ierr)
          call PetscObjectSetName(cpy, trim(vecname), ierr); call CHKERR(ierr)

          if(present(scalevec)) then
            call VecPointwiseMult(cpy, vec, scalevec, ierr); call CHKERR(ierr) ! cpy = v1*v2
          else
            call VecCopy(vec, cpy, ierr); call CHKERR(ierr)
          endif

          call facevec2cellvec(dm, cpy)
          call VecDestroy(cpy, ierr); call CHKERR(ierr)
          call PetscLogEventEnd(solver%logs%debug_output, ierr)
        end subroutine
    end subroutine

    subroutine dump_optical_properties(kabs, ksca, g, albedo, plck, postfix)
      character(len=*), intent(in) :: postfix
      type(tVec), allocatable, intent(in) :: kabs, ksca, g, albedo, plck
      call dump_var(kabs, 'kabs')
      call dump_var(ksca, 'ksca')
      call dump_var(g, 'g')
      call dump_var(plck, 'plck')
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

    subroutine print_optical_properties_summary(comm, Nlay, kabs, ksca, g, albedo, plck)
      integer(mpiint), intent(in) :: comm
      integer(iintegers), intent(in) :: Nlay
      type(tVec), allocatable, intent(in) :: kabs, ksca, g, albedo, plck
      real(ireals), pointer, dimension(:) :: xkabs, xksca, xg, xalbedo, xplck
      real(ireals), pointer, dimension(:,:) :: xxkabs, xxksca, xxg, xxalbedo, xxplck
      ! min mean max vals:
      real(ireals), dimension(3) :: mkabs, mksca, mg, malbedo, mplck
      integer(iintegers) :: k
      integer(mpiint) :: myid, ierr
      logical :: lview, lflg

      lview=.False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_view_optprop",&
            lview, lflg, ierr) ;call CHKERR(ierr)
      if(.not.lview) return

      call get_col_vec(albedo, i1    , xalbedo, xxalbedo)
      call get_col_vec(kabs  , Nlay  , xkabs  , xxkabs  )
      call get_col_vec(ksca  , Nlay  , xksca  , xxksca  )
      call get_col_vec(g     , Nlay  , xg     , xxg     )
      call get_col_vec(plck  , Nlay+1, xplck  , xxplck  )


      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      mkabs=nan
      mksca=nan
      mg=nan
      malbedo=nan
      mplck=nan

      if(allocated(plck)) then
        if(myid.eq.0) print *,'*        k  '// &
          '                                '//cstr('kabs(min/mean/max)', 'blue')//'                            '// &
          '                                '//cstr('ksca              ', 'red' )//'                            '// &
          '                                '//cstr('g                 ', 'blue')//'                            '// &
          '                                '//cstr('plck              ', 'red')
      else
        if(myid.eq.0) print *,'*        k  '// &
          '                                '//cstr('kabs(min/mean/max)', 'blue')//'                            '// &
          '                                '//cstr('ksca              ', 'red' )//'                            '// &
          '                                '//cstr('g                 ', 'blue')
    endif
      if(myid.eq.0) print *,cstr('---------------------------------------------------------------------------------------', 'blue')
      do k=1,Nlay
      if(allocated(kabs)) call imp_min_mean_max(comm, xxkabs(k,:), mkabs)
      if(allocated(ksca)) call imp_min_mean_max(comm, xxksca(k,:), mksca)
      if(allocated(g   )) call imp_min_mean_max(comm, xxg   (k,:), mg   )
      if(allocated(plck)) call imp_min_mean_max(comm, xxplck(k,:), mplck)
      if(allocated(plck)) then
        if(myid.eq.0) print *,k,cstr(ftoa(mkabs),'blue'), cstr(ftoa(mksca),'red'), cstr(ftoa(mg),'blue'), cstr(ftoa(mplck),'red')
      else
        if(myid.eq.0) print *,k,cstr(ftoa(mkabs),'blue'), cstr(ftoa(mksca),'red'), cstr(ftoa(mg),'blue')
      endif
      enddo
      if(myid.eq.0) print *,cstr('---------------------------------------------------------------------------------------', 'blue')
      if(allocated(albedo)) call imp_min_mean_max(comm, xxalbedo(i1,:), malbedo)
      if(myid.eq.0) print *,cstr('Surface Albedo (min,mean,max) '//ftoa(malbedo), 'blue')

      call restore_col_vec(albedo, xalbedo, xxalbedo)
      call restore_col_vec(kabs  , xkabs  , xxkabs  )
      call restore_col_vec(ksca  , xksca  , xxksca  )
      call restore_col_vec(g     , xg     , xxg     )
      call restore_col_vec(plck  , xplck  , xxplck  )
      contains
        subroutine get_col_vec(vec, N, xv, xxv)
          type(tVec), allocatable, intent(in) :: vec
          integer(iintegers), intent(in) :: N ! number of entries per column/stride
          real(ireals), intent(inout), pointer :: xv(:), xxv(:,:)
          integer(iintegers) :: vecsize, Ncol
          if(allocated(vec)) then
            call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
            Ncol = vecsize / N
            call CHKERR(int(Ncol*N-vecsize, mpiint), 'stride size of vec not as expected. '//new_line('')// &
              'Given stride size '//itoa(N)//new_line('')// &
              'Local vector size '//itoa(vecsize)//new_line('')// &
              'expected number of columns/strides '//itoa(Ncol)//' would result in '//itoa(Ncol*N)//'entries')

            call VecGetArrayF90(vec, xv, ierr); call CHKERR(ierr)
            xxv(1:N, 1:Ncol) => xv(:)
          endif
        end subroutine
        subroutine restore_col_vec(vec, xv, xxv)
          type(tVec), allocatable, intent(in) :: vec
          real(ireals), intent(inout), pointer :: xv(:), xxv(:,:)
          if(associated(xxv)) nullify(xxv)
          if(allocated(vec)) call VecRestoreArrayF90(vec, xv, ierr); call CHKERR(ierr)
        end subroutine
    end subroutine

    subroutine create_edir_src_vec(solver, plex, edirdm, E0, kabs, ksca, g, sundir, srcVec)
      class(t_plex_solver), allocatable, intent(in) :: solver
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
      integer(iintegers) :: iface, icell, idof, numdof

      integer(iintegers), pointer :: xx_v(:)
      integer(iintegers), pointer :: cell_support(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: area, mu_top, face_normal(3)

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
          area = geoms(i1+geom_offset) / real(solver%dirtop%area_divider, ireals)

          call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)

          if(is_solar_src(face_normal, sundir)) then
            mu_top = -dot_product(sundir, -face_normal)
            call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
            call PetscSectionGetDof(edirSection, iface, numDof, ierr); call CHKERR(ierr)
            do idof = 1, numDof
              xv(voff+idof) = E0 * area * mu_top
            enddo
          endif
        enddo
        call ISRestoreIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)
      endif ! TOA boundary ids

      call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      !if(.False.) &
        call set_sidewards_direct_fluxes()

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
          real(ireals) :: tau, dz, transport_this_cell, mu_top, mu_side

          type(tVec) :: vedir
          real(ireals), pointer :: xedir(:)

          call DMGetGlobalVector(edirdm, vedir, ierr); call CHKERR(ierr)

          call DMGetStratumIS(edirdm, 'DomainBoundary', SIDEFACE, IS_side_faces, ierr); call CHKERR(ierr)
          if (IS_side_faces.eq.PETSC_NULL_IS) then ! dont have side boundary faces
          else

            call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
            call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
            call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)


            ! First we compute the edir vec on all columns that are on the domain side and the sun is shining on it
            call VecSet(vedir, -one, ierr); call CHKERR(ierr)
            call VecGetArrayF90(vedir, xedir, ierr); call CHKERR(ierr)

            call ISGetIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)

            do o = 1, size(iside_faces)
              iface_side = iside_faces(o)

              call DMPlexGetSupport(edirdm, iface_side, cell_support, ierr); call CHKERR(ierr)
              icell = cell_support(1)
              call DMPlexRestoreSupport(edirdm, iface_side, cell_support, ierr); call CHKERR(ierr)

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
              if(.not.is_solar_src(face_normal, sundir)) cycle ! if sun is not shining on top of TOA face, skip it
              mu_top = -dot_product(sundir, -face_normal)

              ! compute direct radiation in 1D fashion
              call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)
              xedir(i1+topoffset) = E0

              do k = 1, plex%Nlay
                icell = vert_cell_idx(k)

                call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
                dz = geoms(i1+geom_offset)

                call DMPlexGetCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
                topface = faces_of_cell(1)
                botface = faces_of_cell(2)
                call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)

                call get_inward_face_normal(topface, icell, geomSection, geoms, face_normal)
                mu_top = -dot_product(sundir, -face_normal)

                tau = (xkabs((i1+icell)) + xksca((i1+icell))) * dz / mu_top

                call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)
                call PetscSectionGetOffset(edirsection, botface, botoffset, ierr); call CHKERR(ierr)

                xedir(i1+botoffset) = xedir(i1+topoffset) * exp(-tau)

                !print *,icell,'edir', k, xedir(i1+botoffset)
              enddo
            enddo

            ! Now we go ahead and try to estimate the fluxes along the side boundary

            call VecGetArrayF90(localVec, xv, ierr); call CHKERR(ierr)


            do o = 1, size(iside_faces)
              iface_side = iside_faces(o)
              call PetscSectionGetDof(edirSection, iface_side, numDof, ierr); call CHKERR(ierr)
              if(numDof.lt.i1) cycle

              call DMPlexGetSupport(edirdm, iface_side, cell_support, ierr); call CHKERR(ierr)
              icell = cell_support(1)
              call DMPlexRestoreSupport(edirdm, iface_side, cell_support, ierr); call CHKERR(ierr)

              call get_inward_face_normal(iface_side, icell, geomSection, geoms, side_face_normal)
              if(.not.is_solar_src(side_face_normal, sundir)) cycle ! if sun is not shining on this side face, skip this

              call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
              dz = geoms(i1+geom_offset)

              call DMPlexGetCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
              topface = faces_of_cell(1)
              botface = faces_of_cell(2)
              call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)

              call PetscSectionGetOffset(edirsection, topface, topoffset, ierr); call CHKERR(ierr)

              call PetscSectionGetFieldOffset(geomSection, iface_side, i2, geom_offset, ierr); call CHKERR(ierr)
              area = geoms(i1+geom_offset) / real(solver%dirtop%area_divider, ireals)

              call get_inward_face_normal(topface, icell, geomSection, geoms, face_normal)
              mu_top = -dot_product(sundir, -face_normal)

              tau = (xkabs((i1+icell)) + xksca((i1+icell))) * dz/2 / mu_top

              mu_side = -dot_product(sundir, -side_face_normal)

              transport_this_cell = exp(-tau) * mu_side

              call PetscSectionGetOffset(edirsection, iface_side, voff, ierr); call CHKERR(ierr)
              do idof=1,numDof
                xv(voff+idof) = xedir(i1+topoffset) * area * transport_this_cell
              enddo
              !print *,icell,'edir top', xedir(i1+topoffset), 'sidesrc', xv(i1+voff)
            enddo

            call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

            call ISRestoreIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)

            call VecRestoreArrayF90(vedir, xedir, ierr); call CHKERR(ierr)

            call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
            call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
            call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
          endif ! TOA boundary ids
          call DMRestoreGlobalVector(edirdm, vedir, ierr); call CHKERR(ierr)
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
        plckdm, plckVec, &
        edirdm, edirVec)
      class(t_plex_solver), allocatable, intent(in) :: solver
      type(t_plexgrid), intent(in) :: plex
      class(t_optprop_wedge), intent(in) :: OPP
      type(tDM), allocatable, intent(in) :: ediffdm
      type(tVec), allocatable, intent(in) :: kabs, ksca, g ! cell1_dm
      type(tVec), allocatable, intent(in) :: albedo ! srfc_boundary_dm
      logical, intent(in) :: lthermal, lsolar
      type(tVec), allocatable, intent(inout) :: srcVec

      type(tDM), allocatable, intent(in), optional :: plckdm
      type(tVec), allocatable, intent(in), optional :: plckVec ! horizface1_dm

      type(tDM), allocatable, intent(in), optional :: edirdm
      type(tVec), allocatable, intent(in), optional :: edirVec

      type(tVec) :: lsrcVec

      type(tPetscSection) :: geomSection, wedgeSection, ediffSection, srfcSection, plckSection
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
      call DMGetSection(plckdm, plckSection, ierr); CHKERRQ(ierr)
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
          integer(iintegers) :: icell, isrc, isrc_face, idst, icol, in_dof, isrc_side, numSrc
          integer(iintegers) :: bmcsrcdof, bmcsrcside, dof_src_offset

          integer(iintegers) :: face_plex2bmc(5), geom_offset
          integer(iintegers) :: wedge_offset1, wedge_offset2, wedge_offset3, wedge_offset4
          integer(iintegers) :: dir_plex2bmc(solver%dirdof), diff_plex2bmc(solver%diffdof)
          real(ireals), pointer :: xedir(:)
          real(ireals) :: param_phi, param_theta, area_bot, area_top, dz, incsol
          !real(ireals) :: dir2diff(solver%diffdof/2)
          logical :: lsrc(5)

          real(irealLUT) :: coeff(solver%dirdof*(solver%diffdof/2)), sumcoeff
          integer(iintegers), pointer :: xinoutdof(:)

          integer(mpiint) :: myid, cerr

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

            call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
            param_phi   = wedgeorient(wedge_offset1+i3)
            param_theta = wedgeorient(wedge_offset1+i4)

            call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(wedgeSection, icell, i2, wedge_offset3, ierr); call CHKERR(ierr)
            do isrc = 1, size(faces_of_cell)
              face_plex2bmc(int(wedgeorient(wedge_offset2+isrc), iintegers)) = isrc
              lsrc(isrc) = wedgeorient(wedge_offset3+isrc) .gt. zero
            enddo
            call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)
            call face_idx_to_bmc_idx(solver%dirtop, solver%dirside, face_plex2bmc, dir_plex2bmc)
            !print *,'icell', icell, 'face_plex2bmc', face_plex2bmc, 'dir_plex2bmc', dir_plex2bmc

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            call PetscLogEventBegin(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)
            !print *,'get coeff dir2diff', shape(coeff), ':', solver%dirdof, solver%diffdof
            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
              aspect_zx=wedgeorient(wedge_offset1+i5), &
              Cx=wedgeorient(wedge_offset4+i1), &
              Cy=wedgeorient(wedge_offset4+i2), &
              ldir=.False., coeff=coeff, ierr=cerr, &
              angles=[real(param_phi, irealLUT), real(param_theta, irealLUT)])

            !if(cerr.eq.OPP_TINYASPECT_RETCODE) call CHKERR(1_mpiint, 'tiny aspect !DEBUG')
            if(cerr.eq.OPP_1D_RETCODE) then
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
              area_top = geoms(i1+geom_offset)
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
              area_bot = geoms(i1+geom_offset)
              !DEBUG coeff(4*8+1) = coeff(4*8+1) * area_bot / area_top
            endif

            call PetscLogEventEnd(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)

            if(ldebug) sumcoeff=0
            in_dof = 0 ! this counts plex dofs from 1 to solver%dirdof

            do isrc_side = 1, size(faces_of_cell)
              isrc_face = faces_of_cell(isrc_side)
                call PetscSectionGetDof(edirSection, isrc_face, numSrc, ierr); call CHKERR(ierr)
                do isrc = 1, numSrc
                  in_dof = in_dof+1

                  if(.not.lsrc(isrc_side)) cycle

                  bmcsrcdof = dir_plex2bmc(in_dof)
                  call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcsrcdof, &
                    bmcsrcside, dof_src_offset)

                  call PetscSectionGetOffset(edirSection, isrc_face, icol, ierr); call CHKERR(ierr)
                  icol = icol + dof_src_offset
                  incsol = xedir(i1+icol)

                  associate( dir2diff => coeff(bmcsrcdof:size(coeff):solver%dirdof) ) ! dir2diff in src ordering
                    !print *,'icell '//itoa(icell)//' in_dof '//itoa(in_dof)//' bmsrc '//itoa(bmcsrcdof)// &
                    !  'dir2diff', dir2diff, ':', cstr(ftoa(sum(dir2diff)), 'green')

                    do idst = 1, size(outgoing_offsets)
                      if(dir2diff(diff_plex2bmc(idst)).gt.zero) then
                        if(outgoing_offsets(idst).lt.0) call CHKERR(1_mpiint, 'does this happen?')
                      endif
                      if(outgoing_offsets(idst).lt.0) cycle

                      xb(i1+outgoing_offsets(idst)) = xb(i1+outgoing_offsets(idst)) + &
                        incsol * dir2diff(diff_plex2bmc(idst))

                      if(ldebug) sumcoeff = sumcoeff + dir2diff(diff_plex2bmc(idst))
                      if(.False.) then
                        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
                        area_top = geoms(i1+geom_offset)

                        print *,cstr('Setting diffsrc for ', 'blue')// &
                          ' icell '//itoa(icell)//&
                          ' dst '//itoa(outgoing_offsets(idst))//&
                          ' src-face '//itoa(isrc_side)//' edir_off '//itoa(icol)// &
                          ' edir='//ftoa(incsol/area_top)//  &
                          ' * idst ('//itoa(idst)//')'// &
                          ' p2bmc '//itoa(diff_plex2bmc(idst))// &
                          ' coeff:'//ftoa(dir2diff(diff_plex2bmc(idst)))// &
                          ' xb '//ftoa(xb(i1+outgoing_offsets(idst))/area_top)// &
                          ' ('//ftoa(xb(i1+outgoing_offsets(idst)))//')'
                      endif

                    enddo ! outgoing_offsets
                  end associate
!                        todo: sum up the coeffs and check that this is one-sum(dir2dir) if kabs=0, i.e. make sure that the sum of
!                        coeffs that have been used is the sum of dir2diff, this energy conservation check would be useful in dir2dir
!                        and diff2diff as well, just to make sure that we dont neglect any lut coeff contributions
                enddo ! numSrc
            enddo ! srcfaces
            if(ldebug) then
              if(.not.approx(sumcoeff, sum(coeff), sqrt(epsilon(coeff)) )) then
                print *,'faces_of_cell', faces_of_cell, ' lsrc? ', lsrc
                do isrc=1,solver%dirdof
                  call CHKWARN(1_mpiint, 'dir2diff src '//itoa(isrc)//' : '//&
                    ftoa(coeff(isrc:size(coeff):solver%dirdof)))
                enddo
                call CHKERR(1_mpiint, 'have we used all coefficients? '// &
                  'used '//ftoa(sumcoeff)//' out of '//ftoa(sum(coeff)))
              endif
            endif


            call DMPlexRestoreCone(edirdm, icell, faces_of_cell, ierr); call CHKERR(ierr)
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
          integer(iintegers) :: wedge_offset1, wedge_offset2, wedge_offset4, geom_offset, plck_offset

          integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)
          integer(iintegers) :: i, j, icell, icol, iface, idof, numDof
          real(irealLUT) :: coeff((solver%diffdof/2)**2) ! coefficients for each src=[1..8] and dst[1..8]

          integer(iintegers), pointer :: xinoutdof(:)
          real(ireals), pointer :: xplanck(:)
          real(ireals) :: diff2diff(8), dz
          real(ireals) :: emissivity, b0, b1, btop, bbot, bside, Beff

          call ISGetIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
          call VecGetArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)

          do icell = plex%cStart, plex%cEnd-1

            call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
            call get_in_out_dof_offsets(solver%IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)

            call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)

            do i = 1, size(faces_of_cell)
              face_plex2bmc(int(wedgeorient(wedge_offset2+i), iintegers)) = i
            enddo
            call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
              aspect_zx=wedgeorient(wedge_offset1+i5), &
              Cx=wedgeorient(wedge_offset4+i1), &
              Cy=wedgeorient(wedge_offset4+i2), &
              ldir=.False., coeff=coeff, ierr=ierr)

            !DEBUG have to check this again if(ierr.eq.OPP_1D_RETCODE) then
            !DEBUG have to check this again   !call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
            !DEBUG have to check this again   !area_top = geoms(i1+geom_offset)
            !DEBUG have to check this again   !call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
            !DEBUG have to check this again   !area_bot = geoms(i1+geom_offset)
            !DEBUG have to check this again   !if(area_top.gt.area_bot) then
            !DEBUG have to check this again   !  ! send overhanging energy to side faces
            !DEBUG have to check this again   !  coeff(7*8+[2,4,6]) = coeff(7*8+[2,4,6]) + coeff(7*8+1) * real(one - area_bot / area_top, irealLUT) / i3
            !DEBUG have to check this again   !  ! reduce transmission bc receiver is smaller than top face
            !DEBUG have to check this again   !  coeff(7*8+1) = coeff(7*8+1) * real(area_bot / area_top, irealLUT)
            !DEBUG have to check this again   !else
            !DEBUG have to check this again   !  ! send overhanging energy to side faces
            !DEBUG have to check this again   !  coeff([3,5,7]) = coeff([3,5,7]) + coeff(8) * real(one - area_top / area_bot, irealLUT) / i3
            !DEBUG have to check this again   !  ! transmission from bot to top face
            !DEBUG have to check this again   !  coeff(8) = coeff(8) * real(area_top / area_bot, irealLUT)
            !DEBUG have to check this again   !endif
            !DEBUG have to check this again endif
            call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

            call PetscSectionGetFieldOffset(plckSection, faces_of_cell(1), i0, plck_offset, ierr); call CHKERR(ierr)
            b0 = xplanck(i1+plck_offset)  * pi ! top planck value
            b1 = xplanck(i1+plck_offset+1)* pi ! bot planck value
            call B_eff(b1, b0, xkabs(i1+icell), btop)
            call B_eff(b0, b1, xkabs(i1+icell), bbot)
            bside = (btop+bbot)/2

            i = 1
            do j = 1, size(faces_of_cell)
              iface = faces_of_cell(j)
              call PetscSectionGetDof(ediffSection, iface, numDof, ierr); call CHKERR(ierr)

              select case(j)
                case(1) ! top face
                  Beff = btop
                case(2) ! bot face
                  Beff = bbot
                case default !side faces
                  Beff = bside
              end select

              do idof = 1, numDof/2
                icol = outgoing_offsets(i)

                diff2diff = coeff(diff_plex2bmc(i): size(coeff): i8)

                emissivity = max(zero, one - sum(diff2diff))

                xsrc(i1+icol) = Beff * emissivity / real(numDof/2, ireals)
                i = i+1
              enddo
            enddo
            call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
          enddo ! icell
          call VecRestoreArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)
          call ISRestoreIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
        end subroutine

        subroutine thermal_srfc_emission(xsrc)
          real(ireals), pointer :: xplanck(:), xsrc(:)
          real(ireals), pointer :: xalbedo(:)
          type(tIS) :: srfc_ids
          integer(iintegers), pointer :: xi(:)
          integer(iintegers) :: i, ke1, iface, offset_Eup, offset_srfc, plck_offset

          ke1 = solver%plex%Nlay+1

          call VecGetArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)
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

              call PetscSectionGetFieldOffset(plckSection, iface, i0, plck_offset, ierr); call CHKERR(ierr)
              xsrc(i1+offset_Eup) = xplanck(i1+plck_offset) * (one-xalbedo(i1+offset_srfc)) * pi
            enddo
            call ISRestoreIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
          endif

          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
          call VecRestoreArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)
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

      real(ireals),parameter :: rtol=1e-6_ireals, rel_atol=1e-6_ireals
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
        atol = rel_atol*real(Nrows_global, ireals)
        !call imp_allreduce_min(comm, rel_atol*real(Nrows_global, ireals), atol)
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

      call hegedus_trick(ksp, b, x)
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
          logical :: lhandle_diverged, lflg

          call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)
          if(reason.le.0) then
            call PetscObjectViewFromOptions(b  , PETSC_NULL_VEC, '-show_diverged_b', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(x  , PETSC_NULL_VEC, '-show_diverged_x', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(A  , PETSC_NULL_MAT, '-show_diverged_A', ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(ksp, PETSC_NULL_KSP, '-show_diverged_ksp', ierr); call CHKERR(ierr)
          endif

          lhandle_diverged = .False.
          call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-plexrt_handle_diverged_solver",&
            lhandle_diverged, lflg, ierr) ;call CHKERR(ierr)
          if(.not.lhandle_diverged) return

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

    subroutine scale_facevec(plex, face_dm, top_dof, side_dof, globalfaceVec, lW_to_Wm2)
    type(t_plexgrid), intent(in) :: plex
    type(tDM), intent(in) :: face_dm
    type(t_dof), intent(in) :: top_dof, side_dof
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
      if(plex%ltopfacepos(iface)) then
        area = area / real(top_dof%area_divider, ireals)
      else
        area = area / real(side_dof%area_divider, ireals)
      endif

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
    class(t_optprop_wedge), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g
    real(ireals), intent(in) :: sundir(3)
    type(tMat), allocatable, intent(inout) :: A

    type(tPetscSection) :: sec

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: icell, iface, irow, icol, idst, isrc
    integer(iintegers) :: bmcsrcdof, bmcdstdof, bmcsrcside, bmcdstside, dof_src_offset, dof_dst_offset
    integer(iintegers) :: idst_side, isrc_side
    integer(iintegers) :: fStart, fEnd
    integer(iintegers) :: in_dof, out_dof
    real(ireals) :: c

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: geom_offset, wedge_offset1, wedge_offset2, wedge_offset3, wedge_offset4

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: plex2bmc(solver%dirdof)

    real(ireals) :: dz, area_top, area_bot

    real(irealLUT) :: dir2dir(solver%dirdof)
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)
    integer(iintegers) :: numSrc, numDst

    real(ireals) :: aspect_zx, param_phi, param_theta
    real(irealLUT) :: coeff(solver%dirdof**2), sumcoeff ! coefficients for each src=[1..5] and dst[1..5]

    logical :: lflg, ldestroy_mat, l1d

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    !if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...'

    if(.not.allocated(plex%geom_dm)) call CHKERR(1_mpiint, 'geom_dm has to allocated in order to create an Edir Matrix')
    if(.not.allocated(plex%edir_dm)) call CHKERR(1_mpiint, 'edir_dm has to allocated in order to create an Edir Matrix')
    if(.not.allocated(kabs  )) call CHKERR(1_mpiint, 'kabs   has to be allocated')
    if(.not.allocated(ksca  )) call CHKERR(1_mpiint, 'ksca   has to be allocated')
    if(.not.allocated(g     )) call CHKERR(1_mpiint, 'g      has to be allocated')
    if(.not.allocated(plex%wedge_orientation_dm)) &
      call CHKERR(1_mpiint, 'wedge_orientation_dm has to allocated in order to create an Edir Matrix')

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

    if(.False.) call set_side_incoming_boundary_condition()

    ! Set Diagonal Entries
    call DMPlexGetDepthStratum(plex%edir_dm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
    do iface = fStart, fEnd-1
      call PetscSectionGetDof(sec, iface, numDst, ierr); call CHKERR(ierr)
      call PetscSectionGetOffset(sec, iface, irow, ierr); call CHKERR(ierr)
      do idst = 0, numDst-1
        call MatSetValuesLocal(A, i1, irow+idst, i1, irow+idst, [one], INSERT_VALUES, ierr); call CHKERR(ierr)
      enddo
    enddo

    do icell = plex%cStart, plex%cEnd-1
      call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i2, wedge_offset3, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      !do iface=1,size(faces_of_cell)
      !  call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr)
      !  print *,'faces of cell',iface, faces_of_cell(iface),':: offset',icol
      !enddo

      param_phi   = wedgeorient(wedge_offset1+i3)
      param_theta = wedgeorient(wedge_offset1+i4)
      aspect_zx   = wedgeorient(wedge_offset1+i5)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset2+iface), iintegers)) = iface
        lsrc(iface) = wedgeorient(wedge_offset3+iface) .gt. zero
      enddo
      call face_idx_to_bmc_idx(solver%dirtop, solver%dirside, face_plex2bmc, plex2bmc)

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      call PetscLogEventBegin(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
        aspect_zx=wedgeorient(wedge_offset1+i5), &
        Cx=wedgeorient(wedge_offset4+i1), &
        Cy=wedgeorient(wedge_offset4+i2), &
        ldir=.True., coeff=coeff, ierr=ierr, &
        angles=[real(param_phi, irealLUT), real(param_theta, irealLUT)])
      l1d = ierr.eq.OPP_1D_RETCODE

      !print *,'tau',xkabs(i1+icell)+xksca(i1+icell),aspect_zx, ':', param_phi, param_theta, 'l1d', l1d
      !print *,'tau',xkabs(i1+icell)+xksca(i1+icell),dz, cos(zenith), ':', exp(-(xkabs(i1+icell)+xksca(i1+icell))*dz / cos(zenith))

      if(ldebug) then
        do out_dof = 1, solver%dirdof
          bmcdstdof = plex2bmc(out_dof)
          dir2dir = coeff(1+(bmcdstdof-1)*solver%dirdof:bmcdstdof*solver%dirdof)
          !write(*, FMT='("out_dof " I2  " bmc_dst " I2 " : " 18(f10.5))') out_dof, bmcdstdof, dir2dir
          !print *,out_dof, 'bmcdst', bmcdstdof, dir2dir
        enddo
        do in_dof = 1, solver%dirdof
          bmcsrcdof = plex2bmc(in_dof)
          dir2dir = coeff(bmcsrcdof:size(coeff):solver%dirdof)
          !print *, icell, in_dof, 'dir2dir', dir2dir, ':', cstr(ftoa(sum(dir2dir)),'green')
          if(sum(dir2dir).gt.1+sqrt(epsilon(dir2dir))*100) then
            print *,in_dof,': bmcface', bmcsrcdof, 'dir2dir gt one', dir2dir,':',sum(dir2dir)
            call CHKERR(1_mpiint, 'energy conservation violated! '//ftoa(sum(dir2dir)))
          endif
        enddo
        !print *,'-----'
      endif

      if(l1d) then
        ! for the case of 1D spherical radiative transfer,
        ! need to consider the change in area between upper and lower face
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
        area_top = geoms(i1+geom_offset)
        call PetscSectionGetFieldOffset(geomSection, faces_of_cell(2), i2, geom_offset, ierr); call CHKERR(ierr)
        area_bot = geoms(i1+geom_offset)

        ! super compensate direct radiation
        !DEBUG coeff(21) = coeff(21) * (area_top / area_bot)
      endif
      call PetscLogEventEnd(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)

      if(ldebug) sumcoeff=0
      in_dof = 0 ! this counts plex dofs from 1 to solver%dirdof
      do isrc_side = 1, size(faces_of_cell)
        call PetscSectionGetDof(sec, faces_of_cell(isrc_side), numSrc, ierr); call CHKERR(ierr)
        do isrc = 1, numSrc
          in_dof = in_dof+1
          !if(.not.lsrc(isrc_side)) cycle

          bmcsrcdof = plex2bmc(in_dof)
          call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcsrcdof, &
            bmcsrcside, dof_src_offset)

          call PetscSectionGetOffset(sec, faces_of_cell(isrc_side), icol, ierr); call CHKERR(ierr)
          icol = icol + dof_src_offset

          dir2dir = coeff(bmcsrcdof:size(coeff):solver%dirdof)
          if(sum(dir2dir).gt.1+sqrt(epsilon(dir2dir))) dir2dir = dir2dir / sum(dir2dir)

          out_dof = 0
          do idst_side = 1, size(faces_of_cell)
            call PetscSectionGetDof(sec, faces_of_cell(idst_side), numDst, ierr); call CHKERR(ierr)
            do idst = 1, numDst
              out_dof = out_dof+1

              !if(lsrc(idst_side)) cycle
              bmcdstdof = plex2bmc(out_dof)
              call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcdstdof, &
                bmcdstside, dof_dst_offset)

              call PetscSectionGetOffset(sec, faces_of_cell(idst_side), irow, ierr); call CHKERR(ierr)
              irow = irow + dof_dst_offset

              c = -real(dir2dir(bmcdstdof), kind(c))

              if(c.lt.zero .and. .not.lsrc(isrc_side)) then
                call CHKERR(1_mpiint, 'found transport coeff but I thought this incoming side is not a designated src face')
              endif
              if(c.le.-1e-6_ireals.and..not.l1d.and.param_theta.gt.epsilon(zero)) then
                ierr = 0
                if(.not.lsrc(isrc_side)) ierr = 1
                if(     lsrc(idst_side)) ierr = 2
                if(ierr.ne.0) then
                  print *,'numsrc', numSrc, 'kindex cell', plex%zindex(icell)
                  print *,'faces of cell',faces_of_cell
                  print *,'plexface 2 bmc', face_plex2bmc
                  print *,'lsrc(plex)', lsrc
                  print *,'param_phi', param_phi, 'param_theta', param_theta
                  print *,'cell   ('//itoa(icell)//') '// &
                    'incoming dof '//itoa(in_dof)//' -> out_dof '//itoa(out_dof)//' '// &
                    'src_side  '//itoa(isrc_side)//' '// &
                    ' ('//itoa(faces_of_cell(isrc_side))//') '// &
                    'dst_side  '//itoa(idst_side)//' '// &
                    ' ('//itoa(faces_of_cell(idst_side))//') '// &
                    'bmc_src_dof    ('//itoa(bmcsrcdof)//') '// &
                    'bmc_src_side   ('//itoa(bmcsrcside)//') '// &
                    'bmc_dst_side   ('//itoa(bmcdstside)//') '// &
                    'src_offset   ('//itoa(dof_src_offset)//') '// &
                    'dst_offset   ('//itoa(dof_dst_offset)//') '// &
                    'col -> row    '//itoa(icol)//' -> '//itoa(irow)//'   = ', c
                  if(ierr.eq.1) call CHKERR(ierr, 'found coeff, but thought this incoming side is not a src face')
                  if(ierr.eq.2) call CHKERR(ierr, 'have coeff, but the target dst is marked as being a src face')
                endif
              endif
              if(c.ge.zero) cycle

              if(ldebug) sumcoeff = sumcoeff + dir2dir(bmcdstdof)
              call MatSetValuesLocal(A, i1, irow, i1, icol, c, INSERT_VALUES, ierr); call CHKERR(ierr)

            enddo
          enddo

        enddo
      enddo

      if(ldebug) then
        if(.not.approx(sumcoeff, sum(coeff), sqrt(epsilon(coeff)) )) then
          print *,'faces_of_cell', faces_of_cell, ' lsrc? ', lsrc
          do isrc=1,solver%dirdof
          call CHKWARN(1_mpiint, 'dir2dir src '//itoa(isrc)//' : '//&
            ftoa(coeff(isrc:size(coeff):solver%dirdof)))
          enddo
          call CHKERR(1_mpiint, 'have we used all coefficients? '// &
            'used '//ftoa(sumcoeff)//' out of '//ftoa(sum(coeff)))
        endif
      endif

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
        integer(iintegers) :: o, i, icell, iface_top, iface_side
        real(ireals) :: inv_sundir(3), zenith, azimuth
        real(ireals) :: top_face_normal(3), side_face_normal(3), Mrot(3,3), Cx, Cy
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
            Mrot = rotation_matrix_around_axis_vec(deg2rad(180._ireals), top_face_normal)
            inv_sundir = matmul(sundir, Mrot)

            !print *,sundir, norm2(sundir), 'inv_sundir',inv_sundir, norm2(inv_sundir), &
            !  'angle between', rad2deg( angle_between_two_vec(inv_sundir, sundir) )
            !call CHKERR(1_mpiint, 'DEBUG')
            !print *,'faces_of_cell', faces_of_cell, 'sideface', iface_side

            call compute_local_wedge_ordering(plex, icell, &
              geomSection, geoms, inv_sundir, &
              lsrc, zenith, azimuth, param_phi, param_theta, &
              Cx, Cy, aspect_zx, &
              upper_face, bottom_face, &
              base_face, left_face, right_face)

            forient = [upper_face, base_face, left_face, right_face, bottom_face]
            do i=1, size(faces_of_cell)
              face_plex2bmc(forient(i)) = i
            enddo
            call face_idx_to_bmc_idx(solver%dirtop, solver%dirside, face_plex2bmc, plex2bmc)

            call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
            dz = geoms(i1+geom_offset)

            !print *,'az', rad2deg(azimuth), 'zenith', rad2deg(zenith), 'aspect_zx', aspect_zx
            !print *,'upper_face, bottom_face', upper_face, bottom_face
            !print *,'base_face, left_face, right_face', base_face, left_face, right_face
            !print *,'coords2d', coords_2d
            !print *,'face_plex2bmc', face_plex2bmc

            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
              aspect_zx, Cx, Cy, .True., coeff, ierr, &
              angles=[real(param_phi,irealLUT), real(param_theta,irealLUT)])

            in_dof = 0 ! this counts plex dofs from 1 to solver%dirdof
            do isrc_side = 1, size(faces_of_cell)
              call PetscSectionGetDof(sec, faces_of_cell(isrc_side), numSrc, ierr); call CHKERR(ierr)
              do isrc = 1, numSrc
                in_dof = in_dof+1
                if(.not.lsrc(isrc_side)) cycle

                bmcsrcdof = plex2bmc(in_dof)
                call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcsrcdof, &
                  bmcsrcside, dof_src_offset)

                call PetscSectionGetOffset(sec, faces_of_cell(isrc_side), icol, ierr); call CHKERR(ierr)
                icol = icol + dof_src_offset

                dir2dir = coeff(bmcsrcdof:size(coeff):solver%dirdof)

                out_dof = 0
                do idst_side = 1, size(faces_of_cell)
                  call PetscSectionGetDof(sec, faces_of_cell(idst_side), numDst, ierr); call CHKERR(ierr)
                  do idst = 1, numDst
                    out_dof = out_dof+1

                    if(faces_of_cell(idst_side).ne.iface_side) cycle ! only the respective side face needs the boundary condition
                    if(lsrc(idst_side)) cycle
                    bmcdstdof = plex2bmc(out_dof)
                    call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcdstdof, &
                      bmcdstside, dof_dst_offset)

                    c = -real(dir2dir(bmcdstdof), kind(c))

                    call PetscSectionGetOffset(sec, faces_of_cell(idst_side), irow, ierr); call CHKERR(ierr)
                    irow = irow + dof_dst_offset

                    !print *,'cell   ('//itoa(icell)//') '// &
                    !  'incoming dof -> out_dof   '//itoa(in_dof)//' -> '//itoa(out_dof)//' '// &
                    !  'src_side  '//itoa(isrc_side)//' '// &
                    !  ' ('//itoa(faces_of_cell(isrc_side))//') '// &
                    !  'dst_side  '//itoa(idst_side)//' '// &
                    !  ' ('//itoa(faces_of_cell(idst_side))//') '// &
                    !  'bmc_src_side   ('//itoa(bmcsrcside)//') '// &
                    !  'bmc_dst_side   ('//itoa(bmcdstside)//') '// &
                    !  'src_offset   ('//itoa(dof_src_offset)//') '// &
                    !  'dst_offset   ('//itoa(dof_dst_offset)//') '// &
                    !  'col -> row    '//itoa(icol)//' -> '//itoa(irow)//'   = ', c
                    call MatSetValuesLocal(A, i1, irow, i1, icol, c, INSERT_VALUES, ierr); call CHKERR(ierr)


                  enddo
                enddo
              enddo
            enddo

          endif ! boundary side face is sunlit
        enddo
      end subroutine

  end subroutine

  subroutine face_idx_to_bmc_idx(t_top_dof, t_side_dof, face_plex2bmc, plex2bmc)
    type(t_dof), intent(in) :: t_top_dof, t_side_dof
    integer(iintegers), intent(in) :: face_plex2bmc(:)  ! mapping from faces_of_cell to bmc_faces ordering (dim=5)
    integer(iintegers), intent(out) :: plex2bmc(:) ! mapping from dof streams in plex ordering to bmc dof (dim=dof)
    integer(iintegers) :: i, j, idof

    ! basic mapping for 5 dof would be e.g.:
    ! top face  -> 1
    ! base_face -> 2
    ! left_face -> 3
    ! right_face-> 4
    ! bot face  -> 5

    ! or for 18 dof
    ! top face  -> 1,2,3
    ! base_face -> 4,5,6,7
    ! left_face -> 8,9,10,11
    ! right_face-> 12,13,14,15
    ! bot face  -> 16,17,18

    ! however, we also need to consider the side face permuations

    j = 1
    do i = 1, size(face_plex2bmc)
      select case(face_plex2bmc(i)) ! case switch on bmc indices
      case(1)
        do idof=1,t_top_dof%dof
          plex2bmc(j) = idof
          j = j+1
        enddo
      case(5)
        do idof=t_top_dof%dof + 3*t_side_dof%dof+1, 2*t_top_dof%dof + 3*t_side_dof%dof
          plex2bmc(j) = idof
          j = j+1
        enddo
      case(2)
        do idof=t_top_dof%dof + 1, t_top_dof%dof + t_side_dof%dof
          plex2bmc(j) = idof
          j = j+1
        enddo
      case(3)
        do idof=t_top_dof%dof + t_side_dof%dof + 1, t_top_dof%dof + 2*t_side_dof%dof
          plex2bmc(j) = idof
          j = j+1
        enddo
      case(4)
        do idof=t_top_dof%dof + 2*t_side_dof%dof + 1, t_top_dof%dof + 3*t_side_dof%dof
          plex2bmc(j) = idof
          j = j+1
        enddo
      end select
    enddo
  end subroutine

  ! gives the bmc side and bmc offset on this side for a given idof number in bmc ordering
  ! e.g. for dir_18 gives
  ! top face first dof  :: 1 -> (1,0)
  ! top face second dof :: 2 -> (1,1)
  ! base face first dof :: 4 -> (2,0)
  ! bot face second dof :: 17 -> (5,1)
  subroutine get_side_and_offset_from_total_bmc_dofs(t_top_dof, t_side_dof, idof, side, dof_offset)
    type(t_dof), intent(in) :: t_top_dof, t_side_dof
    integer(iintegers), intent(in) :: idof
    integer(iintegers), intent(out) :: side, dof_offset
    integer(iintegers) :: k
    side = 1
    dof_offset = idof-1
    if(dof_offset.lt.t_top_dof%dof) return
    dof_offset = dof_offset - t_top_dof%dof

    do k=1,3
      side = side +1
      if(dof_offset.lt.t_side_dof%dof) return
      dof_offset = dof_offset - t_side_dof%dof
    enddo
    side = side+1
  end subroutine

  subroutine create_ediff_mat(solver, plex, OPP, kabs, ksca, g, albedo, A)
    class(t_plex_solver), intent(in) :: solver
    type(t_plexgrid), intent(in) :: plex
    class(t_optprop_wedge), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g ! cell1_dm
    type(tVec), allocatable, intent(in) :: albedo        ! srfc_boundary_dm
    type(tMat), allocatable, intent(inout) :: A

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: i, j, icell, iface
    integer(iintegers), allocatable :: irows(:), icols(:)
    real(ireals), allocatable :: c(:)

    type(tPetscSection) :: ediffSection, geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: wedge_offset1, wedge_offset2, wedge_offset4, geom_offset

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: diff_plex2bmc(solver%diffdof/2)

    !real(ireals) :: diff2diff(solver%diffdof/2)

    real(irealLUT) :: coeff((solver%diffdof/2)**2), sumcoeff ! coefficients for each src=[1..8] and dst[1..8]
    integer(iintegers), pointer :: xinoutdof(:)
    real(ireals) :: dz, area_top, area_bot
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
      if(.not.allocated(irows)) allocate(irows(size(outgoing_offsets)))
      if(.not.allocated(icols)) allocate(icols(size(incoming_offsets)))
      if(.not.allocated(c    )) allocate(c    (size(outgoing_offsets)*size(incoming_offsets)))
      !print *,'icell', icell, 'faces_of_cell', faces_of_cell
      !print *,'icell', icell, 'incoming_offsets', incoming_offsets
      !print *,'icell', icell, 'outgoing_offsets', outgoing_offsets

      call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset2+iface), iintegers)) = iface
      enddo
      call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      !print *,'icell',icell,': foc',faces_of_cell
      !print *,'icell',icell,':',face_plex2bmc
      call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
        aspect_zx=wedgeorient(wedge_offset1+i5), &
        Cx=wedgeorient(wedge_offset4+i1), &
        Cy=wedgeorient(wedge_offset4+i2), &
        ldir=.False., coeff=coeff, ierr=ierr)

      if(ldebug_optprop) then
        do i = 1, size(incoming_offsets)
          associate(diff2diff => coeff(diff_plex2bmc(i):size(coeff):solver%diffdof/2))
            if(sum(diff2diff).gt.coeff_norm_err_tolerance) then
              print *,i,': bmcface', diff_plex2bmc(i), 'diff2diff gt one', diff2diff, &
                ':', sum(diff2diff), 'l1d', ierr.eq.OPP_1D_RETCODE, 'tol', coeff_norm_err_tolerance
              call CHKERR(1_mpiint, '1 energy conservation violated! '//ftoa(sum(diff2diff)))
            endif
          end associate
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

      if(ldebug) sumcoeff=0
      do j = 1, size(outgoing_offsets)
        do i = 1, size(incoming_offsets)
          if(ldebug) then
            if(outgoing_offsets(j).ge.0 .and. incoming_offsets(i).ge.0) &
              sumcoeff = sumcoeff + coeff((diff_plex2bmc(j)-1)*(solver%diffdof/2) + diff_plex2bmc(i))
          endif
          c((j-1)*size(outgoing_offsets) + i) = -coeff((diff_plex2bmc(j)-1)*(solver%diffdof/2) + diff_plex2bmc(i))
        enddo
      enddo
      call MatSetValuesLocal(A, &
        size(outgoing_offsets, kind=iintegers), outgoing_offsets, &
        size(incoming_offsets, kind=iintegers), incoming_offsets, c, INSERT_VALUES, ierr); call CHKERR(ierr)

      if(ldebug) then
        if(.not.approx(sumcoeff, sum(coeff), sqrt(epsilon(coeff)))) then
          do i = 1, solver%diffdof/2
          call CHKWARN(1_mpiint, 'diff2diff src '//itoa(i)//' : '//&
            ftoa(coeff(i:size(coeff):solver%diffdof/2)))
          enddo
          call CHKERR(1_mpiint, 'have we used all coefficients? '// &
            'used '//ftoa(sumcoeff)//' out of '//ftoa(sum(coeff)))
        endif
      endif
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

    if(myid.eq.0.and.ldebug) print *,'Creating Diffuse Matrix... done'
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

                call MatSetValuesLocal(A, i1, offset_Ein, i1, offset_Eout, -xalbedo(i1+offset_srfc), INSERT_VALUES, ierr)
                call CHKERR(ierr)
              endif
            enddo
          enddo
          call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
        endif


        sideward_bc_coeff = one
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-sideward_bc_coeff", &
          sideward_bc_coeff, lflg, ierr); call CHKERR(ierr)

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

                call MatSetValuesLocal(A, i1, offset_Ein, i1, offset_Eout, -sideward_bc_coeff, INSERT_VALUES, ierr)
                call CHKERR(ierr)
                if(ldebug.and.offset_Ein.eq.offset_Eout) call CHKERR(1_mpiint, &
                  'src and dst are the same :( ... should not happen here'// &
                  ' row '//itoa(offset_Ein)// &
                  ' col '//itoa(offset_Eout))
              elseif (numDof.eq.i0) then
              else
                call CHKERR(1_mpiint, 'dont expect to have this numbering of dof, please check')
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
  subroutine get_coeff(OPP, kabs, ksca, g, dz, aspect_zx, Cx, Cy, ldir, coeff, ierr, angles)
    class(t_optprop_wedge), intent(in) :: OPP
    real(ireals), intent(in)     :: kabs, ksca, g
    real(ireals),intent(in)      :: dz, aspect_zx
    real(ireals), intent(in)     :: Cx, Cy ! coordinates of upper triangle pts A,B,C in in (x,y)
    logical,intent(in)           :: ldir
    real(irealLUT),intent(out)     :: coeff(:)
    integer(mpiint), intent(out) :: ierr

    real(irealLUT),intent(in),optional  :: angles(2)

    real(irealLUT) :: dkabs, dksca, dg
    real(irealLUT) :: tauz, w0, param_phi

    logical :: lflg
    real(ireals) :: max_g_delta

    real(irealLUT), save :: max_g(4)
    integer(iintegers), save :: dimidx(10)
    logical, save :: linit_idx=.False.
    integer, parameter :: itaudir=1, itaudiff=6
    integer, parameter :: iw0dir=2 , iw0diff=7
    integer, parameter :: iphidir=3
    integer, parameter :: iCxdir=4 , iCxdiff=9
    integer, parameter :: iCydir=5 , iCydiff=10

    if(.not.linit_idx) then
      dimidx(itaudir) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'tau')
      dimidx(iw0dir ) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'w0')
      dimidx(iphidir) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'param_phi')
      dimidx(iCxdir ) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cx')
      dimidx(iCydir ) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cy')

      dimidx(itaudiff) = find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'tau')
      dimidx(iw0diff ) = find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'w0')
      dimidx(iCxdiff ) = find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'wedge_coord_Cx')
      dimidx(iCydiff ) = find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'wedge_coord_Cy')

      if(find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'g').ge.1) then
        max_g(1:2) = OPP%OPP_LUT%dirconfig%dims(&
          find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'g'))%vrange
      else
        max_g(1:2) = 0._irealLUT
      endif

      if(find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'g').ge.1) then
        max_g(3:4) = OPP%OPP_LUT%diffconfig%dims(&
          find_lut_dim_by_name(OPP%OPP_LUT%diffconfig, 'g'))%vrange
      else
        max_g(3:4) = 0._irealLUT
      endif

      if(ldebug) then
        print *,'found LUT max_g to be:', max_g
        print *,'found direct LUT range tau to be:', OPP%OPP_LUT%dirconfig%dims(dimidx(itaudir))%vrange
        print *,'found direct LUT range  w0 to be:', OPP%OPP_LUT%dirconfig%dims(dimidx(iw0dir ))%vrange
      endif

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
        "-plexrt_delta_scale_maxg", max_g_delta, lflg , ierr) ;call CHKERR(ierr)
      if(lflg) then
        max_g([2,4]) = max(max_g([1,3]), min(max_g([2,4]), real(max_g_delta, kind(max_g))))
        if(ldebug) print *,'found deltascaling max_g -> setting limits to:', cstr(ftoa(max_g), 'green')
      endif
      max_g([2,4]) = max(max_g([1,3]), max_g([2,4]) - sqrt(epsilon(max_g)))

      linit_idx = .True.
    endif

    dkabs = real(kabs, kind(dkabs))
    dksca = real(ksca, kind(dksca))
    dg    = real(g,    kind(dg   ))

    tauz = real((dkabs+dksca) * dz, irealLUT)
    w0 = real(dksca / max(dkabs+dksca, tiny(dksca)), irealLUT)

    if(present(angles)) then
      tauz = snap_limits(tauz, OPP%OPP_LUT%dirconfig%dims(dimidx(itaudir))%vrange)
      w0   = snap_limits(w0  , OPP%OPP_LUT%dirconfig%dims(dimidx(iw0dir ))%vrange)

      param_phi = max(OPP%OPP_LUT%dirconfig%dims(dimidx(iphidir))%vrange(1), &
        min(OPP%OPP_LUT%dirconfig%dims(dimidx(iphidir))%vrange(2), angles(1)))

      call delta_scale( dkabs, dksca, dg, max_g=max_g(2))

      call OPP%get_coeff(tauz, w0, real(dg, irealLUT), real(aspect_zx, irealLUT), &
        ldir, coeff, ierr, &
        angles=[param_phi, angles(2)], wedge_coords=real([Cx, Cy], irealLUT))

      !if(.not.ldir) then
      !  print *,'dir2diff optprop', &
      !    'kabs, ksca, g', kabs, ksca, g, new_line(''), &
      !    'dkabs, dksca, dg', dkabs, dksca, dg, new_line(''), &
      !    'tau w0, g', tauz, w0, dg, &
      !    '::', dksca / (dkabs+dksca), ':', max(dkabs+dksca, tiny(dksca))
      !endif
    else
      tauz = snap_limits(tauz, OPP%OPP_LUT%diffconfig%dims(dimidx(itaudiff))%vrange)
      w0   = snap_limits(w0  , OPP%OPP_LUT%diffconfig%dims(dimidx(iw0diff ))%vrange)

      call delta_scale( dkabs, dksca, dg, max_g=max_g(4))

      call OPP%get_coeff(tauz, w0, real(dg, irealLUT), real(aspect_zx, irealLUT), &
        ldir, coeff, ierr, &
        wedge_coords=real([Cx, Cy], irealLUT))
    endif

    if(ldebug) then
      if(present(angles)) then
        if(any(coeff.lt.0._irealLUT).or.any(coeff.gt.1._irealLUT+sqrt(epsilon(coeff)))) then
          if(     ldir) print *,'Lookup Coeffs for dir2dir', tauz, w0, aspect_zx, dg, angles,'::', coeff
          if(.not.ldir) print *,'Lookup Coeffs for dir2diff', tauz, w0, aspect_zx, dg, angles,'::', coeff
        endif
      else
        if(tauz.lt.1._irealLUT .and. all(coeff.le.0)) then
          print *,'Lookup Coeffs for diff2diff', tauz, w0, dg, aspect_zx, ':', Cx, Cy, '::', coeff
          call CHKERR(1_mpiint, 'Found all zero entries where I would not expect it')
        endif
      endif

      if(any(coeff.lt.0._irealLUT).or.any(coeff.gt.1._irealLUT+sqrt(epsilon(coeff)))) then
        print *,'Lookup Coeffs for', aspect_zx, tauz, w0, dg, angles,'::', coeff
        call CHKERR(1_mpiint, 'Found corrupted coefficients!')
      endif
    endif

    contains
      function snap_limits(val, vrange) result (res)
        real(irealLUT), intent(in) :: val, vrange(:)
        real(irealLUT) :: res
        res = max(vrange(1), min( val, vrange(2) ))
      end function
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

    if(lcalc_nca) then
      call plexrt_NCA_wrapper(solver, solution, ierr); call CHKERR(ierr)
    endif

    call update_absorption_norms_for_adaptive_spectral_integration()

    contains
      subroutine update_absorption_norms_for_adaptive_spectral_integration()
        real(ireals) :: norm3
        integer(mpiint) :: myid, comm, ierr
        if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
          call VecAXPY(abso_old, -one, solution%abso, ierr); call CHKERR(ierr) ! overwrite abso_old with difference to new one
          call VecNorm(abso_old, NORM_INFINITY, norm3, ierr); call CHKERR(ierr)

          call DMRestoreGlobalVector(solver%plex%abso_dm, abso_old, ierr)   ; call CHKERR(ierr)

          ! Save norm for later analysis
          solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
          solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

          solution%maxnorm( 1 ) = norm3
          solution%time( 1 )    = time

          if(ldebug) then
            call PetscObjectGetComm(solver%plex%abso_dm, comm, ierr); call CHKERR(ierr)
            call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
            if(myid.eq.0) &
              print *,'Updating error statistics for solution',solution%uid,'at time', &
                      solution%time(1),':: inf_norm',norm3,'[W] :: ', &
                      'hr_norm approx:',norm3*86.1,'[K/d]'
          endif
        endif !present(time) .and. solver%lenable_solutions_err_estimates
      end subroutine
  end subroutine

  subroutine compute_absorption(solver, solution)
    class(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution

    logical, parameter :: by_flx_divergence=.True.
    integer(mpiint) :: ierr

    if(.not.solution%lset) call CHKERR(1_mpiint, 'compute_absorption needs to get an initialized solution obj')

    if(.not.allocated(solution%abso)) then
      allocate(solution%abso)
      call DMCreateGlobalVector(solver%plex%abso_dm, solution%abso, ierr); call CHKERR(ierr)
    endif

    if(solution%lchanged) then
      if(solution%lsolar_rad .and. .not.allocated(solution%edir)) &
        call CHKERR(1_mpiint, 'solution%lsolar_rad true but edir not allocated')
      if(.not.allocated(solution%ediff)) call CHKERR(1_mpiint, 'ediff vec not allocated')

      ! Make sure we have the radiation vecs in plain energy units
      call scale_flx(solver, solver%plex, &
        solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
        solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
        solution, lWm2=.False., logevent=solver%logs%scale_flx)

      call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

      if(solution%lsolar_rad) then
        call PetscLogEventBegin(solver%logs%compute_absorption, ierr)
        if(by_flx_divergence) then
          call compute_edir_absorption(solver%plex, solution%edir, solution%abso)
        else
          call compute_edir_absorption_by_coeff_divergence(solver, solution%edir, solution%abso)
        endif
        call PetscLogEventEnd(solver%logs%compute_absorption, ierr)

        call PetscLogEventBegin(solver%logs%debug_output, ierr)
        call PetscObjectSetName(solution%abso, 'abso_direct_'//itoa(solution%uid), ierr);call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_abso_direct', ierr); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%debug_output, ierr)
      endif

      call PetscLogEventBegin(solver%logs%compute_absorption, ierr)
        if(by_flx_divergence) then
          call compute_ediff_absorption(solver%plex, solver%IS_diff_in_out_dof, solution%ediff, solution%abso)
        else
          call compute_ediff_absorption_by_coeff_divergence(solver, solver%plex, solution%ediff, solution%abso)
        endif
      call PetscLogEventEnd(solver%logs%compute_absorption, ierr)
    endif

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

    real(ireals) :: cell_abso, inv_volume

    logical :: lsrc ! is src or destination of solar beam

    if(.not.allocated(plex)) call CHKERR(1_mpiint, 'called compute_edir_absorption but plex is not allocated')
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%edir_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_edir_absorption with an unallocated abso vec')

    !if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%wedge_orientation_dm)) &
      call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated: plex%wedge_orientation_dm?')
    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i2, wedge_offset, ierr); call CHKERR(ierr)
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      cell_abso = 0
      do i = 1, size(faces_of_cell)
        iface = faces_of_cell(i)
        lsrc = wedgeorient(wedge_offset+i) .gt. zero

        call PetscSectionGetDof(edir_section, iface, numDof, ierr); call CHKERR(ierr)
        call PetscSectionGetOffset(edir_section, iface, edir_offset, ierr); call CHKERR(ierr)
        do idof = 1, numDof
          if(lsrc) then ! sun is shining into this face
            cell_abso = cell_abso + xedir(edir_offset+idof)
            !print *,i1+icell,'Dir Abso',abso_offset,'=',cell_abso,'( did +',xedir(edir_offset+i1),')'
          else
            cell_abso = cell_abso - xedir(edir_offset+idof)
            !print *,i1+icell,'Dir Abso',abso_offset,'=',cell_abso,'( did -',xedir(edir_offset+i1),')'
          endif
        enddo
      enddo
      xabso(i1+abso_offset) = xabso(i1+abso_offset) + cell_abso * inv_volume
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end subroutine
  subroutine compute_edir_absorption_by_coeff_divergence(solver, edir, abso)
    class(t_plex_solver), intent(in) :: solver
    type(tVec), intent(in) :: edir
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_edir

    real(ireals), pointer :: xedir(:), xabso(:)
    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

    type(tPetscSection) :: abso_section, edir_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface
    integer(iintegers),pointer :: faces_of_cell(:)
    integer(mpiint) :: myid, ierr

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:), wedgeorient(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset
    integer(iintegers) :: wedge_offset1, wedge_offset2, wedge_offset3, wedge_offset4

    real(ireals) :: cell_abso, inv_volume
    real(ireals) :: aspect_zx, dz, param_theta, param_phi
    integer(iintegers) :: bmcsrcdof, dof_src_offset, icol, in_dof, isrc, isrc_side, numSrc, bmcsrcside
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: plex2bmc(solver%dirdof)

    logical :: lsrc(solver%dirdof) ! is src or destination of solar beam

    real(irealLUT) :: coeff_dir2diff(solver%dirdof*(solver%diffdof/2))
    real(irealLUT) :: coeff_dir2dir(solver%dirdof**2)
    !real(ireals) :: dir2dir (solver%dirdof)
    !real(ireals) :: dir2diff(solver%diffdof/2)
    real(ireals) :: div

    associate(plex => solver%plex)
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%edir_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_edir_absorption with an unallocated abso vec')

    !if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%g   , xg   , ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%wedge_orientation_dm)) &
      call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated: plex%wedge_orientation_dm?')
    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1
      call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i2, wedge_offset3, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      param_phi   = wedgeorient(wedge_offset1+i3)
      param_theta = wedgeorient(wedge_offset1+i4)
      aspect_zx   = wedgeorient(wedge_offset1+i5)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset2+iface), iintegers)) = iface
        lsrc(iface) = wedgeorient(wedge_offset3+iface) .gt. zero
      enddo
      call face_idx_to_bmc_idx(solver%dirtop, solver%dirside, face_plex2bmc, plex2bmc)

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      call get_coeff(solver%OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
        aspect_zx=wedgeorient(wedge_offset1+i5), &
        Cx=wedgeorient(wedge_offset4+i1), &
        Cy=wedgeorient(wedge_offset4+i2), &
        ldir=.True., coeff=coeff_dir2dir, ierr=ierr, &
        angles=[real(param_phi, irealLUT), real(param_theta, irealLUT)])

      call get_coeff(solver%OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
        aspect_zx=wedgeorient(wedge_offset1+i5), &
        Cx=wedgeorient(wedge_offset4+i1), &
        Cy=wedgeorient(wedge_offset4+i2), &
        ldir=.False., coeff=coeff_dir2diff, ierr=ierr, &
        angles=[real(param_phi, irealLUT), real(param_theta, irealLUT)])

      cell_abso = 0

      in_dof = 0 ! this counts plex dofs from 1 to solver%dirdof
      do isrc_side = 1, size(faces_of_cell)
        iface = faces_of_cell(isrc_side)

        call PetscSectionGetDof(edir_section, iface, numSrc, ierr); call CHKERR(ierr)
        do isrc = 1, numSrc
          in_dof = in_dof+1
          bmcsrcdof = plex2bmc(in_dof)

          call get_side_and_offset_from_total_bmc_dofs(solver%dirtop, solver%dirside, bmcsrcdof, &
            bmcsrcside, dof_src_offset)

          associate( &
            dir2dir  => coeff_dir2dir (bmcsrcdof:size(coeff_dir2dir ):solver%dirdof), &
            dir2diff => coeff_dir2diff(bmcsrcdof:size(coeff_dir2diff):solver%dirdof)  &
          )
          div = one-sum(dir2dir)-sum(dir2diff)
          div = max(zero, min(one, div))
          end associate

          call PetscSectionGetOffset(edir_section, iface, icol, ierr); call CHKERR(ierr)
          icol = icol + dof_src_offset

          if(lsrc(isrc_side)) then
            cell_abso = cell_abso + div * xedir(i1+icol)
            !print *,'Adding dir cell abso cell='//itoa(icell)//' src-face='//itoa(iface)//' bmcsrc='//itoa(bmcsrcdof)// &
            !  ' edir='//ftoa(xedir(i1+icol)), 'div', div, 'lsrc', lsrc(isrc_side)
          endif
        enddo
      enddo

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      !print *,icell, cstr('icell abso '//ftoa(xabso(i1+abso_offset))//' += '//ftoa(cell_abso * inv_volume), 'red')
      xabso(i1+abso_offset) = xabso(i1+abso_offset) + cell_abso * inv_volume
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    call VecRestoreArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%g   , xg   , ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end associate
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
    integer(iintegers) :: icell, i
    integer(mpiint) :: myid, ierr

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    integer(iintegers), pointer :: faces_of_cell(:)
    real(ireals) :: area_top

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset

    integer(iintegers), pointer :: xinoutdof(:)

    real(ireals) :: inv_volume, cell_abso

    if(.not.allocated(plex)) stop 'called compute_ediff_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%ediff_dm)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with an unallocated abso vec')

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
      cell_abso = zero
      do i = 1, size(incoming_offsets)
        if(incoming_offsets(i).ge.i0) cell_abso = cell_abso + xediff(i1+incoming_offsets(i))
        if(outgoing_offsets(i).ge.i0) cell_abso = cell_abso - xediff(i1+outgoing_offsets(i))
        if(.False.) then
          call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
          call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
          area_top = geoms(i1+geom_offset)
          call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

          print *,'icell', i1+icell, incoming_offsets(i), outgoing_offsets(i), 'abso', &
            xabso(i1+abso_offset), &
            '+', xediff(i1+incoming_offsets(i))/area_top, &
            '-', xediff(i1+outgoing_offsets(i))/area_top, &
            '=', xabso(i1+abso_offset) + cell_abso * inv_volume
        endif
      enddo

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      xabso(i1+abso_offset) = xabso(i1+abso_offset) + cell_abso * inv_volume
    enddo
    call ISRestoreIndicesF90(IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)
    !if(ldebug) print *,'plex_rt::compute_ediff_absorption....finished'
  end subroutine
  subroutine compute_ediff_absorption_by_coeff_divergence(solver, plex, ediff, abso)
    class(t_plex_solver), intent(in) :: solver
    type(t_plexgrid), allocatable, intent(in) :: plex
    type(tVec), intent(in) :: ediff
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_ediff

    real(ireals), pointer :: xediff(:), xabso(:)
    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers), pointer :: xinoutdof(:)

    type(tPetscSection) :: abso_section, ediff_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, i, iface
    integer(mpiint) :: myid, ierr

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    real(ireals) :: dz, div

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates and orientation vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: geom_offset, abso_offset
    integer(iintegers) :: wedge_offset1, wedge_offset2, wedge_offset4

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: diff_plex2bmc(solver%diffdof/2)

    real(irealLUT) :: coeff((solver%diffdof/2)**2) ! coefficients for each src=[1..8] and dst[1..8]

    real(ireals) :: inv_volume, cell_abso

    if(.not.allocated(plex)) stop 'called compute_ediff_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%ediff_dm)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) &
      call CHKERR(myid+1, 'called compute_ediff_absorption with an unallocated abso vec')

    !if(ldebug) print *,'plex_rt::compute_ediff_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%g   , xg   , ierr); call CHKERR(ierr)

    call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call ISGetIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)
    do icell = cStart, cEnd-1
      call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      call get_in_out_dof_offsets(solver%IS_diff_in_out_dof, icell, incoming_offsets, outgoing_offsets, xinoutdof)

      call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)


      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset2+iface), iintegers)) = iface
      enddo
      call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)

      call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
      call get_coeff(solver%OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), dz, &
        aspect_zx=wedgeorient(wedge_offset1+i5), &
        Cx=wedgeorient(wedge_offset4+i1), &
        Cy=wedgeorient(wedge_offset4+i2), &
        ldir=.False., coeff=coeff, ierr=ierr)

      cell_abso = zero

      do i = 1, size(incoming_offsets)
        if(incoming_offsets(i).ge.0) then
          associate( diff2diff => coeff(diff_plex2bmc(i):size(coeff):solver%diffdof/2) )
            div = one - sum(diff2diff)
            div = max(zero, min(one, one - sum(diff2diff)))
            cell_abso = cell_abso + div * xediff(i1+incoming_offsets(i))
            !print *,'Adding diffuse cell abso cell='//itoa(icell)//' inc_offset='//itoa(incoming_offsets(i))// &
            !  ' ediff='//ftoa(xediff(i1+incoming_offsets(i))), 'div', div
          end associate
        endif
      enddo

      call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr)
      inv_volume = one / geoms(i1+geom_offset)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      xabso(i1+abso_offset) = xabso(i1+abso_offset) + cell_abso * inv_volume
    enddo
    call ISRestoreIndicesF90(solver%IS_diff_in_out_dof, xinoutdof, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call DMRestoreLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)

    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%g   , xg   , ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    !if(ldebug) print *,'plex_rt::compute_ediff_absorption....finished'
  end subroutine

    !> @brief renormalize fluxes with the size of a face(sides or lid)
    subroutine scale_flx(solver, plex, &
        dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, &
        diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, &
        solution, lWm2, logevent)
      class(t_plex_solver), intent(inout) :: solver
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
            call scale_facevec(plex, plex%edir_dm, solver%dirtop, solver%dirside, &
              dir_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
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
          call scale_facevec(plex, plex%ediff_dm, solver%difftop, solver%diffside, &
            diff_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
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
      integer(iintegers) :: Ncol, ke1, uid, i, k, iface, icell, voff, idof
      integer(iintegers), pointer :: xtoa_faces(:), cell_support(:)

      type(tPetscSection) :: abso_section, ediff_section, edir_section
      real(ireals), pointer :: xediff(:), xedir(:), xabso(:)

      integer(mpiint) :: myid, ierr

      call PetscLogEventBegin(solver%logs%get_result, ierr)

      uid = get_arg(0_iintegers, opt_solution_uid)
      if(.not.solver%solutions(uid)%lset) &
        call CHKERR(1_mpiint, 'You tried to retrieve results from a solution uid which has not yet been calculated')

      ke1 = solver%plex%Nlay+1

      call DMGetStratumIS(solver%plex%ediff_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      if (toa_ids.ne.PETSC_NULL_IS) then
        call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)
      else
        Ncol=0
      endif

      call check_size_or_allocate(redn , [ke1, Ncol])
      call check_size_or_allocate(reup , [ke1, Ncol])
      call check_size_or_allocate(rabso, [ke1-1, Ncol])
      if(present(redir)) call check_size_or_allocate(redir, [ke1, Ncol])

      associate(solution => solver%solutions(uid))

        if(.not.allocated(solution%ediff)) call CHKERR(1_mpiint, 'ediff vec not allocated')
        if(.not.allocated(solution%abso)) call CHKERR(1_mpiint, 'abso vec not allocated')
        if(solution%lchanged) call CHKERR(1_mpiint, 'tried to get results from unrestored solution -- call restore_solution first')

        if(present(redir) .and. .not.solution%lsolar_rad) then
          call CHKWARN(1_mpiint, 'you asked for direct radiation but solution does not have it')
          redir(:,:) = 0
        endif

        call scale_flx(solver, solver%plex, &
          solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, &
          solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, &
          solution, lWm2=.True., logevent=solver%logs%scale_flx)

        if(present(redir).and.solution%lsolar_rad) then
          call DMGetSection(solver%plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
        endif
        call DMGetSection(solver%plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
        call DMGetSection(solver%plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

        if(present(redir).and.solution%lsolar_rad) then
          call VecGetArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
        endif
        call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
        call VecGetArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

        if (toa_ids.ne.PETSC_NULL_IS) then
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

            if(present(redir).and.solution%lsolar_rad) then
              do k = 0, ke1-1
                call PetscSectionGetFieldOffset(edir_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
                redir(i1+k,i) = zero
                do idof = 1, solver%dirtop%dof
                  redir(i1+k,i) = redir(i1+k,i) + xedir(voff+idof)
                enddo
                redir(i1+k,i) = redir(i1+k,i) / real(solver%dirtop%dof, ireals)
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
          call ISRestoreIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
        endif

        if(present(redir).and.solution%lsolar_rad) then
          call VecRestoreArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
        endif
        call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

        if(ldebug) then
          call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)
          if(myid.eq.0) then
            if(present(redir).and.solution%lsolar_rad) then
              print *,'Get Result, k     Edir                   Edn                      Eup                  abso'
              do k = 1, ke1-1
                print *,k, meanval(redir(k,:)), meanval(redn(k,:)), meanval(reup(k,:)), meanval(rabso(k,:))!, &
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

  subroutine vacuum_domain_boundary(dm, kabs, ksca, g, fill_kabs, fill_ksca, fill_g)
    type(tDM), intent(in) :: dm
    type(tVec), intent(inout) :: kabs, ksca, g
    real(ireals), intent(in) :: fill_kabs, fill_ksca, fill_g

    type(tIS) :: IS_side_faces
    integer(iintegers), pointer :: iside_faces(:), cell_support(:)
    integer(iintegers) :: o, iface_side

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
    integer(mpiint) :: comm, myid, ierr

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    if(myid.eq.0) then
      print *,'Vacuum shell of domain with ', fill_kabs, fill_ksca, fill_g
    endif

    call DMGetStratumIS(dm, 'DomainBoundary', SIDEFACE, IS_side_faces, ierr); call CHKERR(ierr)
    if (IS_side_faces.eq.PETSC_NULL_IS) then ! dont have side boundary faces
    else
      call VecGetArrayF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayF90(g   , xg   , ierr); call CHKERR(ierr)

      call ISGetIndicesF90(IS_side_faces, iside_faces, ierr); call CHKERR(ierr)
      do o = 1, size(iside_faces)
        iface_side = iside_faces(o)

        call DMPlexGetSupport(dm, iface_side, cell_support, ierr); call CHKERR(ierr)
        xkabs(i1+cell_support(1)) = fill_kabs
        xksca(i1+cell_support(1)) = fill_ksca
        xg   (i1+cell_support(1)) = fill_g
        call DMPlexRestoreSupport(dm, iface_side, cell_support, ierr); call CHKERR(ierr)
      enddo
      call VecRestoreArrayF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(g   , xg   , ierr); call CHKERR(ierr)
    endif
  end subroutine
end module
