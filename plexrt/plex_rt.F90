module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options

  use m_helper_functions, only: CHKERR, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, strF2C, get_arg, &
    vec_proj_on_plane, cross_3d, norm, rotation_matrix_world_to_local_basis, &
    approx, swap, delta_scale, itoa

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, i6, i7, i8, &
    zero, one, pi, EXP_MINVAL, EXP_MAXVAL

  use m_plex_grid, only: t_plexgrid, compute_face_geometry, &
    setup_edir_dmplex, setup_ediff_dmplex, setup_abso_dmplex, &
    orient_face_normals_along_sundir, compute_wedge_orientation, is_solar_src, get_inward_face_normal, &
    facevec2cellvec, icell_icon_2_plex, iface_top_icon_2_plex, get_vertical_cell_idx, &
    get_top_bot_face_of_cell, destroy_plexgrid, determine_diff_incoming_outgoing_offsets, &
    TOAFACE, BOTFACE, SIDEFACE

  use m_optprop, only : t_optprop, t_optprop_wedge_5_8
  use m_optprop_parameters, only : OPP_LUT_ALL_ANGLES

  use m_pprts_base, only : t_state_container, prepare_solution, destroy_solution

  implicit none

  private

  public :: t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, &
    compute_face_geometry, set_plex_rt_optprop, destroy_plexrt_solver, &
    plexrt_get_result

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

    logical :: lenable_solutions_err_estimates=.True.
  end type

  logical, parameter :: ldebug=.False.
  contains
    subroutine init_plex_rt_solver(plex, solver)
      type(t_plexgrid), intent(in) :: plex
      type(t_plex_solver), allocatable, intent(inout) :: solver

      call read_commandline_options(plex%comm)

      if(allocated(solver)) call CHKERR(1_mpiint, 'Should not call init_plex_rt_solver with already allocated solver object')
      allocate(solver)

      allocate(solver%plex)
      solver%plex = plex

      allocate(t_optprop_wedge_5_8::solver%OPP)
      call solver%OPP%init(plex%comm)
    end subroutine

    subroutine destroy_plexrt_solver(solver, lfinalizepetsc)
      type(t_plex_solver), intent(inout) :: solver
      logical, intent(in) :: lfinalizepetsc

      integer(iintegers) :: uid
      integer(mpiint) :: ierr

      if(allocated(solver%plex)) then
        call destroy_plexgrid(solver%plex)
      endif

      if(allocated(solver%OPP)) then
        call solver%OPP%destroy()
        deallocate(solver%OPP)
      endif

      if(allocated(solver%kabs)) deallocate(solver%kabs)
      if(allocated(solver%ksca)) deallocate(solver%ksca)
      if(allocated(solver%g   )) deallocate(solver%g   )
      if(allocated(solver%plck)) deallocate(solver%plck)

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

      if(lfinalizepetsc) then
        call PetscFinalize(ierr) ;call CHKERR(ierr)
      endif
    end subroutine

    subroutine set_plex_rt_optprop(solver, vlwc, viwc)
      use m_helper_functions, only : delta_scale
      type(t_plex_solver), allocatable, intent(inout) :: solver
      type(tVec),intent(in), optional :: vlwc, viwc
      real(ireals), pointer :: xlwc(:), xiwc(:)
      real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

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

      call VecGetArrayF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solver%ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solver%g   , xg   , ierr); call CHKERR(ierr)
      do i = cStart, cEnd-1
        kabs_tot = rayleigh*(one-w0)
        ksca_tot = rayleigh*w0
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

    subroutine run_plex_rt_solver(solver, lthermal, lsolar, sundir, opt_solution_uid, opt_solution_time)
      type(t_plex_solver), allocatable, intent(inout) :: solver
      logical, intent(in) :: lthermal, lsolar
      real(ireals), intent(in) :: sundir(3) ! cartesian direction of sun rays, norm of vector is the energy in W/m2
      integer(iintegers), intent(in), optional :: opt_solution_uid
      real(ireals), intent(in), optional :: opt_solution_time

      integer(iintegers) :: suid

      integer(mpiint) :: myid, ierr

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

      if(.not.allocated(solver%plex%geom_dm))  call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      if(.not.allocated(solver%plex%edir_dm))  call setup_edir_dmplex(solver%plex, solver%plex%edir_dm)
      if(.not.allocated(solver%plex%ediff_dm)) call setup_ediff_dmplex(solver%plex, solver%plex%ediff_dm)
      if(.not.allocated(solver%plex%abso_dm))  call setup_abso_dmplex(solver%plex, solver%plex%abso_dm)

      ! Prepare the space for the solution
      suid = get_arg(i0, opt_solution_uid)

      if(.not.solver%solutions(suid)%lset) then
        call prepare_solution(solver%plex%edir_dm, solver%plex%ediff_dm, solver%plex%abso_dm, &
          lsolar=lsolar, solution=solver%solutions(suid))
      endif


      ! Wedge Orientation is used in solar and thermal case alike
      call orient_face_normals_along_sundir(solver%plex, sundir)
      call compute_wedge_orientation(solver%plex, sundir, solver%plex%wedge_orientation_dm, &
                                     solver%plex%wedge_orientation)

      associate( solution => solver%solutions(suid) )

      if(solution%lsolar_rad) then
        ! Output of wedge_orient vec
        call create_edir_src_vec(solver%plex, solver%plex%edir_dm, norm(sundir), &
                                 solver%kabs, solver%ksca, &
                                 sundir/norm(sundir), solver%dirsrc)

        ! Output of srcVec
        if(ldebug) then
          call scale_facevec(solver%plex, solver%plex%edir_dm, solver%dirsrc, lW_to_Wm2=.True.)
          call facevec2cellvec(solver%plex, solver%plex%edir_dm, solver%dirsrc)
          call scale_facevec(solver%plex, solver%plex%edir_dm, solver%dirsrc, lW_to_Wm2=.False.)
        endif

        ! Create Direct Matrix
        call create_edir_mat(solver%plex, solver%OPP, solver%kabs, solver%ksca, solver%g, solver%Mdir)

        ! Solve Direct Matrix
        call solve_plex_rt(solver%plex%edir_dm, solver%dirsrc, solver%Mdir, solver%kspdir, solution%edir, 'dir_')
        call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, &
                                        '-show_edir_vec_global', ierr); call CHKERR(ierr)
        solution%lWm2_dir = .False.
        solution%lchanged = .True.

        ! Output of Edir
        if(ldebug) then
          call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.True.)
          call facevec2cellvec(solver%plex, solver%plex%edir_dm, solution%edir)
          call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.False.)
        endif
      endif

      call create_ediff_src_vec(solver%plex, solver%OPP, solver%plex%ediff_dm, &
        solver%kabs, solver%ksca, solver%g, solver%plck, solver%albedo, solver%srfc_emission, &
        solver%diffsrc, solver%plex%edir_dm, solution%edir)

      ! Output of Diffuse Src Vec
      if(ldebug) then
        call scale_facevec(solver%plex, solver%plex%ediff_dm, solver%diffsrc, lW_to_Wm2=.True.)
        call facevec2cellvec(solver%plex, solver%plex%ediff_dm, solver%diffsrc)
        call scale_facevec(solver%plex, solver%plex%ediff_dm, solver%diffsrc, lW_to_Wm2=.False.)
      endif

      call create_ediff_mat(solver%plex, solver%OPP, solver%kabs, solver%ksca, solver%g, solver%albedo, solver%Mdiff)
      call solve_plex_rt(solver%plex%ediff_dm, solver%diffsrc, solver%Mdiff, solver%kspdiff, solution%ediff, 'diff_')
      call PetscObjectSetName(solution%ediff, 'ediff', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, &
        '-show_ediff_vec_global', ierr); call CHKERR(ierr)
      solution%lWm2_dir = .False.
      solution%lchanged = .True.

      ! Output of Ediff
      if(ldebug) then
        call scale_facevec(solver%plex, solver%plex%ediff_dm, solution%ediff, lW_to_Wm2=.True.)
        call facevec2cellvec(solver%plex, solver%plex%ediff_dm, solution%ediff)
        call scale_facevec(solver%plex, solver%plex%ediff_dm, solution%ediff, lW_to_Wm2=.False.)
      endif

      ! Bring solution into a coherent state, i.e. update absorption etc.
      call restore_solution(solver, solution, opt_solution_time)
      end associate
    end subroutine

    subroutine create_edir_src_vec(plex, edirdm, E0, kabs, ksca, sundir, srcVec)
      type(t_plexgrid), allocatable, intent(in) :: plex
      type(tDM),allocatable, intent(in) :: edirdm
      real(ireals), intent(in) :: E0, sundir(3)
      type(tVec), allocatable, intent(in) :: kabs, ksca
      type(tVec), allocatable, intent(inout) :: srcVec

      integer(iintegers) :: i, voff
      type(tPetscSection) :: s

      type(tVec) :: localVec, lambertVec
      real(ireals), pointer :: xv(:), xlambert(:)
      real(ireals), pointer :: xkabs(:), xksca(:)

      type(tIS) :: boundary_ids
      integer(iintegers) :: iface, icell, k, labelval

      integer(iintegers), pointer :: xx_v(:)
      integer(iintegers), pointer :: cell_support(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: area, mu, face_normal(3)

      integer(mpiint) :: myid, ierr

      if(.not.allocated(plex)) stop 'called create_src_vec but plex is not allocated'
      call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

      if(ldebug .and. myid.eq.0) print *,myid,'plex_rt::create_src_vec....'

      if(.not.allocated(edirdm)) call CHKERR(myid+1, 'called create_src_vec but edirdm is not allocated')

      if(.not.allocated(plex%geom_dm)) call CHKERR(myid+1, 'get_normal_of_first_TOA_face::needs allocated geom_dm first')

      call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetSection(edirdm, s, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(s, PETSC_NULL_SECTION, '-show_src_section', ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)

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
        call PetscObjectViewFromOptions(boundary_ids, PETSC_NULL_IS, '-show_IS_Boundary_TOA', ierr); call CHKERR(ierr)

        call ISGetIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)

        ! First set the TOA boundary fluxes on faces
        do i = 1, size(xx_v)
          iface = xx_v(i)
          call DMPlexGetSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

          call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
          area = geoms(geom_offset+i7)

          call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)

          if(is_solar_src(face_normal, sundir)) then
            mu = dot_product(sundir, face_normal)
            call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
            xv(voff+1) = E0 * area * mu
          endif
        enddo
        call ISRestoreIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)
      endif ! TOA boundary ids

      ! Then define sideways boundary conditions
      call DMGetLocalVector(edirdm, lambertVec,ierr); call CHKERR(ierr)
      call PetscObjectSetName(lambertVec, 'srcVec_side_lambert', ierr);call CHKERR(ierr)
      call VecSet(lambertVec, zero, ierr); call CHKERR(ierr)

      do iface = plex%fStart, plex%fEnd-1
        k = plex%zindex(iface)
        call DMLabelGetValue(plex%domainboundarylabel, iface, labelval, ierr); call CHKERR(ierr)
        if(labelval.eq.SIDEFACE .and. k.eq.1) then ! Side face at TOA
          ! Then check if the cell below is a sideward boundary cell
          call DMPlexGetSupport(edirdm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          call compute_lambert_beer(plex, icell, E0, sundir, xkabs, xksca, lambertVec)
        endif
      enddo

      call VecGetArrayReadF90(lambertVec, xlambert, ierr); call CHKERR(ierr)
      do iface = plex%fStart, plex%fEnd-1
        call DMLabelGetValue(plex%domainboundarylabel, iface, labelval, ierr); call CHKERR(ierr)
        if(labelval.eq.SIDEFACE)then
          call DMPlexGetSupport(edirdm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          !print *,'Face', iface, 'is SIDEFACE cell:',icell

          call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)

          if(is_solar_src(face_normal, sundir)) then
            call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
            !print *,'Face', iface, 'and Solar SOURCE', xlambert(voff+i1)
            xv(i1+voff) = xlambert(voff+i1)
          endif
        endif
      enddo
      call VecRestoreArrayReadF90(lambertVec, xlambert, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      call DMLocalToGlobalBegin(edirdm, localVec, INSERT_VALUES, srcVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (edirdm, localVec, INSERT_VALUES, srcVec, ierr); call CHKERR(ierr)
      call PetscObjectSetName(srcVec, 'srcVec', ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(srcVec, PETSC_NULL_VEC, '-show_src_vec_global', ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(edirdm, localVec, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)

      if(ldebug .and. myid.eq.0) print *,myid,'plex_rt::create_src_vec....finished'
      if(ldebug) then
        call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
      endif
    end subroutine

    !> @brief setup source term for diffuse radiation
    !> @details this is either direct radiation scattered into one of the diffuse coeffs:
    !> \n direct source term is
    !> \n   direct radiation times the dir2diff coeffs
    !> \n or it may be that we have a source term due to thermal emission --
    !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
    !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
    subroutine create_ediff_src_vec(plex, OPP, ediffdm, kabs, ksca, g, plckVec, &
      albedo, srfc_emission, srcVec, edirdm, edirVec)
      type(t_plexgrid), intent(in) :: plex
      class(t_optprop), intent(in) :: OPP
      type(tDM), allocatable, intent(in) :: ediffdm
      type(tVec), allocatable, intent(in) :: kabs, ksca, g, plckVec ! cell1_dm
      type(tVec), allocatable, intent(in) :: albedo, srfc_emission  ! srfc_boundary_dm
      type(tVec), allocatable, intent(inout) :: srcVec

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

      call set_solar_source()
      call set_thermal_source()

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
          integer(iintegers) :: i, icell, iface, isrc, idst

          integer(iintegers) :: face_plex2bmc(5), wedge_offset, diff_plex2bmc(8)
          real(ireals), pointer :: xedir(:)
          real(ireals) :: zenith, azimuth
          real(ireals) :: dir2diff(8)
          logical :: lsrc(5)

          integer(iintegers) :: zindex
          real(ireals) :: dz, coeff(5*8)

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

          call VecGetArrayReadF90(ledirVec, xedir, ierr); call CHKERR(ierr)

          do icell = plex%cStart, plex%cEnd-1
            call DMPlexGetCone(ediffdm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

            call determine_diff_incoming_outgoing_offsets(plex, ediffdm, icell, incoming_offsets, outgoing_offsets)
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

            zindex = plex%zindex(icell)
            dz = plex%hhl(zindex) - plex%hhl(zindex+1)

            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff, &
              angles=[rad2deg(azimuth), rad2deg(zenith)])

            do i = 1, size(faces_of_cell)
              iface = faces_of_cell(i)
              !print *,'icell '//itoa(icell)//' i '//itoa(i)//' face ',iface, ':plex2bmc', face_plex2bmc

              if(lsrc(i)) then
                dir2diff = coeff(face_plex2bmc(i):size(coeff):i5) ! dir2diff in src ordering
                call PetscSectionGetOffset(edirSection, iface, isrc, ierr); call CHKERR(ierr)

                do idst = 1, size(dir2diff)
                    !print *,'Setting src for dst', outgoing_offsets(idst),'= src, edir', isrc, xedir(i1+isrc),&
                    !  '* idst, p2bmc coeff', idst, diff_plex2bmc(idst), dir2diff(diff_plex2bmc(idst))
                    xb(i1+outgoing_offsets(idst)) = xb(i1+outgoing_offsets(idst)) + &
                      xedir(i1+isrc) * dir2diff(diff_plex2bmc(idst))
                enddo
              endif
            enddo ! enddo iface

            call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
          enddo

          call set_Edir_srfc_reflection(edirSection, xedir)

          call VecRestoreArrayReadF90(ledirVec, xedir, ierr); call CHKERR(ierr)
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

              !print *,'Srfc Reflection of Edir',offset_Eup, offset_Edir, offset_srfc, &
              !  '::', xedir(i1+offset_Edir), '*', xalbedo(i1+offset_srfc)

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

          integer(iintegers) :: idof, iface, num_dof
          real(ireals) :: area

          if(.not.allocated(plckVec)) return

          call DMGetLocalVector(ediffdm, thermal_src_vec, ierr); call CHKERR(ierr)
          call VecSet(thermal_src_vec, zero, ierr); call CHKERR(ierr)
          call VecGetArrayF90(thermal_src_vec, xsrc, ierr); call CHKERR(ierr)

          call set_cell_emissions(xsrc)
          call thermal_srfc_emission(xsrc)
          call thermal_horizontal_boundary_mirror(xsrc)

          ! Scaling from [W/m2] to Energy [W]
          do iface = plex%fStart, plex%fEnd-1
            call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
            area = geoms(geom_offset+i7)

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
          integer(iintegers) :: wedge_offset

          integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)
          integer(iintegers) :: i, zindex, icell, icol
          real(ireals) :: dz, coeff(8**2) ! coefficients for each src=[1..8] and dst[1..8]

          real(ireals), pointer :: xplanck(:)
          real(ireals) :: diff2diff(8), emissivity

          call VecGetArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)

          do icell = plex%cStart, plex%cEnd-1

            call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
            call determine_diff_incoming_outgoing_offsets(plex, plex%ediff_dm, icell, incoming_offsets, outgoing_offsets)

            call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)

            do i = 1, size(faces_of_cell)
              face_plex2bmc(int(wedgeorient(wedge_offset+i3+i), iintegers)) = i
            enddo
            call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

            zindex = plex%zindex(icell)
            dz = plex%hhl(zindex) - plex%hhl(zindex+1)

            call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff)

            do i = 1, size(outgoing_offsets)
              icol = outgoing_offsets(i)

              ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering
              diff2diff = coeff(diff_plex2bmc(i):size(coeff):i8)

              emissivity = one - sum(diff2diff)

              xsrc(i1+icol) = xplanck(i1+icell) * pi * emissivity
            enddo
          enddo ! icell
          call VecRestoreArrayReadF90(plckVec, xplanck, ierr); call CHKERR(ierr)
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

        subroutine thermal_horizontal_boundary_mirror(xsrc)
          real(ireals), pointer :: xsrc(:)
          type(tIS) :: bc_ids
          integer(iintegers), pointer :: xi(:)
          integer(iintegers) :: i, idof, iface, offset_Ein, offset_Eout, numDof, offset_srfc

          call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', SIDEFACE, bc_ids, ierr); call CHKERR(ierr)
          if (bc_ids.eq.PETSC_NULL_IS) then ! dont have boundary points
          else
            call ISGetIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
            do i = 1, size(xi)
              iface = xi(i)
              ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
              ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
              call PetscSectionGetFieldOffset(ediffSection, iface, i0, offset_Ein, ierr); call CHKERR(ierr)
              call PetscSectionGetFieldOffset(ediffSection, iface, i1, offset_Eout, ierr); call CHKERR(ierr)
              call PetscSectionGetFieldDof(ediffSection, iface, i0, numDof, ierr); call CHKERR(ierr)

              call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

              do idof = 1, numDof
                xsrc(offset_Ein+idof) = xsrc(offset_Eout+idof) ! mirror the outgoing emission
              enddo
            enddo
            call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          endif

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

    subroutine compute_lambert_beer(plex, itopcell, E0, sundir, xkabs, xksca, lambertVec)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: itopcell
      real(ireals), intent(in) :: E0, sundir(3) ! E0 in W/m2
      real(ireals), pointer, intent(in) :: xkabs(:), xksca(:)
      type(tVec),intent(inout) :: lambertVec

      type(tPetscSection) :: s
      real(ireals), pointer :: xv(:) ! vertical entries are the direct radiation flux through the tilted plane
      integer(iintegers), pointer :: faces_of_cell(:)
      integer(iintegers) :: i, icell, iface_side, iface_top, iface_bot, k, voff
      integer(iintegers), allocatable :: cell_idx(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: face_normal(3), face_normal_top(3)

      real(ireals) :: dz, kext, dtau, mu_top, mu_side
      real(ireals) :: area, expon, sin_alpha, transport

      integer(mpiint) :: ierr

      call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetSection(plex%edir_dm, s, ierr); call CHKERR(ierr)
      call VecGetArrayF90(lambertVec, xv, ierr); call CHKERR(ierr)

      call get_vertical_cell_idx(plex%geom_dm, itopcell, plex%Nz-i1, cell_idx)
      !print *,''
      !print *,'----------------------------------------------'
      !print *,'itopcell', itopcell, 'cell_idx', cell_idx

      call get_top_bot_face_of_cell(plex%edir_dm, itopcell, iface_top, iface_bot)

      call get_inward_face_normal(iface_top, itopcell, geomSection, geoms, face_normal_top)
      call get_inward_face_normal(iface_bot, itopcell, geomSection, geoms, face_normal)
      !print *,'iface_top, iface_bot', iface_top, iface_bot,'normal_top', face_normal_top,'normal_bot',face_normal

      if(is_solar_src(face_normal_top, sundir)) then ! if sun actually shines on top of the column
        mu_top = dot_product(sundir, face_normal_top)

        dtau = zero

        do k = 1, size(cell_idx)
          icell = cell_idx(k)
          call get_top_bot_face_of_cell(plex%edir_dm, icell, iface_top, iface_bot)

          call get_inward_face_normal(iface_top, icell, geomSection, geoms, face_normal_top)
          call get_inward_face_normal(iface_bot, icell, geomSection, geoms, face_normal)

          if(.not.is_solar_src(-face_normal, sundir)) then
            cycle ! if the bot face is not sunlit... we dont allow the sun to shine from below
          endif

          call PetscSectionGetOffset(s, iface_top, voff, ierr); call CHKERR(ierr)

          dz = abs(plex%hhl(k) - plex%hhl(k+1))
          kext = (xkabs(i1+icell) + xksca(i1+icell))
          !print *,'k',k,'tau increment', icell, optprop(icell)%kabs , optprop(icell)%ksca , dz, &
          !        '->', dtau + kext * dz

          call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

          do i=1,size(faces_of_cell)
            iface_side = faces_of_cell(i)
            if(.not.plex%ltopfacepos(iface_side))then
              call get_inward_face_normal(iface_side, icell, geomSection, geoms, face_normal)

              if(is_solar_src(face_normal, sundir)) then
                mu_side = dot_product(sundir, face_normal)

                call PetscSectionGetOffset(geomSection, iface_side, geom_offset, ierr); call CHKERR(ierr)
                area = geoms(geom_offset+i7)

                ! integrate exp((-k*z)/s)dz for z=0 to h
                sin_alpha = sqrt(one - mu_side**2)
                expon = min(EXP_MAXVAL, max(EXP_MINVAL, dz * kext / sin_alpha))
                transport = (sin_alpha - sin_alpha * exp(-expon)) / kext / dz

                expon = min(EXP_MAXVAL, max(EXP_MINVAL, dtau / sin_alpha))
                transport = transport * exp(-expon)

                call PetscSectionGetOffset(s, iface_side, voff, ierr); call CHKERR(ierr)
                xv(voff+i1) = E0 * transport * area * mu_top * mu_side
                !print *,'solar_src', iface_side, E0, dtau, transport, mu_side, area, &
                !  '->', E0 * transport * mu_side, xv(voff+i1)
              endif
            endif
          enddo
          call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

          dtau = dtau + kext * dz
        enddo
      endif

      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(lambertVec, xv, ierr); call CHKERR(ierr)
    end subroutine

    subroutine solve_plex_rt(dm, b, A, ksp, x, prefix)
      type(tDM), intent(inout) :: dm
      type(tVec), allocatable, intent(in) :: b
      type(tMat), allocatable, intent(in) :: A
      type(tKSP), allocatable, intent(inout) :: ksp
      type(tVec), allocatable, intent(inout) :: x
      character(len=*),optional :: prefix

      type(tMatNullSpace) :: nullspace
      type(tVec) :: nullvecs(0)

      integer(mpiint) :: comm, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)

      !if(.not.allocated(ksp)) call CHKERR(1_mpiint, 'KSP has to be allocated before running solve')
      if(.not.allocated(b)) call CHKERR(1_mpiint, 'Src Vector has to be allocated before running solve')
      if(.not.allocated(A)) call CHKERR(1_mpiint, 'System Matrix has to be allocated before running solve')

      if(ldebug) print *,'plex_rt::solve Matrix...'
      if(.not.allocated(ksp)) then
        allocate(ksp)
        call KSPCreate(comm, ksp, ierr); call CHKERR(ierr)
        if(present(prefix)) then
          call KSPAppendOptionsPrefix(ksp, trim(prefix), ierr); call CHKERR(ierr)
        endif
      endif

      call MatNullSpaceCreate(comm, PETSC_TRUE, i0, nullvecs, nullspace, ierr) ; call CHKERR(ierr)
      call MatSetNearNullSpace(A, nullspace, ierr);call CHKERR(ierr)
      call MatNullSpaceDestroy(nullspace, ierr); call CHKERR(ierr)

      call KSPSetDM(ksp, dm, ierr); call CHKERR(ierr)
      call KSPSetDMActive(ksp, PETSC_FALSE, ierr); call CHKERR(ierr)
      call KSPSetOperators(ksp, A, A, ierr); CHKERRQ(ierr)

      call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
      call KSPSetUp(ksp, ierr); call CHKERR(ierr)

      if(.not.allocated(x)) then
        allocate(x)
        call VecDuplicate(b, x, ierr); CHKERRQ(ierr)
      endif
      call KSPSolve(ksp, b, x, ierr); CHKERRQ(ierr)
      if(ldebug) print *,'plex_rt::solve Matrix...finished'
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
    integer(iintegers) :: geom_offset, face_offset, iface, idof, num_dof

    real(ireals) :: area

    integer(mpiint) :: myid, ierr

    if(ldebug) print *,'plex_rt::scale_facevec...'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(face_dm, faceSection, ierr); call CHKERR(ierr)

    call DMGetLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)

    call DMGlobalToLocalBegin(face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

    call VecGetArrayF90(faceVec, xv, ierr); call CHKERR(ierr)

    do iface = plex%fStart, plex%fEnd-1
      call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
      area = geoms(geom_offset+i7)

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
    if(ldebug) print *,'plex_rt::scale_facevec... finished'
    end subroutine

  subroutine create_edir_mat(plex, OPP, kabs, ksca, g, A)
    type(t_plexgrid), intent(inout) :: plex
    class(t_optprop), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g
    type(tMat), allocatable, intent(inout) :: A

    type(tPetscSection) :: sec

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: icell, iface, irows(1), icols(1), idst
    real(ireals) :: coeffs(1)

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: wedge_offset

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)

    real(ireals) :: zenith, azimuth

    real(ireals) :: dir2dir(5), T(5), S(8)
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: zindex
    real(ireals) :: dz, coeff(5**2) ! coefficients for each src=[1..5] and dst[1..5]

    logical, parameter :: lonline=.False.

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...'

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

    if(.not.allocated(A)) then
      allocate(A)
      call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)
      !call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)
    endif

    do icell = plex%cStart, plex%cEnd-1
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)
      zenith  = wedgeorient(wedge_offset+i1)
      azimuth = wedgeorient(wedge_offset+i2)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
        lsrc(iface) = nint(wedgeorient(wedge_offset+i8+iface), iintegers) .eq. i1
      enddo
      !print *,'face_plex2bmc', face_plex2bmc

      zindex = plex%zindex(icell)
      dz = plex%hhl(zindex) - plex%hhl(zindex+1)

      !print *,'icell',icell,': foc',faces_of_cell
      !print *,'icell',icell,':',face_plex2bmc,'lsrc',lsrc
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
        dz, wedgeorient(wedge_offset+14:wedge_offset+19), .True., coeff, angles=[rad2deg(azimuth), rad2deg(zenith)])

      !do isrc = 1, 5
      !  print *,'LUT',rad2deg(azimuth), rad2deg(zenith), 'src '//itoa(isrc), coeff(isrc:size(coeff):i5)
      !enddo
      do iface = 1, size(faces_of_cell)

        if(lsrc(iface)) then
          ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering

          !dir2dir([upper_face, base_face, left_face, right_face, bottom_face]) = coeff((iwedgeface-1)*i5+i1:iwedgeface*i5)
          dir2dir = coeff(face_plex2bmc(iface):size(coeff):i5)
          if(lonline) then
            !print *,'lonline LUT0', coeff(face_plex2bmc(iface):size(coeff):i5)
            !call compute_dir2dir_coeff(face_plex2bmc(iface), dz, &
            !  xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
            !  rad2deg(azimuth), rad2deg(zenith), S, T)
            !print *,'lonline BOX1', T
            call compute_dir2dir_coeff(face_plex2bmc(iface), dz, &
              xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
              rad2deg(azimuth), rad2deg(zenith), S, T, wedgeorient(wedge_offset+14:wedge_offset+19))
            print *,'lonline BOX0', dir2dir
            print *,'lonline BOX2', T
            dir2dir = T
          endif

          call PetscSectionGetOffset(sec, faces_of_cell(iface), icols(1), ierr); call CHKERR(ierr)

          do idst = 1, size(faces_of_cell)
            if(.not.lsrc(idst)) then
              coeffs(1) = -dir2dir(face_plex2bmc(idst))
              if(coeffs(1).lt.zero) then
                call PetscSectionGetOffset(sec, faces_of_cell(idst), irows(1), ierr); call CHKERR(ierr)
                !print *,'isrc', face_plex2bmc(iface), 'idst', face_plex2bmc(idst), &
                !  'if', iface, 'id', idst, &
                !  'srcface->dstface', faces_of_cell(iface),faces_of_cell(idst), &
                !  'col -> row', icol, irow, dir2dir(face_plex2bmc(idst))
                call MatSetValuesLocal(A, i1, irows, i1, icols, coeffs, INSERT_VALUES, ierr); call CHKERR(ierr)
              endif
            endif
          enddo

        else ! This is not a src face, no radiation comes from here
        endif

      enddo ! enddo iface

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    !call MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
    !call MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)

    ! Set Diagonal Entries
    do iface = plex%fStart, plex%fEnd-1
      call PetscSectionGetOffset(sec, iface, irows(1), ierr); call CHKERR(ierr)
      call MatSetValuesLocal(A, i1, irows, i1, irows, [one], INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Medir', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...finished'
  end subroutine

  subroutine create_ediff_mat(plex, OPP, kabs, ksca, g, albedo, A)
    type(t_plexgrid), intent(in) :: plex
    class(t_optprop), intent(in) :: OPP
    type(tVec), allocatable, intent(in) :: kabs, ksca, g ! cell1_dm
    type(tVec), allocatable, intent(in) :: albedo        ! srfc_boundary_dm
    type(tMat), allocatable, intent(inout) :: A

    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: xkabs(:), xksca(:), xg(:)
    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: i, j, icell, iface, irows(1), icols(1)
    real(ireals) :: coeffs(1)

    type(tPetscSection) :: ediffSection, geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: wedge_offset

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    ! face_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. face_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: face_plex2bmc(5)
    integer(iintegers) :: diff_plex2bmc(8)

    real(ireals) :: diff2diff(8)

    integer(iintegers) :: zindex
    real(ireals) :: dz, coeff(8**2) ! coefficients for each src=[1..8] and dst[1..8]

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_ediff_mat...'

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

    if(.not.allocated(A)) then
      allocate(A)
      call DMCreateMatrix(plex%ediff_dm, A, ierr); call CHKERR(ierr)
      !call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)
    endif

    do icell = plex%cStart, plex%cEnd-1
      call DMPlexGetCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      call determine_diff_incoming_outgoing_offsets(plex, plex%ediff_dm, icell, incoming_offsets, outgoing_offsets)
      !print *,'icell', icell, 'faces_of_cell', faces_of_cell
      !print *,'icell', icell, 'incoming_offsets', incoming_offsets
      !print *,'icell', icell, 'outgoing_offsets', outgoing_offsets

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)

      do iface = 1, size(faces_of_cell)
        face_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
      enddo
      call face_idx_to_diff_bmc_idx(face_plex2bmc, diff_plex2bmc)

      zindex = plex%zindex(icell)
      dz = plex%hhl(zindex) - plex%hhl(zindex+1)

      !print *,'icell',icell,': foc',faces_of_cell
      !print *,'icell',icell,':',face_plex2bmc
      call get_coeff(OPP, xkabs(i1+icell), xksca(i1+icell), xg(i1+icell), &
        dz, wedgeorient(wedge_offset+14:wedge_offset+19), .False., coeff)

      !do isrc = 1, 8
      !  print *,'LUT src '//itoa(isrc), coeff(isrc:size(coeff):i8)
      !enddo
      do i = 1, size(incoming_offsets)
        icols = incoming_offsets(i)

        ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering
        diff2diff = coeff(diff_plex2bmc(i):size(coeff):i8)

        do j = 1, size(outgoing_offsets)
          coeffs(1) = -diff2diff(diff_plex2bmc(j))
          if(coeffs(1).lt.zero) then
            irows(1) = outgoing_offsets(j)
            !print *,'icell',icell,'i,j',i,j,'icol', icol, 'irow', irow, '=>', diff2diff(diff_plex2bmc(j))
            call MatSetValuesLocal(A, i1, irows, i1, icols, coeffs, INSERT_VALUES, ierr); call CHKERR(ierr)
          endif
        enddo

      enddo ! enddo iface

      call DMPlexRestoreCone(plex%ediff_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    call set_boundary_conditions()

    call set_diagonal_entries()

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Mediff', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_ediff_mat...finished'
    contains
      subroutine set_boundary_conditions()
        type(tIS) :: bc_ids
        integer(iintegers), pointer :: xi(:)
        integer(iintegers) :: i, iface, offset_Ein, offset_Eout, offset_srfc, idof, numDof
        type(tPetscSection) :: srfcSection

        real(ireals), pointer :: xalbedo(:)

        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', BOTFACE, bc_ids, ierr); call CHKERR(ierr)
        if (bc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
        else
          call DMGetSection(plex%srfc_boundary_dm, srfcSection, ierr); CHKERRQ(ierr)
          call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
          call ISGetIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          do i = 1, size(xi)
            iface = xi(i)
            ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
            ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
            call PetscSectionGetFieldOffset(ediffSection, iface, i0, offset_Ein, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(ediffSection, iface, i1, offset_Eout, ierr); call CHKERR(ierr)

            call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

            call MatSetValuesLocal(A, i1, [offset_Ein], i1, [offset_Eout], [-xalbedo(i1+offset_srfc)], INSERT_VALUES, ierr); call CHKERR(ierr)
          enddo
          call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
        endif

        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', SIDEFACE, bc_ids, ierr); call CHKERR(ierr)
        if (bc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
        else
          call ISGetIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
          do i = 1, size(xi)
            iface = xi(i)
            ! field offset for field 0 gives Eup because field 0 is flux from lower cell id to higher cell id
            ! in case of boundary faces: from cell_id -1 (boundary face) to some icell
            call PetscSectionGetFieldOffset(ediffSection, iface, i0, offset_Ein, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldOffset(ediffSection, iface, i1, offset_Eout, ierr); call CHKERR(ierr)
            call PetscSectionGetFieldDof(ediffSection, iface, i0, numDof, ierr); call CHKERR(ierr)

            call PetscSectionGetOffset(srfcSection, iface, offset_srfc, ierr); call CHKERR(ierr)

            do idof = 0, numDof-1
              call MatSetValuesLocal(A, i1, [offset_Ein+idof], i1, [offset_Eout+idof], [-one], INSERT_VALUES, ierr); call CHKERR(ierr)
            enddo
          enddo
          call ISRestoreIndicesF90(bc_ids, xi, ierr); call CHKERR(ierr)
        endif
      end subroutine
      subroutine set_diagonal_entries()
        integer(iintegers) :: irow
        call PetscSectionGetOffsetRange(ediffSection, i, j, ierr); call CHKERR(ierr)
        do irow = i, j-1
          call MatSetValuesLocal(A, i1, [irow], i1, [irow], [one], INSERT_VALUES, ierr); call CHKERR(ierr)
        enddo
      end subroutine
  end subroutine
  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(OPP, kabs, ksca, g, dz, wedge_coords, ldir, coeff, angles)
    class(t_optprop), intent(in) :: OPP
    real(ireals), intent(in)     :: kabs, ksca, g
    real(ireals),intent(in)      :: dz
    real(ireals), intent(in)     :: wedge_coords(:) ! coordinates of upper triangle pts A,B,C in in (x,y)
    logical,intent(in)           :: ldir
    real(ireals),intent(out)     :: coeff(:)

    real(ireals),intent(in),optional  :: angles(2)

    real(ireals) :: aspect, tauz, w0, dx, relcoords(6)
    integer, parameter :: iC1=4, iC2=5

    dx = wedge_coords(3)

    aspect = dz / dx
    tauz = (kabs+ksca) * dz
    if(approx(tauz,zero)) then
      w0 = zero
    else
      w0 = ksca / (kabs+ksca)
    endif

    relcoords = wedge_coords / dx

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

    !relcoords(5) = max(OPP%OPP_LUT%diffconfig%dims(iC1)%vrange(1), &
    !               min(OPP%OPP_LUT%diffconfig%dims(iC1)%vrange(2), relcoords(5)))
    !relcoords(6) = max(OPP%OPP_LUT%diffconfig%dims(iC2)%vrange(1), &
    !               min(OPP%OPP_LUT%diffconfig%dims(iC2)%vrange(2), relcoords(6)))

    !print *,'DEBUG Lookup Coeffs for', tauz, w0, g, aspect, angles, ':', norm(wedge_coords(3:4)-wedge_coords(1:2)), ':', relcoords
    call OPP%get_coeff(tauz, w0, g, aspect, ldir, coeff, angles=angles, wedge_coords=relcoords)
    !print *,'DEBUG Lookup Coeffs for', tauz, w0, g, aspect, angles, ':', norm(wedge_coords(3:4)-wedge_coords(1:2)), ':', relcoords, '::', coeff
    if(ldebug) then
      if(any(coeff.lt.zero).or.any(coeff.gt.one)) then
        print *,'Lookup Coeffs for', aspect, tauz, w0, g, angles,'::', coeff
        call CHKERR(1_mpiint, 'Found corrupted coefficients!')
      endif
    endif

  end subroutine

  subroutine compute_dir2dir_coeff(src, dz, kabs, ksca, g, phi, theta, S, T, coord2d)
    use m_boxmc, only : t_boxmc, t_boxmc_wedge_5_5, t_boxmc_wedge_5_8
    use m_boxmc_geometry, only : setup_default_wedge_geometry
    integer(iintegers), intent(in) :: src
    real(ireals), intent(in) :: dz, kabs, ksca, g, phi, theta
    real(ireals), intent(out) :: S(:),T(:)
    real(ireals), intent(in), optional :: coord2d(:)

    type(t_boxmc_wedge_5_8) :: bmc_wedge
    real(ireals) :: bg(3), dx
    real(ireals), allocatable :: vertices(:)
    real(ireals) :: S_tol(size(S)),T_tol(size(T))

    call bmc_wedge%init(PETSC_COMM_SELF)
    if(present(coord2d)) then
      print *,'computing coeffs for src/phi/theta',src,phi,theta,':',coord2d
    else
      print *,'computing coeffs for src/phi/theta',src,phi,theta
    endif

    !bg  = [kabs, ksca, g]

    call setup_default_wedge_geometry(coord2d(1:2), coord2d(3:4), coord2d(5:6), dz, vertices)

    dx = vertices(4)
    vertices = vertices/dx
    bg  = [kabs*dx, ksca*dx, g]

    call bmc_wedge%get_coeff(PETSC_COMM_SELF, bg, src, .True., &
      phi, theta, vertices, S, T, S_tol, T_tol, inp_atol=2e-4_ireals, inp_rtol=1e-2_ireals)
  end subroutine

  subroutine restore_solution(solver, solution, time)
    type(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution
    real(ireals),intent(in),optional :: time

    if(present(time)) then
      print *,'Adaptive Time Integration not yet fully implemented'
    endif
    call compute_absorption(solver, solution)
    solution%lchanged = .False.
  end subroutine

  subroutine compute_absorption(solver, solution)
    type(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution

    integer(mpiint) :: ierr

    if(.not.solution%lset) call CHKERR(1_mpiint, 'compute_absorption needs to get an initialized solution obj')

    if(.not.solution%lchanged) return
    if(solution%lsolar_rad .and. .not.allocated(solution%edir)) &
      call CHKERR(1_mpiint, 'solution%lsolar_rad true but edir not allocated')
    if(.not.allocated(solution%ediff)) call CHKERR(1_mpiint, 'ediff vec not allocated')

    ! Make sure we have the radiation vecs in plain energy units
    call scale_flx(solver, solution, lWm2=.False.)

    if(.not.allocated(solution%abso)) then
      allocate(solution%abso)
      call DMCreateGlobalVector(solver%plex%abso_dm, solution%abso, ierr); call CHKERR(ierr)
      call PetscObjectSetName(solution%abso, 'absoVecGlobal', ierr);call CHKERR(ierr)
    endif
    call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad) then
      call compute_edir_absorption(solver%plex, solution%edir, solution%abso)
    endif
    call compute_ediff_absorption(solver%plex, solution%ediff, solution%abso)

    call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
  end subroutine

  subroutine compute_edir_absorption(plex, edir, abso)
    type(t_plexgrid), allocatable, intent(inout) :: plex
    type(tVec), intent(in) :: edir
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_edir, local_abso

    real(ireals), pointer :: xedir(:), xabso(:)

    type(tPetscSection) :: abso_section, edir_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface
    integer(iintegers),pointer :: faces_of_cell(:)
    integer(mpiint) :: myid, ierr

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:), wedgeorient(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset, edir_offset, wedge_offset

    real(ireals) :: volume

    logical :: lsrc ! is src or destination of solar beam

    if(.not.allocated(plex)) stop 'called compute_edir_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%edir_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_edir_absorption with an unallocated abso vec')

    if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMGetSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!

    call DMGetLocalVector(plex%abso_dm, local_abso, ierr); call CHKERR(ierr)
    call VecSet(local_abso, zero, ierr); call CHKERR(ierr)

    call DMGetLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
      volume = geoms(geom_offset+i4)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      xabso(abso_offset+i1) = zero

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)

      do iface = 1, size(faces_of_cell)
        lsrc = int(wedgeorient(wedge_offset+i8+iface), iintegers) .eq. i1

        call PetscSectionGetOffset(edir_section, faces_of_cell(iface), edir_offset, ierr); call CHKERR(ierr)

        if(lsrc) then ! sun is shining into this face
          xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1)
        else
          xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1)
        endif
      enddo
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      xabso(abso_offset+i1) = xabso(abso_offset+i1) / volume
    enddo

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%edir_dm, local_edir, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%abso_dm, local_abso, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_edir_abso', ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end subroutine

  subroutine compute_ediff_absorption(plex, ediff, abso)
    type(t_plexgrid), allocatable, intent(inout) :: plex
    type(tVec), intent(in) :: ediff
    type(tVec), allocatable, intent(inout) :: abso

    type(tVec) :: local_ediff, local_abso

    real(ireals), pointer :: xediff(:), xabso(:)

    type(tPetscSection) :: abso_section, ediff_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface
    integer(mpiint) :: myid, ierr

    integer(iintegers), allocatable :: incoming_offsets(:), outgoing_offsets(:)

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates and orientation vec
    integer(iintegers) :: geom_offset, abso_offset

    real(ireals) :: volume

    if(.not.allocated(plex)) stop 'called compute_ediff_absorption but plex is not allocated'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    if(.not.allocated(plex%ediff_dm)) call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+1, 'called compute_ediff_absorption with a dm which is not allocated?')
    if(.not.allocated(abso)) call CHKERR(myid+1, 'called compute_ediff_absorption with an unallocated abso vec')

    if(ldebug) print *,'plex_rt::compute_ediff_absorption....'

    call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
    call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetLocalVector(plex%abso_dm, local_abso, ierr); call CHKERR(ierr)

    call DMGetLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)
    call VecSet(local_ediff, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%ediff_dm, ediff, INSERT_VALUES, local_ediff, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecGetArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call determine_diff_incoming_outgoing_offsets(plex, plex%ediff_dm, icell, incoming_offsets, outgoing_offsets)

      call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
      volume = geoms(geom_offset+i4)

      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
      xabso(i1+abso_offset) = zero

      do iface = 1, size(incoming_offsets)
        xabso(i1+abso_offset) = xabso(i1+abso_offset) + xediff(i1+incoming_offsets(iface))
        xabso(i1+abso_offset) = xabso(i1+abso_offset) + xediff(i1+outgoing_offsets(iface))
      enddo

      xabso(abso_offset+i1) = xabso(abso_offset+i1) / volume
    enddo

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%ediff_dm, local_ediff, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%abso_dm, local_abso, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_ediff_abso', ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::compute_ediff_absorption....finished'
  end subroutine


  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solver, solution, lWm2)
    class(t_plex_solver), intent(inout)   :: solver
    type(t_state_container),intent(inout) :: solution  !< @param solution container with computed fluxes
    logical,intent(in)                    :: lWm2      !< @param determines direction of scaling, if true, scale to W/m**2
    integer(mpiint) :: ierr

    if(solution%lsolar_rad) then
      if(.not.allocated(solver%dir_scalevec_Wm2_to_W)) then
        allocate(solver%dir_scalevec_Wm2_to_W)
        call VecDuplicate(solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        call VecSet(solver%dir_scalevec_Wm2_to_W, one, ierr); call CHKERR(ierr)
        call scale_facevec(solver%plex, solver%plex%edir_dm, solver%dir_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
      endif
      if(.not.allocated(solver%dir_scalevec_W_to_Wm2)) then
        allocate(solver%dir_scalevec_W_to_Wm2)
        call VecDuplicate(solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        call VecSet(solver%dir_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
        call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      endif
      if(solution%lWm2_dir .neqv. lWm2) then
        if(lWm2) then
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        else
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        endif
        solution%lWm2_dir = lWm2
      endif
    endif

    if(.not.allocated(solver%diff_scalevec_Wm2_to_W)) then
      allocate(solver%diff_scalevec_Wm2_to_W)
      call VecDuplicate(solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      call VecSet(solver%diff_scalevec_Wm2_to_W, one, ierr); call CHKERR(ierr)
      call scale_facevec(solver%plex, solver%plex%ediff_dm, solver%diff_scalevec_Wm2_to_W, lW_to_Wm2=.False.)
    endif
    if(.not.allocated(solver%diff_scalevec_W_to_Wm2)) then
      allocate(solver%diff_scalevec_W_to_Wm2)
      call VecDuplicate(solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      call VecSet(solver%diff_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
      call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
    endif

    if(solution%lWm2_diff .neqv. lWm2) then
      if(lWm2) then
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      else
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      endif
      solution%lWm2_diff = lWm2
    endif
  end subroutine

  subroutine plexrt_get_result(solver, redn, reup, rabso, redir, opt_solution_uid)
    type(t_plex_solver), intent(inout) :: solver
    real(ireals), allocatable, dimension(:,:), intent(inout) :: redn,reup,rabso
    real(ireals), allocatable, dimension(:,:), intent(inout), optional :: redir
    integer(iintegers),optional,intent(in) :: opt_solution_uid

    type(tIS) :: toa_ids
    integer(iintegers) :: Ncol, ke1, uid, i, k, iface, icell, voff
    integer(iintegers), pointer :: xtoa_faces(:), cell_support(:)

    type(tPetscSection) :: abso_section, ediff_section, edir_section
    real(ireals), pointer :: xediff(:), xedir(:), xabso(:)

    integer(mpiint) :: myid, ierr

    uid = get_arg(0_iintegers, opt_solution_uid)
    if(.not.solver%solutions(uid)%lset) &
      call CHKERR(1_mpiint, 'You tried to retrieve results from a solution uid which has not yet been calculated')

    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

    ke1 = solver%plex%Nz

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

    call scale_flx(solver, solution, lWm2=.True.)

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
        call PetscSectionGetFieldOffset(ediff_section, iface+k, i1, voff, ierr); call CHKERR(ierr)
        reup(i1+k, i) = xediff(i1+voff)
      enddo
      ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
      call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i1, voff, ierr); call CHKERR(ierr)
      redn(i1+k, i) = xediff(i1+voff)
      call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
      reup(i1+k, i) = xediff(i1+voff)

      ! Fill Aborption Vec
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
      endif
    enddo

    if(present(redir)) then
      call VecRestoreArrayF90(solution%edir , xedir , ierr); call CHKERR(ierr)
    endif
    call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)
    end associate
    if(ldebug) then
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)
      if(myid.eq.0) then
        if(present(redir)) then
          print *,'Get Result, k        Edir                Edn                 Eup                abso'
          do k = 1, ke1-1
            print *,k, redir(k,1), redn(k,1), reup(k,1), rabso(k,1)
          enddo
          print *,k, redir(k,1), redn(k,1), reup(k,1)
        else
          print *,'Get Result, k        Edn                 Eup                abso'
          do k = 1, ke1-1
            print *,k, redn(k,1), reup(k,1), rabso(k,1)
          enddo
          print *,k, redn(k,1), reup(k,1)
        endif
      endif
    endif

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
