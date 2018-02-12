module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options

  use m_helper_functions, only: CHKERR, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, strF2C, &
    vec_proj_on_plane, cross_3d, norm, rotation_matrix_world_to_local_basis, &
    imp_bcast, approx

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, i6, i7, &
    zero, one, pi, EXP_MINVAL, EXP_MAXVAL

  use m_plex_grid, only: t_plexgrid, compute_face_geometry, setup_edir_dmplex, setup_abso_dmplex, &
    facevec2cellvec, icell_icon_2_plex, iface_top_icon_2_plex

  use m_optprop, only : t_optprop, t_optprop_wedge_5_8
  use m_optprop_parameters, only : OPP_LUT_ALL_ANGLES

  implicit none

  private

  public :: t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, &
    get_normal_of_first_TOA_face, compute_face_geometry, &
    set_plex_rt_optprop

  type t_opticalprops
    real(ireals) :: kabs,ksca,g
  end type

  type t_plex_solver
    type(t_plexgrid), allocatable :: plex
    class(t_optprop), allocatable :: OPP
    type(t_opticalprops), allocatable, dimension(:) :: optprop ! in each cell [pStart+1..pEnd]
  end type

  logical, parameter :: ldebug=.True.
  contains
    subroutine init_plex_rt_solver(plex, solver)
      type(t_plexgrid), intent(in) :: plex
      type(t_plex_solver), allocatable, intent(inout) :: solver

      call read_commandline_options()

      if(allocated(solver)) stop 'Should not call init_plex_rt_solver with already allocated solver object'
      allocate(solver)

      allocate(solver%plex)
      solver%plex = plex

      allocate(t_optprop_wedge_5_8::solver%OPP)
      call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], [real(OPP_LUT_ALL_ANGLES, kind=ireals)], solver%plex%comm)
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], real([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70],ireals), plex%comm)
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], [zero, 5*one], plex%comm)
    end subroutine

    subroutine set_plex_rt_optprop(solver, vlwc, viwc)
      use m_helper_functions, only : delta_scale
      type(t_plex_solver), allocatable, intent(inout) :: solver
      type(tVec),intent(in), optional :: vlwc, viwc
      real(ireals), pointer :: xlwc(:), xiwc(:)

      real(ireals),parameter :: w0 = .95, reff_w=10, reff_i=30, rayleigh=1e-5
      real(ireals) :: ksca_tot, kabs_tot, g_tot
      real(ireals) :: ksca_cld, kabs_cld, g_cld

      integer(iintegers) :: i
      integer(mpiint) :: ierr

      if(.not.allocated(solver)) stop 'set_plex_rt_optprop::solver has to be allocated'

      if(.not.allocated(solver%optprop)) allocate(solver%optprop(solver%plex%cStart+i1:solver%plex%cEnd))

      if(present(vlwc)) call VecGetArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      if(present(viwc)) call VecGetArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)

      do i = 1,size(solver%optprop)
        kabs_tot = rayleigh/10
        ksca_tot = rayleigh
        g_tot    = zero

        if(present(vlwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xlwc(i) * 1e3) / reff_w * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xlwc(i) * 1e3) / reff_w * (w0)
          g_cld    = .85_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = g_cld * ksca_cld / (ksca_cld + ksca_tot)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        if(present(viwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xiwc(i) * 1e3) / reff_i * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xiwc(i) * 1e3) / reff_i * (w0)
          g_cld    = .85_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = (g_tot * ksca_tot + g_cld * ksca_cld) / (ksca_tot + ksca_cld)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        solver%optprop(i)%kabs = kabs_tot
        solver%optprop(i)%ksca = ksca_tot
        solver%optprop(i)%g    = g_tot
      enddo

      if(present(vlwc)) call VecRestoreArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      if(present(viwc)) call VecRestoreArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)

      print *,'Min/Max of kabs', minval(solver%optprop(:)%kabs), maxval(solver%optprop(:)%kabs)
      print *,'Min/Max of ksca', minval(solver%optprop(:)%ksca), maxval(solver%optprop(:)%ksca)
      print *,'Min/Max of g   ', minval(solver%optprop(:)%g   ), maxval(solver%optprop(:)%g   )

    end subroutine

    subroutine run_plex_rt_solver(solver, sundir)
      type(t_plex_solver), allocatable, intent(inout) :: solver
      real(ireals), intent(in) :: sundir(3) ! cartesian direction of sun rays, norm of vector is the energy in W/m2

      type(tVec) :: b, edir, abso
      type(tMat) :: Mdir
      integer(mpiint) :: ierr

      if(.not.allocated(solver)) stop 'run_plex_rt_solver::solver has to be allocated'

      if(.not.allocated(solver%plex%geom_dm)) stop 'run_plex_rt_solver::geom_dm has to be allocated first'
      if(.not.allocated(solver%optprop)) stop 'run_plex_rt_solver::optprop has to be allocated first'
      if(.not.allocated(solver%plex%geom_dm)) call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      if(.not.allocated(solver%plex%edir_dm)) call setup_edir_dmplex(solver%plex, solver%plex%edir_dm)

      call create_src_vec(solver%plex, solver%plex%edir_dm, norm(sundir), solver%optprop, sundir, b)

      ! Output of srcVec
      call scale_facevec(solver%plex, solver%plex%edir_dm, b, lW_to_Wm2=.True.)
      call facevec2cellvec(solver%plex, solver%plex%edir_dm, b)
      call scale_facevec(solver%plex, solver%plex%edir_dm, b, lW_to_Wm2=.False.)

      call create_edir_mat(solver%plex, solver%OPP, solver%optprop, sundir, Mdir)

      call solve_plex_rt(solver%plex, b, Mdir, edir)
      call PetscObjectSetName(edir, 'edir', ierr); call CHKERR(ierr)

      ! Output of Edir
      call scale_facevec(solver%plex, solver%plex%edir_dm, edir, lW_to_Wm2=.True.)
      call facevec2cellvec(solver%plex, solver%plex%edir_dm, edir)
      call scale_facevec(solver%plex, solver%plex%edir_dm, edir, lW_to_Wm2=.False.)

      call setup_abso_dmplex(solver%plex, solver%plex%abso_dm)
      call compute_edir_absorption(solver%plex, edir, sundir, abso)
    end subroutine

    subroutine create_src_vec(plex, edirdm, E0, optprop, sundir, srcVec)
      type(t_plexgrid), intent(in) :: plex
      type(tDM),allocatable, intent(in) :: edirdm
      real(ireals), intent(in) :: E0, sundir(3)
      type(t_opticalprops), intent(in) :: optprop(:)
      type(tVec),intent(out) :: srcVec

      integer(iintegers) :: i, voff
      type(tPetscSection) :: s

      type(tVec) :: localVec, lambertVec
      real(ireals), pointer :: xv(:), xlambert(:)

      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: fStart, fEnd

      type(tIS) :: boundary_ids
      logical :: ltopface
      integer(iintegers) :: iface, icell, k

      integer(iintegers), pointer :: xx_v(:)
      integer(iintegers), pointer :: cell_support(:)

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: area, mu, cell_center(3), face_center(3), face_normal(3)

      integer(mpiint) :: comm, myid, ierr

      if(.not.allocated(edirdm)) stop 'called create_src_vec but face_dm is not allocated'

      if(.not.allocated(plex%geom_dm)) stop 'get_normal_of_first_TOA_face::needs allocated geom_dm first'

      call PetscObjectGetComm(edirdm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if(ldebug .and. myid.eq.0) print *,myid,'plex_rt::create_src_vec....'


      call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetDefaultSection(edirdm, s, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(s, PETSC_NULL_SECTION, '-show_src_section', ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(edirdm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call DMPlexGetHeightStratum(edirdm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges


      ! Now lets get vectors!
      call DMCreateGlobalVector(edirdm, srcVec,ierr); call CHKERR(ierr)
      call PetscObjectSetName(srcVec, 'srcVecGlobal', ierr);call CHKERR(ierr)
      call VecSet(srcVec, zero, ierr); call CHKERR(ierr)

      call DMGetLocalVector(edirdm, localVec,ierr); call CHKERR(ierr)
      call VecSet(localVec, zero, ierr); call CHKERR(ierr)

      call VecGetArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      call DMGetStratumIS(edirdm, 'DomainBoundary', i1, boundary_ids, ierr); call CHKERR(ierr)
      if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
      else
        call PetscObjectViewFromOptions(boundary_ids, PETSC_NULL_IS, '-show_IS_Boundary_TOA', ierr); call CHKERR(ierr)

        call ISGetIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)

        ! First set the TOA boundary fluxes on faces
        do i = 1, size(xx_v)
            iface = xx_v(i)
            k = plex%zindex(iface)
            ltopface = plex%ltopfacepos(iface)
            if(k.eq.1 .and. ltopface) then
                call DMPlexGetSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
                icell = cell_support(1)
                call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
                call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
                cell_center = geoms(i1+geom_offset:i3+geom_offset)

                call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
                area = geoms(geom_offset+i7)
                face_center = geoms(geom_offset+i1:geom_offset+i3)
                face_normal = geoms(geom_offset+i4:geom_offset+i6)
                face_normal = face_normal * & ! Determine the inward normal vec for the face
                  real(determine_normal_direction(face_normal, face_center, cell_center), ireals)

                mu = dot_product(sundir, face_normal)
                if(mu.gt.zero) then
                  call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
                  xv(voff+1) = E0 * area * mu
                endif
            endif ! k.eq.1
        enddo
        call ISRestoreIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)
      endif ! TOA boundary ids

      ! Then define sideways boundary conditions
      call DMGetStratumIS(edirdm, 'DomainBoundary', i2, boundary_ids, ierr); call CHKERR(ierr)
      if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have side boundary faces
      else
        call PetscObjectViewFromOptions(boundary_ids, PETSC_NULL_IS, '-show_IS_Boundary_SIDE', ierr); call CHKERR(ierr)

        call ISGetIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)

        call DMGetLocalVector(edirdm, lambertVec,ierr); call CHKERR(ierr)
        call PetscObjectSetName(lambertVec, 'srcVec_side_lambert', ierr);call CHKERR(ierr)
        call VecSet(lambertVec, zero, ierr); call CHKERR(ierr)

        do i = 1, size(xx_v)
            iface = xx_v(i)
            k = plex%zindex(iface)
            ltopface = plex%ltopfacepos(iface)
            if(k.eq.1.and..not.ltopface) then
                ! Then check if the cell below is a sideward boundary cell
                call DMPlexGetSupport(edirdm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell
                icell = cell_support(1)
                call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
                call compute_lambert_beer(plex, icell, E0, sundir, optprop, lambertVec)
            endif ! k.eq.1
        enddo

        call VecGetArrayReadF90(lambertVec, xlambert, ierr); call CHKERR(ierr)
        do i = 1, size(xx_v)
            iface = xx_v(i)
            ltopface = plex%ltopfacepos(iface)
            if(.not.ltopface) then
                call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
                area = geoms(geom_offset+7)

                call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
                xv(i1+voff) = xlambert(i1+voff) * area
            endif
        enddo
        call VecRestoreArrayReadF90(lambertVec, xlambert, ierr); call CHKERR(ierr)

        call ISRestoreIndicesF90(boundary_ids, xx_v, ierr); call CHKERR(ierr)
      endif

      call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

      call VecSet(srcVec, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(edirdm, localVec, ADD_VALUES, srcVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (edirdm, localVec, ADD_VALUES, srcVec, ierr); call CHKERR(ierr)
      call PetscObjectSetName(srcVec, 'srcVec', ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(srcVec, PETSC_NULL_VEC, '-show_src_vec_global', ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(edirdm, localVec, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      if(ldebug .and. myid.eq.0) print *,myid,'plex_rt::create_src_vec....finished'
      if(ldebug) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    end subroutine

    subroutine compute_lambert_beer(plex, itopcell, E0, sundir, optprop, lambertVec)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: itopcell
      real(ireals), intent(in) :: E0, sundir(3)
      type(t_opticalprops), intent(in) :: optprop(:)
      type(tVec),intent(in) :: lambertVec

      type(tPetscSection) :: s
      real(ireals), pointer :: xv(:)
      integer(iintegers), pointer :: faces_of_cell(:)
      integer(iintegers) :: i, iconcell, icell, iface, iface_top, iface_bot, k, voff

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: cell_center(3), face_normal(3), face_center(3)

      real(ireals) :: area, sidearea, dz, kext, dtau, m, mu, mu_side, Edir_top, transport

      integer(mpiint) :: ierr

      call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetDefaultSection(plex%edir_dm, s, ierr); call CHKERR(ierr)
      call VecGetArrayF90(lambertVec, xv, ierr); call CHKERR(ierr)

      iconcell = plex%localiconindex(itopcell)
      k = plex%zindex(itopcell)

      iface = iface_top_icon_2_plex(plex, iconcell, k)

      call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
      face_normal = geoms(geom_offset+4:geom_offset+6)
      area = geoms(geom_offset+7)

      mu = abs(dot_product(sundir, face_normal))

      Edir_top = E0

      call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
      xv(voff+1) = Edir_top

      do k = 1, plex%Nz
        iface_top = iface_top_icon_2_plex(plex, iconcell, k)
        iface_bot = iface_top_icon_2_plex(plex, iconcell, k+1)
        call PetscSectionGetOffset(s, iface_top, voff, ierr); call CHKERR(ierr)
        Edir_top = xv(voff+1)

        icell = icell_icon_2_plex(plex, iconcell, k)

        call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
        cell_center = geoms(1+geom_offset:3+geom_offset)

        call PetscSectionGetOffset(geomSection, iface_top, geom_offset, ierr); call CHKERR(ierr)
        face_normal = geoms(4+geom_offset: 6+geom_offset)
        mu = max(epsilon(mu), abs(dot_product(sundir, face_normal)))
        dz = plex%hhl(k) - plex%hhl(k+1)
        kext = (optprop(icell+1)%kabs + optprop(icell+1)%ksca)
        dtau = max(EXP_MINVAL, min(EXP_MAXVAL, kext * dz / mu))

        call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

        do i=1,size(faces_of_cell)
          iface = faces_of_cell(i)
          if(iface.eq.iface_top) cycle

          call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
          if(plex%ltopfacepos(iface))then
            xv(voff+1) = Edir_top * exp(-dtau)
            !print *,'lambert: itopcell',itopcell,':',iconcell, '::', k, xv(voff+1)
          else
            call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
            face_center = geoms(1+geom_offset: 3+geom_offset)
            face_normal = geoms(4+geom_offset: 6+geom_offset)
            sidearea = geoms(geom_offset+7)

            ! Determine the inward normal vec for the face
            face_normal = face_normal * &
              real(determine_normal_direction(face_normal, face_center, cell_center), ireals)

            ! average transmission on the slanted paths
            ! (integrate exp(-k * x / m) for x=0 to Z)/ Z
            mu_side = dot_product(sundir, face_normal)
            if(mu_side.gt.zero) then
              m = sqrt(one-mu_side**2)
              dtau = max(EXP_MINVAL, min(EXP_MAXVAL, kext * dz / m))
              transport = (one - exp(-dtau)) / ( dtau )
              xv(voff+1) = Edir_top * mu_side * transport
            endif
          endif
        enddo
        call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      enddo
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(lambertVec, xv, ierr); call CHKERR(ierr)
    end subroutine

    subroutine solve_plex_rt(plex, b, A, x)
      type(t_plexgrid) :: plex
      type(tVec), intent(in) :: b
      type(tVec), intent(inout) :: x
      type(tMat), intent(in) :: A

      type(tKSP) :: ksp
      integer(mpiint) :: ierr

      if(ldebug) print *,'plex_rt::solve Matrix...'
      call KSPCreate(plex%comm, ksp, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, ierr); CHKERRQ(ierr)

      call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
      call KSPSetUp(ksp, ierr); CHKERRQ(ierr)

      call VecDuplicate(b, x, ierr); CHKERRQ(ierr)
      call KSPSolve(ksp, b, x, ierr); CHKERRQ(ierr)
      call KSPDestroy(ksp, ierr); CHKERRQ(ierr)
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
    integer(iintegers) :: geom_offset, face_offset, iface

    real(ireals) :: area

    integer(mpiint) :: ierr

    if(ldebug) print *,'plex_rt::scale_facevec...'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(face_dm, faceSection, ierr); call CHKERR(ierr)

    call DMGetLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)

    call DMGlobalToLocalBegin(face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

    call VecGetArrayF90(faceVec, xv, ierr); call CHKERR(ierr)

    do iface = plex%fStart, plex%fEnd-1
      call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
      area = geoms(geom_offset+7)

      call PetscSectionGetOffset(faceSection, iface, face_offset, ierr); call CHKERR(ierr)
      if(lW_to_Wm2) then
        xv(i1+face_offset) = xv(i1+face_offset) / area
      else
        xv(i1+face_offset) = xv(i1+face_offset) * area
      endif
    enddo

    call VecRestoreArrayF90(faceVec, xv, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call VecSet(globalfaceVec, zero, ierr); call CHKERR(ierr)
    call DMLocalToGlobalBegin(face_dm, faceVec, ADD_VALUES, globalfaceVec, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (face_dm, faceVec, ADD_VALUES, globalfaceVec, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::scale_facevec... finished'
    end subroutine

  subroutine compute_edir_absorption(plex, edir, sundir, abso)
    type(t_plexgrid), intent(inout) :: plex
    type(tVec), intent(in) :: edir
    real(ireals), intent(in) :: sundir(3)
    type(tVec), intent(out) :: abso

    type(tVec) :: local_edir

    real(ireals), pointer :: xedir(:), xabso(:)

    type(tPetscSection) :: abso_section, edir_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface
    integer(iintegers),pointer :: faces_of_cell(:)
    integer(mpiint) :: ierr

    type(tPetscSection) :: geomSection
    real(ireals) :: cell_center(3), face_normal(3), face_center(3)
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    integer(iintegers) :: geom_offset, abso_offset, edir_offset

    real(ireals) :: mu, volume

    if(.not.allocated(plex%edir_dm) .or. .not.allocated(plex%abso_dm)) stop 'called compute_edir_absorption with a dm which is not allocated?'

    if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetDefaultSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetGlobalVector(plex%abso_dm, abso,ierr); call CHKERR(ierr)
    call PetscObjectSetName(abso, 'absoVecGlobal', ierr);call CHKERR(ierr)
    call VecSet(abso, zero, ierr); call CHKERR(ierr)

    call DMGetLocalVector(plex%edir_dm, local_edir,ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1
      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)

      call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
      cell_center = geoms(geom_offset+i1:geom_offset+i3)
      volume = geoms(geom_offset+i4)

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(edir_section, faces_of_cell(iface), edir_offset, ierr); call CHKERR(ierr)

        call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); call CHKERR(ierr)
        face_center = geoms(geom_offset+i1: geom_offset+i3)
        face_normal = geoms(geom_offset+i4: geom_offset+i6)

        ! Determine the inward normal vec for the face
        face_normal = face_normal * &
          real(determine_normal_direction(face_normal, face_center, cell_center), ireals)

        ! Then determine if the face is src(in line with the sun vec) or if it is destination(contra sun direction)
        mu = dot_product(face_normal, sundir)

        if(mu.gt.zero) then ! sun is shining into this face
          xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1)
        else
          xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1)
        endif
      enddo
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      xabso(abso_offset+i1) = xabso(abso_offset+i1) / volume
    enddo

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end subroutine

  subroutine create_edir_mat(plex, OPP, optprop, sundir, A)
    type(t_plexgrid), intent(inout) :: plex
    real(ireals), intent(in) :: sundir(3)
    class(t_optprop), intent(in) :: OPP
    type(t_opticalprops), intent(in) :: optprop(:)
    type(tMat), intent(out) :: A

    type(tPetscSection) :: sec

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(mpiint) :: myid, ierr

    integer(iintegers), pointer :: faces_of_cell(:), edges_of_face(:)
    integer(iintegers) :: i, icell, iface, iedge, irow, iwedgeface

    integer(iintegers) :: irows(5), icols(5)

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    integer(iintegers) :: geom_offset

    real(ireals) :: zenith, azimuth

    integer(iintegers) :: base_face   ! index of face which is closest to sun angle, index regarding faces_of_cell
    integer(iintegers) :: left_face   ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: right_face  ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: upper_face  ! index of face which is top/bot of base face, index regarding faces_of_cell
    integer(iintegers) :: bottom_face ! index of face which is top/bot of base face, index regarding faces_of_cell

    real(ireals),target :: dir2dir(5,5) ! Nface**2, dim=(ncol, nrow), i.e. dim=(isrc, idst), isrc changing quickest
    real(ireals),pointer :: pdir2dir(:)
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: zindex
    real(ireals) :: dx, dz, coeff(5**2) ! coefficients for each src=[1..5] and dst[1..5]

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, sec, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(plex%edir_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(plex%edir_dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges

    call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)
    !call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)

    pdir2dir(1:size(dir2dir)) => dir2dir

    do icell = cStart, cEnd-1

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      call compute_local_wedge_ordering(plex, icell, faces_of_cell, geomSection, geoms, sundir, &
        zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)

      zindex = plex%zindex(icell)
      dz = plex%hhl(zindex) - plex%hhl(zindex+1)

      call DMPlexGetCone(plex%edir_dm, faces_of_cell(upper_face), edges_of_face, ierr); call CHKERR(ierr) ! Get Faces of cell
      dx = 0
      do i=1,size(edges_of_face)
        iedge = edges_of_face(i)
        call PetscSectionGetOffset(geomSection, iedge, geom_offset, ierr); call CHKERR(ierr)
        dx = dx + geoms(i1+geom_offset)
      enddo
      dx = dx / size(edges_of_face)
      call DMPlexRestoreCone(plex%edir_dm, faces_of_cell(upper_face), edges_of_face, ierr); call CHKERR(ierr)

      call get_coeff(OPP, optprop(icell+1)%kabs, optprop(icell+1)%ksca, optprop(icell+1)%g, &
        dz, dx, .True., coeff, angles=[rad2deg(azimuth), rad2deg(zenith)])

      do iface = 1, size(faces_of_cell)

        ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering
        if (iface.eq.upper_face) then
          iwedgeface = 1
        else if (iface.eq.bottom_face) then
          iwedgeface = 5
        else if (iface.eq.base_face) then
          iwedgeface = 2
        else if (iface.eq.left_face) then
          iwedgeface = 3
        else if (iface.eq.right_face) then
          iwedgeface = 4
        else
          stop 'iface is not in local wedge face numbering... something must have gone terribly wrong in compute_local_wedge_ordering'
        endif

        if(.not.lsrc(iface)) then
          dir2dir([upper_face,base_face,left_face,right_face,bottom_face], iface) = -coeff((iwedgeface-1)*5+1:iwedgeface*5)
        else ! This is a src face, no radiation gets here unless we have a boundary condition to apply
          dir2dir(:, iface) = zero
        endif !lsrc
        dir2dir(iface, iface) = one ! setting diagonal entry

        call PetscSectionGetOffset(sec, faces_of_cell(iface), irow, ierr); call CHKERR(ierr)
        irows(iface) = irow
        icols(iface) = irow
      enddo ! enddo iface
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)

      !print *,'Setting Direct Matrix values rows/cols:', irows, '::', icols
      !print *,'face order:', [upper_face,base_face,left_face,right_face,bottom_face]
      !print *,'Setting Direct Matrix values coeffs 1 :', dir2dir(:,1)
      !print *,'Setting Direct Matrix values coeffs 2 :', dir2dir(:,2)
      !print *,'Setting Direct Matrix values coeffs 3 :', dir2dir(:,3)
      !print *,'Setting Direct Matrix values coeffs 4 :', dir2dir(:,4)
      !print *,'Setting Direct Matrix values coeffs 5 :', dir2dir(:,5)

      !if(ldebug) then
      !  do iface = 1, size(faces_of_cell)
      !    print *,'Sum for iface', iface,'::', sum(dir2dir(iface,:))
      !  enddo
      !endif

      call MatSetValuesLocal(A, i5, irows, i5, icols, pdir2dir, INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Medir', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...finished'
    !stop 'debug'
  end subroutine

  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(OPP, kabs, ksca, g, dz, dx, ldir, coeff, angles)
    class(t_optprop), intent(in) :: OPP
    real(ireals), intent(in)     :: kabs, ksca, g
    real(ireals),intent(in)      :: dz, dx
    logical,intent(in)           :: ldir
    real(ireals),intent(out)     :: coeff(:)

    real(ireals),intent(in),optional  :: angles(2)

    real(ireals) :: aspect, tauz, w0

    aspect = dz / dx
    tauz = (kabs+ksca) * dz
    if(approx(tauz,zero)) then
      w0 = zero
    else
      w0 = ksca / (kabs+ksca)
    endif

    !print *,'Lookup Coeffs for', aspect, tauz, w0, g,'::',angles
    call OPP%get_coeff(aspect, tauz, w0, g,ldir, coeff, inp_angles=angles)
    !print *,'Coeffs', coeff
  end subroutine

  subroutine compute_local_wedge_ordering(plex, icell, faces_of_cell, geomSection, geoms, sundir, &
      zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)
    type(t_plexgrid), intent(in) :: plex
    integer(iintegers), intent(in) :: icell
    integer(iintegers), intent(in), pointer :: faces_of_cell(:)
    type(tPetscSection) :: geomSection
    real(ireals), intent(in), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), intent(in) :: sundir(3)
    real(ireals), intent(out) :: zenith, azimuth
    integer(iintegers), intent(out) :: upper_face, bottom_face, base_face, left_face, right_face
    logical, intent(out) :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: iface, geom_offset, izindex(2)

    real(ireals) :: cell_center(3)
    real(ireals) :: face_normals(3,5), face_centers(3,5)
    real(ireals) :: side_faces_angles_to_sun(5), proj_angles_to_sun(3), proj_normal(3)
    real(ireals) :: e_x(3), e_y(3), e_z(3) ! unit vectors of local coord system in which we compute the transfer coefficients
    real(ireals) :: projected_sundir(3)

    integer(iintegers) :: side_faces(3), top_faces(2) ! indices in faces_of_cell which give the top/bot and side faces via labeling
    integer(iintegers) :: iside_faces, itop_faces, ibase_face ! indices to fill above arrays
    real(ireals) :: MrotWorld2Local(3,3) ! Rotation Matrix into the local wedge space, (ex, ey, ez)

    real(ireals) :: local_normal_left(3), local_normal_right(3) ! normal vectors in local wedgemc geometry. For even sided triangles, this is smth. like left: [.5, -.8] or right [-.5, -.8]

    integer(mpiint) :: ierr

    call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
    cell_center = geoms(1+geom_offset:3+geom_offset)

    do iface = 1, size(faces_of_cell)
      call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); call CHKERR(ierr)
      face_centers(:,iface) = geoms(1+geom_offset: 3+geom_offset)
      face_normals(:,iface) = geoms(4+geom_offset: 6+geom_offset)
      face_normals(:,iface) = face_normals(:,iface) * &
        real(determine_normal_direction(face_normals(:,iface), face_centers(:, iface), cell_center), ireals)
      ! to find the face whose normal is closest to sun vector:
      side_faces_angles_to_sun(iface) = angle_between_two_vec(face_normals(:,iface), sundir) ! in [rad]
      lsrc(iface) = side_faces_angles_to_sun(iface).le.pi/2-deg2rad(.01_ireals) ! dont propagate energy along edges where sun is only .1 degrees off
    enddo

    ! get numbering for 2 top/bot faces and 3 side faces which indicate the position in the faces_of_cell vec
    iside_faces = 1; itop_faces = 1
    do iface = 1, size(faces_of_cell)
      if(plex%ltopfacepos(faces_of_cell(iface))) then
        top_faces(itop_faces) = iface
        itop_faces = itop_faces+1
      else
        side_faces(iside_faces) = iface
        iside_faces = iside_faces+1
      endif
    enddo

    izindex(1) = plex%zindex(faces_of_cell(top_faces(1)))
    izindex(2) = plex%zindex(faces_of_cell(top_faces(2)))

    if(izindex(1).gt.izindex(2)) then
      upper_face = top_faces(2)
      bottom_face = top_faces(1)
    else
      upper_face = top_faces(1)
      bottom_face = top_faces(2)
    endif

    do iface=1,size(side_faces)
      if(all(approx(sundir, face_normals(:,upper_face)))) then
        proj_angles_to_sun(iface) = 0
      else
        proj_normal = vec_proj_on_plane(sundir, face_normals(:,upper_face))
        !print *,'vec_proj_on_plane',iface, '::', sundir, 'vs', face_normals(:,upper_face), '->', proj_normal
        proj_angles_to_sun(iface) = angle_between_two_vec(proj_normal, face_normals(:,side_faces(iface)))
      endif
    enddo

    ibase_face = minloc(proj_angles_to_sun,dim=1)
    base_face = side_faces(ibase_face)

    e_y = face_normals(:, base_face)   ! inward facing normal -> in local wedgemc coordinates
    e_z = -face_normals(:, upper_face) ! outward facing normal with respect to the top plate
    e_x = cross_3d(e_y, e_z)           ! in local wedge_coords, this is y=0 coordinate

    MrotWorld2Local = rotation_matrix_world_to_local_basis(e_x, e_y, e_z)

    left_face  = side_faces(modulo(ibase_face,size(side_faces, kind=iintegers))+i1)
    right_face = side_faces(modulo(ibase_face+i1,size(side_faces, kind=iintegers))+i1)

    local_normal_left  = matmul(MrotWorld2Local, face_normals(:,left_face))
    local_normal_right = matmul(MrotWorld2Local, face_normals(:,right_face))

    if(local_normal_left(1).lt.local_normal_right(1)) then ! switch right and left face
      iface = right_face
      right_face = left_face
      left_face = iface
    endif

    ! Now we all the info for the local wedge calculations as we do em with the MonteCarlo raytracer

    zenith = angle_between_two_vec(sundir, -e_z)

    projected_sundir = vec_proj_on_plane(matmul(MrotWorld2Local, sundir), [zero,zero,one])
    if(norm(projected_sundir).le.epsilon(zero)) then
      azimuth = 0
    else
      projected_sundir = projected_sundir
      azimuth = angle_between_two_vec([zero,one,zero], projected_sundir) * sign(one, projected_sundir(1))
    endif

    if(ldebug .and. norm(face_centers(:,upper_face)) .le. norm(face_centers(:,bottom_face))) then ! we expect the first face to be the upper one
      print *,'norm upper_face ', norm(face_centers(:,upper_face))
      print *,'norm bottom_face', norm(face_centers(:,bottom_face))
      print *,'we expect the first face to be the upper one but found:',icell, faces_of_cell(1), faces_of_cell(2)
      stop 'create_edir_mat() :: wrong zindexlabel'
    endif

    if(azimuth.lt.-60 .or. azimuth.gt.60) stop 'local azimuth greater than 60 deg. something must have gone wrong with the base face selection!'
    end subroutine

    !subroutine compute_dir2dir_coeff(src, phi, theta, S,T)
    !  use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_5
    !  integer(iintegers), intent(in) :: src
    !  real(ireals), intent(in) :: phi, theta
    !  real(ireals), intent(out) :: S(5),T(5)

    !  type(t_boxmc_wedge_5_5) :: bmc_wedge_5_5
    !  real(ireals) :: bg(3), dx,dy,dz
    !  real(ireals) :: S_tol(5),T_tol(5)

    !  call bmc_wedge_5_5%init(PETSC_COMM_SELF)
    !  print *,'computing coeffs for src/phi/theta',src,phi,theta

    !  bg  = [1e-3_ireals, zero, one/2 ]

    !  !phi   =  0
    !  !theta = 45

    !  dx = 100
    !  dy = dx
    !  dz = 200

    !  call bmc_wedge_5_5%get_coeff(PETSC_COMM_SELF, bg, src, .True., &
    !    phi, theta, dx, dy, dz, S, T, S_tol, T_tol, inp_atol=5e-2_ireals, inp_rtol=1e-1_ireals)
    !end subroutine

    function get_normal_of_first_TOA_face(plex)
      type(t_plexgrid) :: plex
      real(ireals) :: get_normal_of_first_TOA_face(3)
      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      real(ireals) :: cell_center(3), face_center(3)
      real(ireals),allocatable :: face_normal(:)

      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xitoa(:), cell_support(:)
      integer(iintegers) :: geom_offset, iface, icell

      integer(mpiint) :: myid, ierr

      call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        if(.not.allocated(plex%geom_dm)) stop 'get_normal_of_first_TOA_face::needs allocated geom_dm first'
        call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
        !call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        !call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        call DMGetStratumIS(plex%geom_dm, 'TOA', i1, toa_ids, ierr); call CHKERR(ierr)

        if (toa_ids.eq.PETSC_NULL_IS) then ! dont have TOA points
          stop 'This didnt work, we tried to set the sundir according to first face on rank 0 but it seems he does not have TOA faces'
        else
          call ISGetIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
          iface = xitoa(1) ! first face of TOA faces

          call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
          call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)

          face_center = geoms(geom_offset+1:geom_offset+3)
          allocate(face_normal(3))
          face_normal = geoms(geom_offset+4:geom_offset+6)

          call DMPlexGetSupport(plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

          call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
          cell_center = geoms(1+geom_offset:3+geom_offset)

          ! Determine the inward normal vec for the face
          face_normal = face_normal * real(determine_normal_direction(face_normal, face_center, cell_center), kind=ireals)

          call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

          call ISRestoreIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
        endif
      endif

      call imp_bcast(plex%comm, face_normal, 0_mpiint)
      get_normal_of_first_TOA_face = face_normal

    end function
end module
