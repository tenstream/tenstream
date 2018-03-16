module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options

  use m_helper_functions, only: CHKERR, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, strF2C, get_arg, &
    vec_proj_on_plane, cross_3d, norm, rotation_matrix_world_to_local_basis, &
    imp_bcast, approx, swap

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, i6, i7, i8, &
    zero, one, pi, EXP_MINVAL, EXP_MAXVAL

  use m_plex_grid, only: t_plexgrid, compute_face_geometry, setup_edir_dmplex, setup_abso_dmplex, &
    orient_face_normals_along_sundir, compute_wedge_orientation, is_solar_src, get_inward_face_normal, &
    facevec2cellvec, icell_icon_2_plex, iface_top_icon_2_plex, &
    TOAFACE, BOTFACE, SIDEFACE

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


  type t_state_container
    integer(iintegers)  :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    type(tVec), allocatable :: edir, ediff, abso

    logical             :: lset        = .False. ! initialized?
    logical             :: lsolar_rad  = .False. ! direct radiation calculated?
    logical             :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](true) or [W/m**2](false)
    logical             :: lWm2_dir=.False. , lWm2_diff=.False.
  end type

  type t_plex_solver
    type(t_plexgrid), allocatable :: plex
    class(t_optprop), allocatable :: OPP

    type(t_opticalprops), allocatable, dimension(:) :: optprop ! in each cell [pStart..pEnd-1]
    !real(ireals)        , allocatable, dimension(:) :: planck  ! in each cell [pStart+1..pEnd]
    !real(ireals)        , allocatable, dimension(:) :: albedo ! on each surface face [?TODO]

    type(t_state_container) :: solutions(-1000:1000)

    type(tVec), allocatable :: incSolar
    type(tMat) :: Mdir
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
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], [real(OPP_LUT_ALL_ANGLES, kind=ireals)], solver%plex%comm)
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], real([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70],ireals), plex%comm)
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], real([0, 5, 10, 15],ireals), plex%comm)
      !call solver%OPP%init([real(OPP_LUT_ALL_ANGLES, kind=ireals)], real([0],ireals), plex%comm)

      call solver%OPP%init(plex%comm)
    end subroutine

    subroutine set_plex_rt_optprop(solver, vlwc, viwc)
      use m_helper_functions, only : delta_scale
      type(t_plex_solver), allocatable, intent(inout) :: solver
      type(tVec),intent(in), optional :: vlwc, viwc
      real(ireals), pointer :: xlwc(:), xiwc(:)

      real(ireals),parameter :: w0 = .8, reff_w=10, reff_i=80, rayleigh=1e-8
      real(ireals) :: ksca_tot, kabs_tot, g_tot
      real(ireals) :: ksca_cld, kabs_cld, g_cld

      integer(iintegers) :: i
      integer(mpiint) :: myid, ierr

      if(.not.allocated(solver)) call CHKERR(1_mpiint, 'set_plex_rt_optprop::solver has to be allocated')
      if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'set_plex_rt_optprop::solver%plex has to be allocated')
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(solver%optprop)) allocate(solver%optprop(solver%plex%cStart:solver%plex%cEnd-1))
      associate( optprop => solver%optprop )

      if(present(vlwc)) call VecGetArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      if(present(viwc)) call VecGetArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)

      do i = 1,size(optprop)
        kabs_tot = rayleigh*(one-w0)
        ksca_tot = rayleigh*w0
        g_tot    = zero

        if(present(vlwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xlwc(i) * 1e3) / reff_w * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xlwc(i) * 1e3) / reff_w * (w0)
          g_cld    = .5_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = g_cld * ksca_cld / (ksca_cld + ksca_tot)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        if(present(viwc)) then
          kabs_cld = 3._ireals / 2._ireals * (xiwc(i) * 1e3) / reff_i * (one-w0)
          ksca_cld = 3._ireals / 2._ireals * (xiwc(i) * 1e3) / reff_i * (w0)
          g_cld    = .3_ireals
          call delta_scale( kabs_cld, ksca_cld, g_cld )
          g_tot    = (g_tot * ksca_tot + g_cld * ksca_cld) / (ksca_tot + ksca_cld)
          kabs_tot = kabs_tot + kabs_cld
          ksca_tot = ksca_tot + ksca_cld
        endif

        optprop(i-1)%kabs = kabs_tot
        optprop(i-1)%ksca = ksca_tot
        optprop(i-1)%g    = g_tot
      enddo

      if(present(vlwc)) call VecRestoreArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
      if(present(viwc)) call VecRestoreArrayReadF90(viwc, xiwc, ierr); call CHKERR(ierr)

      print *,'Min/Max of kabs', minval(optprop(:)%kabs), maxval(optprop(:)%kabs)
      print *,'Min/Max of ksca', minval(optprop(:)%ksca), maxval(optprop(:)%ksca)
      print *,'Min/Max of g   ', minval(optprop(:)%g   ), maxval(optprop(:)%g   )

      end associate
    end subroutine

    subroutine run_plex_rt_solver(solver, sundir, opt_solutions_uid)
      type(t_plex_solver), allocatable, intent(inout) :: solver
      real(ireals), intent(in) :: sundir(3) ! cartesian direction of sun rays, norm of vector is the energy in W/m2
      integer(iintegers), intent(in), optional :: opt_solutions_uid

      integer(iintegers) :: suid

      integer(mpiint) :: myid, ierr

      if(.not.allocated(solver)) call CHKERR(1_mpiint, 'run_plex_rt_solver::solver has to be allocated')

      if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'run_plex_rt_solver::plex has to be allocated first')
      call mpi_comm_rank(solver%plex%comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(solver%plex%geom_dm)) call CHKERR(myid+i1, 'run_plex_rt_solver::geom_dm has to be allocated first')
      if(.not.allocated(solver%optprop)) call CHKERR(myid+i1, 'run_plex_rt_solver::optprop has to be allocated first')
      if(.not.allocated(solver%plex%geom_dm)) call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      if(.not.allocated(solver%plex%edir_dm)) call setup_edir_dmplex(solver%plex, solver%plex%edir_dm)
      if(.not.allocated(solver%plex%abso_dm)) call setup_abso_dmplex(solver%plex, solver%plex%abso_dm)

      ! Prepare the space for the solution
      suid = get_arg(i0, opt_solutions_uid)
      associate( solution => solver%solutions(suid) )

      solution%lsolar_rad = norm(sundir).gt.zero

      if(solution%lsolar_rad) then
        call orient_face_normals_along_sundir(solver%plex, sundir)
        call compute_wedge_orientation(solver%plex, sundir, solver%plex%wedge_orientation_dm, solver%plex%wedge_orientation)
        ! Output of wedge_orient vec
        call create_edir_src_vec(solver%plex, solver%plex%edir_dm, norm(sundir), solver%optprop, sundir/norm(sundir), solver%incSolar)

        ! Output of srcVec
        if(ldebug) then
          call scale_facevec(solver%plex, solver%plex%edir_dm, solver%incSolar, lW_to_Wm2=.True.)
          call facevec2cellvec(solver%plex, solver%plex%edir_dm, solver%incSolar)
          call scale_facevec(solver%plex, solver%plex%edir_dm, solver%incSolar, lW_to_Wm2=.False.)
        endif

        !return!DEBUG
        ! Create Direct Matrix
        call create_edir_mat(solver%plex, solver%OPP, solver%optprop, solver%Mdir)

        ! Solve Direct Matrix
        call solve_plex_rt(solver%plex, solver%incSolar, solver%Mdir, solution%edir)
        call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-show_edir_vec_global', ierr); call CHKERR(ierr)
        solution%lWm2_dir = .False.

        ! Output of Edir
        if(ldebug) then
          call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.True.)
          call facevec2cellvec(solver%plex, solver%plex%edir_dm, solution%edir)
          call scale_facevec(solver%plex, solver%plex%edir_dm, solution%edir, lW_to_Wm2=.False.)
        endif
      endif

      call compute_edir_absorption(solver%plex, solution%edir, solution%abso)

      end associate
    end subroutine

    !> @brief setup source term for diffuse radiation
    !> @details this is either direct radiation scattered into one of the diffuse coeffs:
    !> \n direct source term is
    !> \n   direct radiation times the dir2diff coeffs
    !> \n or it may be that we have a source term due to thermal emission --
    !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
    !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
    !subroutine setup_b(plex, OPP, optprop, solution, srcVec)
    !  type(t_plexgrid), intent(in) :: plex
    !  class(t_optprop), intent(in) :: OPP
    !  type(t_opticalprops), intent(in) :: optprop(:)
    !  type(t_state_container), intent(in) :: solution
    !  type(tVec), allocatable, intent(inout) :: srcVec
    !end subroutine

    subroutine create_edir_src_vec(plex, edirdm, E0, optprop, sundir, srcVec)
      type(t_plexgrid), allocatable, intent(in) :: plex
      type(tDM),allocatable, intent(in) :: edirdm
      real(ireals), intent(in) :: E0, sundir(3)
      type(t_opticalprops), allocatable, intent(in) :: optprop(:)
      type(tVec), allocatable, intent(inout) :: srcVec

      integer(iintegers) :: i, voff
      type(tPetscSection) :: s

      type(tVec) :: localVec, lambertVec
      real(ireals), pointer :: xv(:), xlambert(:)

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

      if(.not.allocated(edirdm)) call CHKERR(myid+i1, 'called create_src_vec but edirdm is not allocated')

      if(.not.allocated(plex%geom_dm)) call CHKERR(myid+i1, 'get_normal_of_first_TOA_face::needs allocated geom_dm first')

      call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetDefaultSection(edirdm, s, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(s, PETSC_NULL_SECTION, '-show_src_section', ierr); call CHKERR(ierr)

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
          call compute_lambert_beer(plex, icell, E0, sundir, optprop, lambertVec)
        endif
      enddo

      call VecGetArrayReadF90(lambertVec, xlambert, ierr); call CHKERR(ierr)
      do iface = plex%fStart, plex%fEnd-1
        call DMLabelGetValue(plex%domainboundarylabel, iface, labelval, ierr); call CHKERR(ierr)
        if(labelval.eq.SIDEFACE)then
          call DMPlexGetSupport(edirdm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(edirdm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

          call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)

          if(is_solar_src(face_normal, sundir)) then
            call PetscSectionGetOffset(s, iface, voff, ierr); call CHKERR(ierr)
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

      if(ldebug .and. myid.eq.0) print *,myid,'plex_rt::create_src_vec....finished'
      if(ldebug) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    end subroutine

    subroutine compute_lambert_beer(plex, itopcell, E0, sundir, optprop, lambertVec)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: itopcell
      real(ireals), intent(in) :: E0, sundir(3) ! E0 in W/m2
      type(t_opticalprops), allocatable, intent(in) :: optprop(:)
      type(tVec),intent(inout) :: lambertVec

      type(tPetscSection) :: s
      real(ireals), pointer :: xv(:) ! vertical entries are the direct radiation flux through the tilted plane
      integer(iintegers), pointer :: faces_of_cell(:)
      integer(iintegers) :: i, iconcell, icell, iface_side, iface_top, iface_bot, k, voff

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: geom_offset
      real(ireals) :: face_normal(3), face_normal_top(3)

      real(ireals) :: dz, kext, dtau, mu, mu_top, mu_side
      real(ireals) :: area, expon, sin_alpha, transport

      integer(mpiint) :: ierr

      call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call DMGetDefaultSection(plex%edir_dm, s, ierr); call CHKERR(ierr)
      call VecGetArrayF90(lambertVec, xv, ierr); call CHKERR(ierr)

      iconcell = plex%localiconindex(itopcell)
      k = plex%zindex(itopcell)

      iface_top = iface_top_icon_2_plex(plex, iconcell, k)

      call get_inward_face_normal(iface_top, itopcell, geomSection, geoms, face_normal_top)

      if(is_solar_src(face_normal_top, sundir)) then ! if sun actually shines on top of the column
        mu_top = dot_product(sundir, face_normal_top)

        dtau = zero

        !call PetscSectionGetOffset(s, iface_top, voff, ierr); call CHKERR(ierr)
        !xv(voff+1) = E0 * mu_top

        do k = 1, plex%Nz
          iface_top = iface_top_icon_2_plex(plex, iconcell, k)
          iface_bot = iface_top_icon_2_plex(plex, iconcell, k+1)

          icell = icell_icon_2_plex(plex, iconcell, k)

          call get_inward_face_normal(iface_top, icell, geomSection, geoms, face_normal_top)
          call get_inward_face_normal(iface_bot, icell, geomSection, geoms, face_normal)

          if(.not.is_solar_src(-face_normal, sundir)) cycle ! if the bot face is not sunlit... we dont allow the sun to shine from below

          call PetscSectionGetOffset(s, iface_top, voff, ierr); call CHKERR(ierr)

          !mu_top = dot_product(sundir, face_normal_top)
          dz = plex%hhl(k) - plex%hhl(k+1)
          kext = (optprop(icell)%kabs + optprop(icell)%ksca)
          !print *,'k',k,'tau increment', icell, optprop(icell)%kabs , optprop(icell)%ksca , dz , ':', dtau, '->', dtau + kext * dz

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
                !print *,'solar_src', iface_side, E0, dtau, transport, mu_side, area, '->', E0 * transport * mu_side, xv(voff+i1)
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

    subroutine solve_plex_rt(plex, b, A, x)
      type(t_plexgrid) :: plex
      type(tVec), intent(in) :: b
      type(tVec), allocatable, intent(inout) :: x
      type(tMat), intent(in) :: A

      type(tKSP) :: ksp
      integer(mpiint) :: ierr

      if(ldebug) print *,'plex_rt::solve Matrix...'
      call KSPCreate(plex%comm, ksp, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, ierr); CHKERRQ(ierr)

      call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
      call KSPSetUp(ksp, ierr); CHKERRQ(ierr)

      if(.not.allocated(x)) then
        allocate(x)
        call VecDuplicate(b, x, ierr); CHKERRQ(ierr)
      endif
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

    integer(mpiint) :: myid, ierr

    if(ldebug) print *,'plex_rt::scale_facevec...'
    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(face_dm, faceSection, ierr); call CHKERR(ierr)

    call DMGetLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)

    call DMGlobalToLocalBegin(face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (face_dm, globalfaceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

    call VecGetArrayF90(faceVec, xv, ierr); call CHKERR(ierr)

    if(ldebug .and. myid.ne.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    do iface = plex%fStart, plex%fEnd-1
      call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
      area = geoms(geom_offset+i7)

      call PetscSectionGetOffset(faceSection, iface, face_offset, ierr); call CHKERR(ierr)
      if(lW_to_Wm2) then
        xv(face_offset+i1) = xv(face_offset+i1) / area
      else
        xv(face_offset+i1) = xv(face_offset+i1) * area
      endif
    enddo
    if(ldebug .and. myid.eq.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)

    call VecRestoreArrayF90(faceVec, xv, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(face_dm, faceVec, INSERT_VALUES, globalfaceVec, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (face_dm, faceVec, INSERT_VALUES, globalfaceVec, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(face_dm, faceVec, ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::scale_facevec... finished'
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

    if(.not.allocated(plex%edir_dm)) call CHKERR(myid+i1, 'called compute_edir_absorption with a dm which is not allocated?')
    if(.not.allocated(plex%abso_dm)) call CHKERR(myid+i1, 'called compute_edir_absorption with a dm which is not allocated?')

    if(ldebug) print *,'plex_rt::compute_edir_absorption....'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetDefaultSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!

    call DMGetLocalVector(plex%abso_dm, local_abso,ierr); call CHKERR(ierr)
    call VecSet(local_abso, zero, ierr); call CHKERR(ierr)

    call DMGetLocalVector(plex%edir_dm, local_edir,ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    !if(ldebug .and. myid.ne.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    !do iface = plex%fStart, plex%fEnd-1
    !  call PetscSectionGetOffset(edir_section, iface, edir_offset, ierr); call CHKERR(ierr)
    !  print *,myid,'edir for ', iface, '::', edir_offset, xedir(edir_offset+i1)
    !enddo
    !if(ldebug .and. myid.eq.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)

    if(ldebug .and. myid.eq.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
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
          !print *,myid,'absorption for ', icell, iface, faces_of_cell(iface), '+',&
          !  xedir(edir_offset+i1),'::',xabso(abso_offset+i1) / volume, '->', (xabso(abso_offset+i1) + xedir(edir_offset+i1))/volume
          xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1)
        else
          !print *,myid,'absorption for ', icell, iface, faces_of_cell(iface), '-',&
          !  xedir(edir_offset+i1),'::',xabso(abso_offset+i1) / volume, '->', (xabso(abso_offset+i1) - xedir(edir_offset+i1))/volume
          xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1)
        endif
      enddo
      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      xabso(abso_offset+i1) = xabso(abso_offset+i1) / volume
    enddo
    if(ldebug .and. myid.ne.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)

    !if(ldebug .and. myid.eq.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    !do icell = cStart, cEnd-1
    !  call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)
    !  print *,myid,'absorption for ', icell, '::', xabso(abso_offset+i1)
    !enddo
    !if(ldebug .and. myid.ne.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(local_abso, xabso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%abso_dm, local_edir, ierr); call CHKERR(ierr)

    if(.not.allocated(abso)) then
      allocate(abso)
      call DMCreateGlobalVector(plex%abso_dm, abso, ierr); call CHKERR(ierr)
      call PetscObjectSetName(abso, 'absoVecGlobal', ierr);call CHKERR(ierr)
    endif
    call VecSet(abso, zero, ierr); call CHKERR(ierr)
    call DMLocalToGlobalBegin(plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd  (plex%abso_dm, local_abso, ADD_VALUES, abso, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(plex%abso_dm, local_abso, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
    if(ldebug) print *,'plex_rt::compute_edir_absorption....finished'
  end subroutine

  subroutine create_edir_mat(plex, OPP, optprop, A)
    type(t_plexgrid), intent(inout) :: plex
    class(t_optprop), intent(in) :: OPP
    type(t_opticalprops), allocatable, intent(in) :: optprop(:)
    type(tMat), intent(out) :: A

    type(tPetscSection) :: sec

    integer(mpiint) :: myid, ierr

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: icell, iface, irow, icol, isrc, idst

    type(tPetscSection) :: geomSection, wedgeSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    real(ireals), pointer :: wedgeorient(:) ! pointer to orientation vec
    integer(iintegers) :: wedge_offset

    ! iwedge_plex2bmc :: mapping from plex_indices, i.e. iface(1..5) to boxmc wedge numbers,
    ! i.e. iwedge_plex2bmc(1) gives the wedge boxmc position of first dmplex face
    integer(iintegers) :: iwedge_plex2bmc(5)

    ! iwedge_bmc2plex :: mapping from boxmc wedge numbers to plex indices,
    ! i.e. iwedge_bmc2plex(1) gives the plex index for top face
    integer(iintegers) :: iwedge_bmc2plex(5)

    real(ireals) :: zenith, azimuth

    real(ireals) :: dir2dir(5), T(5), S(8)
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: zindex
    real(ireals) :: dx, dz, coeff(5**2) ! coefficients for each src=[1..5] and dst[1..5]

    logical, parameter :: lonline=.True.

    call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, sec, ierr); call CHKERR(ierr)

    call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)
    !call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)

    !if(ldebug .and. myid.ne.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
    do icell = plex%cStart, plex%cEnd-1
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      !call compute_local_wedge_ordering(plex, icell, geomSection, geoms, sundir, &
      !  lsrc, zenith, azimuth, dx, &
      !  upper_face, bottom_face, base_face, left_face, right_face)

      call PetscSectionGetOffset(wedgeSection, icell, wedge_offset, ierr); call CHKERR(ierr)
      zenith  = wedgeorient(wedge_offset+i1)
      azimuth = wedgeorient(wedge_offset+i2)
      dx      = wedgeorient(wedge_offset+i3)

      do iface = 1, size(faces_of_cell)
        iwedge_plex2bmc(int(wedgeorient(wedge_offset+i3+iface), iintegers)) = iface
        lsrc(iface) = int(wedgeorient(wedge_offset+i8+iface), iintegers) .gt. one/i2
      enddo
      !print *,'iwedge_plex2bmc', iwedge_plex2bmc

      zindex = plex%zindex(icell)
      dz = plex%hhl(zindex) - plex%hhl(zindex+1)

      !print *,'icell',icell,': foc',faces_of_cell
      !print *,'icell',icell,':',iwedge_plex2bmc,'lsrc',lsrc
      call get_coeff(OPP, optprop(icell)%kabs, optprop(icell)%ksca, optprop(icell)%g, &
        dz, dx, .True., coeff, angles=[rad2deg(azimuth), rad2deg(zenith)])

      !do isrc = 1, 5
      !  print *,'LUT',rad2deg(azimuth), rad2deg(zenith), isrc, coeff(isrc:size(coeff):i5)
      !enddo
      do iface = 1, size(faces_of_cell)

        !if(lsrc(iface)) then
          ! we have to reorder the coefficients to their correct position from the local LUT numbering into the petsc face numbering

          !dir2dir([upper_face, base_face, left_face, right_face, bottom_face]) = coeff((iwedgeface-1)*i5+i1:iwedgeface*i5)
          dir2dir = coeff(iwedge_plex2bmc(iface):size(coeff):i5)
          if(lonline) then
            !print *,'lonline LUT0', coeff(iwedge_plex2bmc(iface):size(coeff):i5)
            !call compute_dir2dir_coeff(iwedge_plex2bmc(iface), dz, dx, &
            !  optprop(icell)%kabs, optprop(icell)%ksca, optprop(icell)%g, &
            !  rad2deg(azimuth), rad2deg(zenith), S, T)
            !print *,'lonline BOX1', T
            call compute_dir2dir_coeff(iwedge_plex2bmc(iface), dz, dx, &
              optprop(icell)%kabs, optprop(icell)%ksca, optprop(icell)%g, &
              rad2deg(azimuth), rad2deg(zenith), S, T, wedgeorient(wedge_offset+14:wedge_offset+19))
            !print *,'lonline BOX2', T
            dir2dir = T
          endif

          call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr)

          do idst = 1, size(faces_of_cell)
            !if(.not.lsrc(idst)) then
              call PetscSectionGetOffset(sec, faces_of_cell(idst), irow, ierr); call CHKERR(ierr)
              !print *,'isrc', iwedge_plex2bmc(iface), 'idst', iwedge_plex2bmc(idst), &
              !  'if', iface, 'id', idst, &
              !  'srcface->dstface', faces_of_cell(iface),faces_of_cell(idst), &
              !  'col -> row', icol, irow, dir2dir(iwedge_plex2bmc(idst))
              call MatSetValuesLocal(A, i1, [irow], i1, [icol], [-dir2dir(iwedge_plex2bmc(idst))], ADD_VALUES, ierr); call CHKERR(ierr)
            !endif
          enddo

        !else ! This is not a src face, no radiation comes from here
        !endif

      enddo ! enddo iface

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo
    !if(ldebug .and. myid.eq.0) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)

    call MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)

    ! Set Diagonal Entries
    do iface = plex%fStart, plex%fEnd-1
      call PetscSectionGetOffset(sec, iface, irow, ierr); call CHKERR(ierr)
      call MatSetValuesLocal(A, i1, [irow], i1, [irow], [one], INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, '-show_Medir', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%wedge_orientation, wedgeorient, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,'plex_rt::create_edir_mat...finished'
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

    call OPP%get_coeff(aspect, tauz, w0, g, ldir, coeff, inp_angles=angles)
    !print *,'DEBUG Lookup Coeffs for', aspect, tauz, w0, g, angles,'::', coeff
    if(ldebug) then
      if(any(coeff.lt.zero).or.any(coeff.gt.one)) then
        print *,'Lookup Coeffs for', aspect, tauz, w0, g, angles,'::', coeff
        call CHKERR(i1, 'Found corrupted coefficients!')
      endif
    endif

  end subroutine

  subroutine compute_dir2dir_coeff(src, dz, dx, kabs, ksca, g, phi, theta, S, T, coord2d)
    use m_boxmc, only : t_boxmc, t_boxmc_wedge_5_5, t_boxmc_wedge_5_8
    integer(iintegers), intent(in) :: src
    real(ireals), intent(in) :: dz, dx, kabs, ksca, g, phi, theta
    real(ireals), intent(out) :: S(:),T(:)
    real(ireals), intent(in), optional :: coord2d(:)

    type(t_boxmc_wedge_5_8) :: bmc_wedge
    real(ireals) :: bg(3), dy
    real(ireals) :: S_tol(size(S)),T_tol(size(T))

    call bmc_wedge%init(PETSC_COMM_SELF, coord2d)
    if(present(coord2d)) then
      print *,'computing coeffs for src/phi/theta',src,phi,theta,':',coord2d
    else
      print *,'computing coeffs for src/phi/theta',src,phi,theta
    endif

    bg  = [kabs, ksca, g]

    dy = dx

    call bmc_wedge%get_coeff(PETSC_COMM_SELF, bg, src, .True., &
      phi, theta, dx, dy, dz, S, T, S_tol, T_tol, inp_atol=1e-3_ireals, inp_rtol=5e-2_ireals)
  end subroutine

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
        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)

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
