module m_plexrt_external_solvers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    i0, i1, i2, i3, zero
  use m_helper_functions, only: CHKERR, angle_between_two_vec, delta_scale

  use m_pprts_base, only : t_state_container
  use m_plex_rt_base, only: t_plex_solver

  use m_plex_grid, only: t_plexgrid, &
    get_inward_face_normal, &
    get_consecutive_vertical_cell_idx, &
    TOAFACE

  use m_schwarzschild, only: schwarzschild, B_eff
  use m_twostream, only: delta_eddington_twostream
  use m_tenstr_disort, only: default_flx_computation

  implicit none

contains
  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine plexrt_schwarz(solver, solution)
    class(t_plex_solver)    :: solver
    type(t_state_container) :: solution

    real(ireals),allocatable :: dtau(:), Blev(:), Edn(:),Eup(:)

    type(tIS) :: boundary_ids
    integer(iintegers), pointer :: xitoa(:), cell_support(:)
    integer(iintegers), allocatable :: cell_idx(:)
    integer(iintegers) :: i, k, icell, iface, voff, ke1, geom_offset
    real(ireals) :: dz
    real(ireals), pointer :: xkabs(:), xalbedo(:), xplck(:), xediff(:), xabso(:), xgeoms(:)
    type(tPetscSection) :: ediff_section, abso_section, plck_section, geom_section

    integer(iintegers) :: Nmu
    logical :: lflg
    integer(mpiint) :: ierr

    if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
    if( .not. allocated(solver%plck) ) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

    Nmu = 10
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
      "-schwarzschild_Nmu" , Nmu, lflg , ierr) ;call CHKERR(ierr)

    associate( plex => solver%plex )

    call DMGetStratumIS(plex%ediff_dm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
    if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
    else
      allocate(dtau(plex%Nlay))
      allocate(Blev(plex%Nlay+1))
      allocate(Edn (plex%Nlay+1))
      allocate(Eup (plex%Nlay+1))

      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%horizface1_dm, plck_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geom_section, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)

      call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

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
        enddo

        do k=1,ke1
          call PetscSectionGetFieldOffset(plck_section, iface+k-1, i0, voff, ierr); call CHKERR(ierr)
          Blev(k) = xplck(i1+voff)
        enddo

        call schwarzschild(Nmu, dtau, xalbedo(i), Edn, Eup, Blev)

        do k = 0, ke1-2
          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          xediff(i1+voff) = Edn(i1+k)
          xediff(i2+voff) = Eup(i1+k)
        enddo
        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
        call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
        xediff(i1+voff) = Eup(i1+k)
        xediff(i2+voff) = Edn(i1+k)

        ! compute absorption as flux divergence
        do k = 1, ke1-1
          icell = cell_idx(k)
          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
          xabso(i1+voff) = Edn(k) - Edn(k+1) - Eup(k) + Eup(k+1)
          xabso(i1+voff) = xabso(i1+voff) / dz
        enddo
      enddo
      call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)

    endif ! TOA boundary ids

    !Schwarzschild solver returns fluxes as [W/m^2]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is not up to date
    solution%lchanged  = .True.

    end associate
  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine plexrt_twostream(solver, plex, kabs, ksca, g, albedo, sundir, solution, plck)
    class(t_plex_solver), intent(inout) :: solver
    type(t_plexgrid), intent(in) :: plex
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), intent(in), optional :: plck
    real(ireals), intent(in) :: sundir(:)
    type(t_state_container) :: solution

    real(ireals),allocatable :: vdtau(:), vw0(:), vg(:), Blev(:), Edir(:), Edn(:),Eup(:)

    type(tIS) :: boundary_ids
    integer(iintegers), pointer :: xitoa(:), cell_support(:)
    integer(iintegers), allocatable :: cell_idx(:)
    integer(iintegers) :: i, k, icell, iface, voff, ke1, geom_offset, idof
    real(ireals) :: dz, theta0, mu0
    real(ireals), pointer :: xksca(:), xkabs(:), xg(:), xalbedo(:), xplck(:)
    real(ireals), pointer :: xedir(:), xediff(:), xabso(:), xgeoms(:)
    type(tPetscSection) :: edir_section, ediff_section, abso_section, plck_section, geom_section

    real(ireals) :: dkabs, dksca, dg
    real(ireals) :: face_normal(3)
    logical :: lthermal, lsolar

    integer(mpiint) :: ierr

    if(solution%lsolar_rad) then
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad .and. norm2(sundir).le.zero) then
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
      call DMGetSection(plex%horizface1_dm, plck_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geom_section, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)

      call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

      if(lthermal) then
        if(.not.present(plck)) call CHKERR(1_mpiint, 'have to provide planck vec to compute thermal rad with twostr')
        allocate(Blev(plex%Nlay+1))
        call VecGetArrayReadF90(plck, xplck, ierr); call CHKERR(ierr)
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
        mu0 = cos(theta0)

        call get_consecutive_vertical_cell_idx(plex, icell, cell_idx)
        ke1 = size(cell_idx)+1
        do k=1,ke1-1
          icell = cell_idx(k)

          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          dkabs = xkabs(i1+icell)
          dksca = xksca(i1+icell)
          dg    = xg(i1+icell)

          call delta_scale( dkabs, dksca, dg, max_g=.65_ireals)

          vdtau(k) = (dkabs + dksca) * dz
          vw0(k)   = dksca / max(tiny(dksca), dkabs + dksca)
          vg(k)    = dg

          !call delta_scale_optprop(vdtau(k), vw0(k), vg(k), vg(k))
        enddo
        if(lthermal) then
          do k=1,ke1
            call PetscSectionGetFieldOffset(plck_section, iface+k-1, i0, voff, ierr); call CHKERR(ierr)
            Blev(k) = xplck(i1+voff)
          enddo
        endif

        if(solution%lsolar_rad) then
          call delta_eddington_twostream(vdtau,vw0,vg,&
            mu0,norm2(sundir)*mu0,xalbedo(i), &
            Edir, Edn, Eup )
        else
          call delta_eddington_twostream(vdtau,vw0,vg,&
            mu0,norm2(sundir)*mu0,xalbedo(i), &
            Edir, Edn, Eup, Blev)
        endif

        if(lsolar) then
          do k = 0, ke1-1
            call PetscSectionGetOffset(edir_section, iface+k, voff, ierr); call CHKERR(ierr)
            do idof = 1, solver%dirtop%dof
              xedir(voff+idof) = edir(i1+k)
            enddo
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

        ! compute absorption as flux divergence
        do k = 1, ke1-1
          icell = cell_idx(k)
          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
          if(lsolar) then
            xabso(i1+voff) = edn(k) - edn(k+1) - eup(k) + eup(k+1) + edir(k) - edir(k+1)
          else
            xabso(i1+voff) = edn(k) - edn(k+1) - eup(k) + eup(k+1)
          endif
          xabso(i1+voff) = xabso(i1+voff) / dz
        enddo
      enddo
      call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
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
      endif
    endif ! TOA boundary ids

    !Twostream solver returns fluxes as [W/m^2]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is not up to date
    solution%lchanged  = .False.
  end subroutine

  !> @brief wrapper for the disort solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine plexrt_disort(solver, plex, kabs, ksca, g, albedo, sundir, solution, plck)
    class(t_plex_solver), intent(inout) :: solver
    type(t_plexgrid), intent(in) :: plex
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), intent(in), optional :: plck
    real(ireals), intent(in) :: sundir(:)
    type(t_state_container) :: solution

    real,allocatable :: vdtau(:), vw0(:), vg(:), col_Bfrac(:) ! size nlay
    real,allocatable :: col_temper(:) ! size nlay+1
    real,allocatable :: RFLDIR(:), RFLDN(:), FLUP(:), DFDT(:), UAVG(:) ! size nlev
    real :: col_tskin, col_albedo, mu0

    type(tIS) :: boundary_ids
    integer(iintegers), pointer :: xitoa(:), cell_support(:)
    integer(iintegers), allocatable :: cell_idx(:)
    integer(iintegers) :: i, k, icell, iface, voff, ke1, geom_offset, idof
    real(ireals) :: dz, theta0
    real(ireals), pointer :: xksca(:), xkabs(:), xg(:), xalbedo(:), xplck(:)
    real(ireals), pointer :: xedir(:), xediff(:), xabso(:), xgeoms(:)
    type(tPetscSection) :: edir_section, ediff_section, abso_section, plck_section, geom_section

    real(ireals) :: dkabs, dksca, dg
    real(ireals) :: face_normal(3)
    logical :: lthermal, lsolar

    integer(mpiint) :: ierr

    integer(iintegers) :: nstreams
    logical :: ldelta_scale, lflg

    nstreams = 16
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
      "-disort_streams" , nstreams , lflg , ierr) ;call CHKERR(ierr)

    ldelta_scale = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
      "-disort_delta_scale" , ldelta_scale , lflg , ierr) ;call CHKERR(ierr)

    if(solution%lsolar_rad) then
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad .and. norm2(sundir).le.zero) then
      return
    endif

    lsolar = solution%lsolar_rad
    lthermal = .not.solution%lsolar_rad

    if(lthermal) call CHKERR(1_mpiint, 'currently cannot use the plexrt disort solver for thermal computations '// &
      'because we dont put the wavenumbers through')

    call DMGetStratumIS(plex%edir_dm, 'DomainBoundary', TOAFACE, boundary_ids, ierr); call CHKERR(ierr)
    if (boundary_ids.eq.PETSC_NULL_IS) then ! dont have TOA boundary faces
    else
      allocate(vdtau(plex%Nlay))
      allocate(vw0  (plex%Nlay))
      allocate(vg   (plex%Nlay))
      allocate(RFLDIR(plex%Nlay+1))
      allocate(RFLDN (plex%Nlay+1))
      allocate(FLUP  (plex%Nlay+1))
      allocate(DFDT  (plex%Nlay+1))
      allocate(UAVG  (plex%Nlay+1))

      allocate(col_temper(plex%Nlay+1), source=0.)
      allocate(col_Bfrac (plex%Nlay), source=1.)
      col_tskin = 0


      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%horizface1_dm, plck_section, ierr); call CHKERR(ierr)
      call DMGetSection(plex%geom_dm, geom_section, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, xgeoms, ierr); call CHKERR(ierr)

      call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%abso , xabso , ierr); call CHKERR(ierr)

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
        mu0 = real(cos(theta0))

        call get_consecutive_vertical_cell_idx(plex, icell, cell_idx)
        ke1 = size(cell_idx)+1
        do k=1,ke1-1
          icell = cell_idx(k)

          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          dkabs = xkabs(i1+icell)
          dksca = xksca(i1+icell)
          dg    = xg(i1+icell)

          if(ldelta_scale) call delta_scale( dkabs, dksca, dg, max_g=.65_ireals)

          vdtau(k) = real((dkabs + dksca) * dz)
          if(dkabs.le.tiny(dkabs)) then
            vw0 = 1
          else
            vw0(k) = real(dksca / (dkabs + dksca))
          endif
          vg(k)    = real(dg)
        enddo

        col_albedo = real(xalbedo(i))

        if(solution%lsolar_rad) then
            call default_flx_computation(&
              mu0, &
              real(norm2(sundir)), &
              col_albedo, &
              col_tskin, &
              .False., [0., 0.], col_Bfrac, &
              real(vdtau), &
              real(vw0),   &
              real(vg),    &
              col_temper, &
              RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
              int(nstreams), lverbose=.False.)
        endif

        if(lsolar) then
          do k = 0, ke1-1
            call PetscSectionGetOffset(edir_section, iface+k, voff, ierr); call CHKERR(ierr)
            do idof = 1, solver%dirtop%dof
              xedir(voff+idof) = RFLDIR(i1+k)
            enddo
          enddo
        endif

        do k = 0, ke1-2
          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          xediff(i1+voff) = RFLDN(i1+k)
          xediff(i1+voff+i1) = FLUP(i1+k)
        enddo
        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
        call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
        xediff(i1+voff) = FLUP(i1+k)
        xediff(i1+voff+i1) = RFLDN(i1+k)

        ! compute absorption as flux divergence
        do k = 1, ke1-1
          icell = cell_idx(k)
          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
          if(lsolar) then
            xabso(i1+voff) = RFLDN(k) - RFLDN(k+1) - FLUP(k) + FLUP(k+1) + RFLDIR(k) - RFLDIR(k+1)
          else
            xabso(i1+voff) = RFLDN(k) - RFLDN(k+1) - FLUP(k) + FLUP(k+1)
          endif
          xabso(i1+voff) = xabso(i1+voff) / dz
        enddo
      enddo
      call ISRestoreIndicesF90(boundary_ids, xitoa, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
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
      endif
    endif ! TOA boundary ids

    !Twostream solver returns fluxes as [W/m^2]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is not up to date
    solution%lchanged  = .False.
  end subroutine

end module
