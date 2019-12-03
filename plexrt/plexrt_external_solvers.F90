module m_plexrt_external_solvers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, zero
  use m_helper_functions, only: CHKERR, CHKWARN, &
    angle_between_two_vec, delta_scale, &
    cstr, itoa, ftoa, approx

  use m_pprts_base, only : t_state_container
  use m_plex_rt_base, only: t_plex_solver

  use m_plex_grid, only: t_plexgrid, &
    get_inward_face_normal, &
    get_consecutive_vertical_cell_idx, &
    TOAFACE, INNERSIDEFACE, dmplex_set_new_section

  use m_schwarzschild, only: schwarzschild, B_eff
  use m_twostream, only: delta_eddington_twostream
  use m_tenstr_disort, only: default_flx_computation

  use m_plexrt_nca, only: plexrt_nca_init, plexrt_nca
  use m_icon_plex_utils, only: dump_ownership, dmplex_2d_to_3d

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
    if(.not.allocated(solver%plck)) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

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
    solution%lchanged  = .False.

    end associate
  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine plexrt_twostream(solver, plex, lthermal, lsolar, kabs, ksca, g, albedo, sundir, solution, plck)
    class(t_plex_solver), intent(inout) :: solver
    type(t_plexgrid), intent(in) :: plex
    logical, intent(in) :: lthermal, lsolar
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), allocatable, intent(in), optional :: plck
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

    integer(mpiint) :: ierr

    if(lthermal) then
      if(.not.present(plck)) call CHKERR(1_mpiint, 'planck vec has to be present if thermal computations are wanted')
      if(.not.allocated(plck)) call CHKERR(1_mpiint, 'planck vec is present but not allocated')
    endif

    if(allocated(solution%edir)) then
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)


    if(.not.lthermal .and. (lsolar .and. norm2(sundir).le.zero)) then
      return
    endif


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
        allocate(Blev(plex%Nlay+1))
        call VecGetArrayReadF90(plck, xplck, ierr); call CHKERR(ierr)
      endif

      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
      if(lsolar) then
        call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
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
            call PetscSectionGetOffset(plck_section, iface+k-1, voff, ierr); call CHKERR(ierr)
            Blev(k) = xplck(i1+voff)
          enddo
        endif

        if(lthermal) then
          call delta_eddington_twostream(vdtau,vw0,vg,&
            mu0,norm2(sundir)*mu0,xalbedo(i), &
            Edir, Edn, Eup, Blev)
        else
          call delta_eddington_twostream(vdtau,vw0,vg,&
            mu0,norm2(sundir)*mu0,xalbedo(i), &
            Edir, Edn, Eup )
        endif

        if(lsolar) then
          do k = 0, ke1-1
            call PetscSectionGetOffset(edir_section, iface+k, voff, ierr); call CHKERR(ierr)
            do idof = 1, solver%dirtop%dof
              xedir(voff+idof) = Edir(i1+k)
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


        if(lsolar .and. (.not.approx( (Edir(ke1)+Edn(ke1)) * xalbedo(i), Eup(ke1), sqrt(epsilon(Eup)) ))) &
          call CHKERR(1_mpiint, 'Reflected Radiation at the surface does not match the given '//new_line('')// &
          'albedo: '//ftoa(xalbedo(i))//new_line('')// &
          'Edir: '//ftoa(Edir(ke1))//new_line('')// &
          'Edn: '//ftoa(Edn(ke1))//new_line('')// &
          'Eup: '//ftoa(Eup(ke1))//new_line('')// &
          'should be: '//ftoa((Edir(ke1)+Edn(ke1))*xalbedo(i))//new_line('') )

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
    ! and mark solution that it is up to date
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

  !> @brief wrapper to apply NCA on (preferably 1D) solutions
  subroutine plexrt_NCA_wrapper(solver, solution, ierr)
    class(t_plex_solver), intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution
    integer(mpiint) :: comm, myid, ierr


    type(tPetscSection) :: geom_section, ediff_section, plck_section, abso_section, nca_section

    integer(iintegers) :: iface

    real(ireals), pointer :: xgeoms(:), xkabs(:), xplck(:), xediff(:), xabso(:), xnca(:)

    type(tVec) :: lnca ! has 5 dof on cells for ( Edn_top, Eup_top, Edn_bot, Eup_bot, kabs )

    if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling NCA for solar calculation!')
    if(.not.solution%lWm2_diff) call CHKERR(1_mpiint, 'Tried to compute NCA on solution which is not in [W per m2]')
    ierr = 0

    call PetscObjectGetComm(solver%plex%abso_dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call plexrt_nca_init(comm)

    if(.not.allocated(solver%plex%nca_dm)) call generate_NCA_dm(solver%plex)
    call DMGetLocalVector(solver%plex%nca_dm, lnca, ierr); call CHKERR(ierr)

    call DMGetSection(solver%plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
    call DMGetSection(solver%plex%horizface1_dm, plck_section, ierr); call CHKERR(ierr)
    call DMGetSection(solver%plex%geom_dm, geom_section, ierr); call CHKERR(ierr)
    call DMGetSection(solver%plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
    call DMGetSection(solver%plex%nca_dm, nca_section, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(solver%plex%geomVec, xgeoms, ierr); call CHKERR(ierr)
    call VecGetArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)

    call fill_nca_vec(solver%plex, lnca)

    call VecGetArrayReadF90(lnca, xnca, ierr); call CHKERR(ierr)
    call compute_nca()
    call VecRestoreArrayReadF90(lnca, xnca, ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(solver%kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%plck, xplck, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%plex%geomVec, xgeoms, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)

    contains
      subroutine compute_nca()
        ! NCA Variables
        real(ireals) :: dx(3)                  ! edge lengths of triangle: dx1, dx2, dx3
        real(ireals) :: areas(5)               ! area of side faces and top and bot faces
        real(ireals) :: dz, volume             ! height of grid box and volume
        real(ireals) :: base_info(7)           ! current voxel ( kabs, kabs_top, Ldn_top, Btop, kabs_bot, Lup_bot, Bbot )
        real(ireals) :: side_info(15)          ! side voxels dim(3 * 5) ( kabs, Edn_top, Eup_top, Edn_bot, Eup_bot )
        real(ireals) :: nca_info(5)
        real(ireals) :: hr ! new 3D heating rate in the voxel

        integer(iintegers) :: geom_offset
        integer(iintegers) :: icell, iface, iside, iedge, isup, nca_iface, neigh_cell, cStart, cEnd
        integer(iintegers) :: top_face, bot_face, voff
        integer(iintegers), pointer :: faces_of_cell(:), nca_faces_of_cell(:), cell_support(:)

        integer(iintegers), target :: points(2)
        integer(iintegers), pointer :: ppoints(:), coveredPoints(:)

        call DMPlexGetHeightStratum(solver%plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
        do icell = cStart, cEnd-1
          call PetscSectionGetFieldOffset(geom_section, icell, i2, geom_offset, ierr); call CHKERR(ierr)
          volume = xgeoms(i1+geom_offset)
          call PetscSectionGetFieldOffset(geom_section, icell, i3, geom_offset, ierr); call CHKERR(ierr)
          dz = xgeoms(i1+geom_offset)

          call DMPlexGetCone(solver%plex%abso_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
          call DMPlexGetCone(solver%plex%nca_dm, icell, nca_faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
          top_face = faces_of_cell(1)
          bot_face = faces_of_cell(2)

          call PetscSectionGetFieldOffset(plck_section, top_face, i0, voff, ierr); call CHKERR(ierr)
          base_info(4) = xplck(i1+voff)
          call PetscSectionGetFieldOffset(plck_section, bot_face, i0, voff, ierr); call CHKERR(ierr)
          base_info(7) = xplck(i1+voff)

          call PetscSectionGetFieldOffset(plck_section, bot_face, i0, voff, ierr); call CHKERR(ierr)

          do iside = 1, size(faces_of_cell)
            call PetscSectionGetFieldOffset(geom_section, faces_of_cell(iside), i2, geom_offset, ierr); call CHKERR(ierr)
            areas(iside) = xgeoms(i1+geom_offset)
          enddo

          ! Find side_info
          do iside=1,3
            iface = faces_of_cell(2+iside)
            nca_iface = nca_faces_of_cell(2+iside)

            ! Determine edge length of edge between base face and upper face triangle
            ppoints => points
            points = [top_face, iface]
            call DMPlexGetMeet(solver%plex%abso_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
            iedge = coveredPoints(1)
            call DMPlexRestoreMeet(solver%plex%abso_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
            call PetscSectionGetOffset(geom_section, iedge, geom_offset, ierr); call CHKERR(ierr)
            dx(iside) = xgeoms(i1+geom_offset)

            ! Determine edge length of edge between base face and lower face triangle
            ppoints => points
            points = [bot_face, iface]
            call DMPlexGetMeet(solver%plex%abso_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
            iedge = coveredPoints(1)
            call DMPlexRestoreMeet(solver%plex%abso_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
            call PetscSectionGetOffset(geom_section, iedge, geom_offset, ierr); call CHKERR(ierr)
            dx(iside) = ( dx(iside) + xgeoms(i1+geom_offset) ) / 2

            call get_neigh_face_info(icell, nca_iface, side_info( (iside-1)*5 + i1 : iside*5))
          enddo

          ! Find base_info
          call get_neigh_voxel_info(icell, nca_info)
          base_info(1) = nca_info(i1)
          base_info(3) = nca_info(i2)
          base_info(6) = nca_info(i5)

          neigh_cell = icell ! if we dont have a neighbor cell above(TOA), use this cell`s properties
          call DMPlexGetSupport(solver%plex%nca_dm, nca_faces_of_cell(1), cell_support, ierr); call CHKERR(ierr)
          do isup=1,size(cell_support)
            if(cell_support(isup).ne.icell) neigh_cell = cell_support(isup)
          enddo
          base_info(2) = nca_info(i1)

          neigh_cell = icell ! if we dont have a neighbor cell below(surface), use this cell`s properties
          call DMPlexGetSupport(solver%plex%nca_dm, nca_faces_of_cell(2), cell_support, ierr); call CHKERR(ierr)
          do isup=1,size(cell_support)
            if(cell_support(isup).ne.icell) neigh_cell = cell_support(isup)
          enddo
          base_info(5) = nca_info(i1)

          call DMPlexRestoreCone(solver%plex%abso_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
          call DMPlexRestoreCone(solver%plex%nca_dm, icell, nca_faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

          call plexrt_nca (dx(1), dx(2), dx(3), dz, &
            areas(1), areas(2), areas(3), areas(4), areas(5), volume, &
            base_info, side_info, hr)

          call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
          !print *, cstr(itoa(myid)//' icell '//itoa(icell)//' before: '//ftoa(xabso(i1+voff))//' after '//ftoa(hr), 'red')
          xabso(i1+voff) = hr
        enddo
      end subroutine

      subroutine get_neigh_face_info(base_cell, iface, info)
        integer(iintegers), intent(in) :: base_cell, iface
        real(ireals), intent(out) :: info(5) ! ( kabs, Edn_top, Eup_top, Edn_bot, Eup_bot )

        integer(iintegers), pointer :: cell_support(:)
        integer(iintegers) :: isup, neigh_cell, ibnd_type, owner, nca_offset

        info = -9
        neigh_cell = base_cell ! if we dont have a neighbor, use this cell`s properties
        call DMPlexGetSupport(solver%plex%nca_dm, iface, cell_support, ierr); call CHKERR(ierr)
        do isup=1,size(cell_support)
          if(cell_support(isup).ne.base_cell) neigh_cell = cell_support(isup)
        enddo
        call DMPlexRestoreSupport(solver%plex%nca_dm, iface, cell_support, ierr); call CHKERR(ierr)
        if(neigh_cell.eq.base_cell) then ! if we dont have a neighboring cell, check if we are at a boundary
          call DMLabelGetValue(solver%plex%domainboundarylabel, iface, ibnd_type, ierr); call CHKERR(ierr)
          if(ibnd_type.eq.INNERSIDEFACE) then
            call DMLabelGetValue(solver%plex%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
            if(myid.eq.owner) then
              call PetscSectionGetFieldOffset(nca_section, iface, i1, nca_offset, ierr); call CHKERR(ierr)
            else
              call PetscSectionGetFieldOffset(nca_section, iface, i0, nca_offset, ierr); call CHKERR(ierr)
            endif
            info = xnca(nca_offset+i1:nca_offset+i5)
            !print *,myid, 'Found innsersideface boundary, retrieving data from nca_vec', iface, info
            return
          endif
        endif
        call get_neigh_voxel_info(neigh_cell, info)
        if(any(info.lt.0)) call CHKERR(1_mpiint, 'found negative val in neigh_info')
      end subroutine

      subroutine get_neigh_voxel_info(icell, info)
        integer(iintegers), intent(in) :: icell
        real(ireals), intent(out) :: info(5) ! ( kabs, Edn_top, Eup_top, Edn_bot, Eup_bot )
        integer(iintegers),pointer :: faces_of_cell(:)
        integer(iintegers) :: top_face, bot_face, voff

        info(i1) = xkabs(i1+icell)

        call DMPlexGetCone(solver%plex%nca_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
        top_face = faces_of_cell(1)
        bot_face = faces_of_cell(2)
        call DMPlexRestoreCone(solver%plex%nca_dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

        call PetscSectionGetFieldOffset(ediff_section, top_face, i0, voff, ierr); call CHKERR(ierr)
        info(i2) = xediff(i1+voff) ! Edn_top
        info(i3) = xediff(i2+voff) ! Eup_top

        call PetscSectionGetFieldOffset(ediff_section, bot_face, i0, voff, ierr); call CHKERR(ierr)
        if(solver%plex%zindex(bot_face).eq.solver%plex%Nlay+1) then ! surface flux is reversed edn/eup idx
          info(i4) = xediff(i2+voff) ! Edn_bot
          info(i5) = xediff(i1+voff) ! Eup_bot
        else
          info(i4) = xediff(i1+voff) ! Edn_bot
          info(i5) = xediff(i2+voff) ! Eup_bot
        endif
      end subroutine

      subroutine fill_nca_vec(plex, lnca)
        class(t_plexgrid), intent(in) :: plex
        type(tVec), intent(inout) :: lnca ! has 5 dof on cells for ( kabs, Edn_top, Eup_top, Edn_bot, Eup_bot )
        ! Set the side information into the nca_dm global vec and copy that stuff over to neighboring processes
        type(tVec) :: gnca ! global nca vec
        type(tIS) :: boundary_ids
        integer(iintegers), pointer :: xbndry_iface(:), cell_support(:)
        integer(iintegers) :: i, icell
        integer(iintegers) :: owner, nca_offset

        call DMGetLocalVector(solver%plex%nca_dm, lnca, ierr); call CHKERR(ierr)
        call VecSet(lnca, zero, ierr); call CHKERR(ierr)
        call VecGetArrayF90(lnca, xnca, ierr); call CHKERR(ierr)

        call DMGetStratumIS(solver%plex%nca_dm, 'DomainBoundary', INNERSIDEFACE, boundary_ids, ierr); call CHKERR(ierr)
        if (boundary_ids.ne.PETSC_NULL_IS) then
          call ISGetIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)
          do i = 1, size(xbndry_iface)
            iface = xbndry_iface(i)

            call DMLabelGetValue(plex%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
            if(myid.eq.owner) then
              call PetscSectionGetFieldOffset(nca_section, iface, i0, nca_offset, ierr); call CHKERR(ierr)
            else
              call PetscSectionGetFieldOffset(nca_section, iface, i1, nca_offset, ierr); call CHKERR(ierr)
            endif
            !print *, myid, 'boundary face', iface, owner,'nca_offset', nca_offset, 'shape', shape(xnca)

            call DMPlexGetSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr)
            icell = cell_support(1)
            call DMPlexRestoreSupport(plex%ediff_dm, iface, cell_support, ierr); call CHKERR(ierr)

            call get_neigh_voxel_info(icell, xnca(nca_offset+i1:nca_offset+i5))
          enddo
          call ISRestoreIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)
        endif

      call VecRestoreArrayF90(lnca, xnca, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(solver%plex%nca_dm, gnca, ierr); call CHKERR(ierr)
      call VecSet(gnca, zero, ierr); call CHKERR(ierr)

      call DMLocalToGlobalBegin(solver%plex%nca_dm, lnca, ADD_VALUES, gnca, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (solver%plex%nca_dm, lnca, ADD_VALUES, gnca, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(solver%plex%nca_dm, gnca, INSERT_VALUES, lnca, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (solver%plex%nca_dm, gnca, INSERT_VALUES, lnca, ierr); call CHKERR(ierr)
      call DMRestoreGlobalVector(solver%plex%nca_dm, gnca, ierr); call CHKERR(ierr)
      end subroutine

      subroutine generate_NCA_dm(plex)
        class(t_plexgrid), intent(inout) :: plex
        integer(mpiint) :: ierr
        !type(tDM) :: dm2d, dm2d_redist
        !integer(iintegers), allocatable :: zindex(:)

        if(.not.allocated(plex%dm)) call CHKERR(1_mpiint, 'parent dm not allocated')
        if(plex%dm.eq.PETSC_NULL_DM) call CHKERR(1_mpiint, 'parent dm is null_dm')
        if(allocated(plex%nca_dm)) call CHKERR(1_mpiint, 'nca_dm already allocated')

        allocate(plex%nca_dm)
        ! TODO: dmplex partitioner redistributes the underlying mesh.
        !       Therefore, the cell ownerships may change and consequently the meshes dont have the same cell ids
        !       This is not clear how to manage this. i.e. at the moment NCA uses something like a hack...
        !       until I come up with a good strategy to do so...
        !       Instead of using a proper overlap dmplex, we use a face section with 2 fields on the inner domain boundary
        !       and copy the needed values through that, i.e. with a LocalToLocal update
        !       Field 0 is for values from owner -> non-owner and Field 1 vice versa
        call DMClone(plex%dm, plex%nca_dm, ierr); call CHKERR(ierr)

        call dmplex_set_new_section(plex%nca_dm, 'NCA Section', i2, &
          [i0,i0], [i5,i5], [i0,i0], [i0,i0])
      end subroutine
  end subroutine
end module
