module m_plex2rayli

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i8, zero, one, default_str_len

  use m_helper_functions, only: CHKERR, deg2rad, itoa

  use m_netcdfio, only: ncwrite

  use m_plex_grid, only: t_plexgrid, &
    TOAFACE

  use m_pprts_base, only : t_state_container
  use m_f2c_rayli, only: rfft_wedgeF90, rpt_img_wedgeF90

  implicit none

  logical, parameter :: ldebug=.False.

  contains

    subroutine dm3d_to_rayli_dmplex(dm, dmrayli)
      type(tDM), intent(in) :: dm
      type(tDM), intent(out) :: dmrayli

      integer(iintegers) :: dm_dim, Nedges, Nsidefaces, chartsize, pStart, pEnd
      integer(iintegers) :: cStart, cEnd, fStart, fEnd, eStart, eEnd, vStart, vEnd
      integer(iintegers) :: fxStart, fxEnd, exStart, exEnd, vxStart, vxEnd ! start and end indices of new dm
      integer(iintegers) :: i, k, iface, new_side_edge, kside
      integer(iintegers) :: faces(8), edges(3), verts(2)
      integer(iintegers), pointer :: faces_of_cell(:), edges_of_face(:), verts_of_edge(:)
      integer(iintegers), pointer :: transclosure(:)
      integer(iintegers), allocatable :: f2fx(:) ! offsets from original face into new faces dm
      integer(mpiint) :: comm, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)

      call DMGetDimension(dm, dm_dim, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(dm, 3_iintegers, cStart, cEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, 2_iintegers, fStart, fEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, 1_iintegers, eStart, eEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, 0_iintegers, vStart, vEnd, ierr); call CHKERR(ierr)
      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)

      allocate(f2fx(fStart:fEnd-1))
      k=fStart
      Nsidefaces = 0
      do iface = fStart, fEnd-1
          call DMPlexGetConeSize(dm, iface, Nedges, ierr); call CHKERR(ierr)
          f2fx(iface) = k
          if(ldebug) print *,'f2fx',iface,':',f2fx(iface)
          if(Nedges.eq.4) then
            Nsidefaces = Nsidefaces + 1
            k=k+2
          else
            k=k+1
          endif
      enddo

      fxStart = fStart
      fxEnd   = fEnd + Nsidefaces
      exStart = fxEnd
      exEnd   = eEnd + 2*Nsidefaces
      vxStart = exEnd
      vxEnd   = pEnd + 2*Nsidefaces

      chartsize = vxEnd ! have one face per side face more and one more edge
      if(ldebug) print *,'chartsize old:', pEnd, '=>', chartsize

      call DMPlexCreate(comm, dmrayli, ierr); call CHKERR(ierr)
      call DMSetDimension(dmrayli, dm_dim, ierr); call CHKERR(ierr)
      call DMPlexSetChart(dmrayli, i0, chartsize, ierr); call CHKERR(ierr)

      k = 0 ! k is running index for elements in DAG

      ! Every cell has 8 faces
      do i = cStart, cEnd-1
        call DMPlexSetConeSize(dmrayli, k, i8, ierr); call CHKERR(ierr)
        k = k+1
      enddo

      ! Every face has 3 edges
      do i = fxStart, fxEnd-1
        call DMPlexSetConeSize(dmrayli, k, i3, ierr); call CHKERR(ierr)
        k = k+1
      enddo

      ! Edges have 2 vertices
      do i = exStart, exEnd-1
        call DMPlexSetConeSize(dmrayli, k, i2, ierr); call CHKERR(ierr)
        k = k+1
      enddo

      call CHKERR(int((chartsize-(vEnd-vStart))-k, mpiint), 'This does not add up, we forgot something?')

      call DMSetUp(dmrayli, ierr); call CHKERR(ierr) ! Allocate space for cones

      do i = cStart, cEnd-1
        call DMPlexGetCone(dm, i, faces_of_cell, ierr); call CHKERR(ierr)
        faces(1) = f2fx(faces_of_cell(1))
        faces(2) = f2fx(faces_of_cell(2))

        faces(3) = f2fx(faces_of_cell(3))
        faces(4) = f2fx(faces_of_cell(3))+1

        faces(5) = f2fx(faces_of_cell(4))
        faces(6) = f2fx(faces_of_cell(4))+1

        faces(7) = f2fx(faces_of_cell(5))
        faces(8) = f2fx(faces_of_cell(5))+1
        call DMPlexRestoreCone(dm, i, faces_of_cell, ierr); call CHKERR(ierr)

        call DMPlexSetCone(dmrayli, i, faces, ierr); call CHKERR(ierr)
        if(ldebug) print *,'cell',i,'->',faces
      enddo

      k = fStart ! k is running index for faces in new dm
      kside = 0  ! is running index for sidefaces on original dm

      new_side_edge = 0
      do i = fStart, fEnd-1

        call DMPlexGetConeSize(dm, i, Nedges, ierr); call CHKERR(ierr)
        call DMPlexGetCone(dm, i, edges_of_face, ierr); call CHKERR(ierr)

        if(Nedges.eq.3) then
          edges = edges_of_face + Nsidefaces
          call DMPlexSetCone(dmrayli, k, edges, ierr); call CHKERR(ierr)
          if(ldebug) print *,'face',i,'=>', k, 'edges', edges

          k = k+1
        elseif(Nedges.eq.4) then
          new_side_edge = eEnd + kside ! we append cross face edges after the original-dm edges
          ! First subface
          edges = [edges_of_face(1), edges_of_face(4), new_side_edge] + Nsidefaces
          call DMPlexSetCone(dmrayli, k, edges, ierr); call CHKERR(ierr)
          if(ldebug) print *,'face',i,'=>', k, 'edges', edges

          k = k+1

          ! Second subface
          edges = [edges_of_face(2), edges_of_face(3), new_side_edge] + Nsidefaces
          call DMPlexSetCone(dmrayli, k, edges, ierr); call CHKERR(ierr)
          if(ldebug) print *,'face',i,'=>', k, 'edges', edges

          k = k+1
          kside = kside +1

          ! lets immediately setup the vertices of the new edge
          call DMPlexGetTransitiveClosure(dm, i, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
          verts = transclosure([11,17]) - vStart + vxStart
          call DMPlexSetCone(dmrayli, new_side_edge+Nsidefaces, verts, ierr); call CHKERR(ierr)
          if(ldebug) print *,'new edge', new_side_edge+Nsidefaces, 'verts', verts
          call DMPlexRestoreTransitiveClosure(dm, i, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)

        endif

        call DMPlexRestoreCone(dm, i, edges_of_face, ierr); call CHKERR(ierr)
      enddo

      k = exStart
      do i = eStart, eEnd-1
        call DMPlexGetCone(dm, i, verts_of_edge, ierr); call CHKERR(ierr)
        verts = verts_of_edge - vStart + vxStart
        call DMPlexRestoreCone(dm, i, verts_of_edge, ierr); call CHKERR(ierr)
        call DMPlexSetCone(dmrayli, k, verts, ierr); call CHKERR(ierr)
        if(ldebug) print *,'edge',k, 'verts', verts
        k = k+1
      enddo

      call DMPlexSymmetrize(dmrayli, ierr); call CHKERR(ierr)
      call DMPlexStratify(dmrayli, ierr); call CHKERR(ierr)
      if(ldebug) print *,'rayli dm connectivity done'

      call setup_coords()

      contains
        subroutine setup_coords()
          real(ireals), pointer :: coords_dm(:), coords_rayli(:)
          type(tVec)            :: vec_coord_dm, vec_coord_rayli
          type(tPetscSection)   :: coordSection_rayli
          integer(iintegers)    :: coordsize

          ! Create Coordinate stuff for rayli DM
          call DMGetCoordinateSection(dmrayli, coordSection_rayli, ierr); call CHKERR(ierr)
          call PetscSectionSetNumFields(coordSection_rayli, i1, ierr); call CHKERR(ierr) ! just 1 field for spherical coords
          call PetscSectionSetUp(coordSection_rayli, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldComponents(coordSection_rayli, i0, i3, ierr); call CHKERR(ierr)

          call PetscSectionSetChart(coordSection_rayli, vxStart, vxEnd, ierr);call CHKERR(ierr)
          do i = vxStart, vxEnd-1
            call PetscSectionSetDof(coordSection_rayli, i, i3, ierr); call CHKERR(ierr)
            call PetscSectionSetFieldDof(coordSection_rayli, i, i0, i3, ierr); call CHKERR(ierr)
          enddo
          call PetscSectionSetUp(coordSection_rayli, ierr); call CHKERR(ierr)
          call PetscSectionGetStorageSize(coordSection_rayli, coordSize, ierr); call CHKERR(ierr)

          ! Create New Vec to hold rayli coordinates
          call VecCreate(PETSC_COMM_SELF, vec_coord_rayli, ierr); call CHKERR(ierr)
          call VecSetSizes(vec_coord_rayli, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
          call VecSetBlockSize(vec_coord_rayli, i3, ierr);call CHKERR(ierr)
          call VecSetType(vec_coord_rayli, VECSTANDARD, ierr);call CHKERR(ierr)

          call PetscObjectSetName(vec_coord_rayli, "coordinates", ierr); call CHKERR(ierr)

          ! Fill rayli Coord Vec
          call VecGetArrayF90(vec_coord_rayli, coords_rayli, ierr); call CHKERR(ierr)

          call DMGetCoordinatesLocal(dm, vec_coord_dm, ierr); call CHKERR(ierr)
          call VecGetArrayReadF90(vec_coord_dm, coords_dm, ierr); call CHKERR(ierr)

          coords_rayli(:) = coords_dm(:)

          call VecRestoreArrayReadF90(vec_coord_dm, coords_dm, ierr); call CHKERR(ierr)
          call VecRestoreArrayF90(vec_coord_rayli, coords_rayli, ierr); call CHKERR(ierr)

          call DMSetCoordinatesLocal(dmrayli, vec_coord_rayli, ierr);call CHKERR(ierr)
          call VecDestroy(vec_coord_rayli, ierr); call CHKERR(ierr)
        end subroutine
    end subroutine

  !> @brief wrapper for the rayli MonteCarlo solver
  !> @details solve the radiative transfer equation with RayLi, currently only works for single task mpi runs
  subroutine rayli_wrapper(plex, kabs, ksca, g, albedo, sundir, solution, plck)
    use iso_c_binding
    type(t_plexgrid), intent(inout) :: plex
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), intent(in), optional :: plck
    real(ireals), intent(in) :: sundir(:)
    type(t_state_container), intent(inout) :: solution

    real(ireals), pointer :: xksca(:), xkabs(:), xg(:)

    logical :: lthermal, lsolar

    integer(mpiint) :: comm, myid, numnodes, ierr

    integer(iintegers) :: fStart, fEnd, cStart, cEnd, vStart, vEnd
    integer(iintegers) :: icell, iface, ivert, voff, idof
    integer(iintegers), pointer :: trans_closure(:), faces_of_cell(:)

    type(tVec) :: coordinates
    real(ireals), pointer :: coords(:)
    type(tPetscSection) :: coord_section

    integer(c_size_t), allocatable :: verts_of_face(:,:)
    integer(c_size_t), allocatable :: faces_of_wedges(:,:)
    real(c_double),    allocatable :: vert_coords(:,:)
    real(c_float),    allocatable :: rkabs(:), rksca(:), rg(:)
    real(c_float),    allocatable :: ralbedo_on_faces(:)
    real(c_float)                 :: rsundir(3)
    real(c_float),    allocatable :: flx_through_faces_edir(:)
    real(c_float),    allocatable :: flx_through_faces_ediff(:)
    real(c_float),    allocatable :: abso_in_cells(:)

    real(ireals) :: diffuse_point_origin(3)
    integer(c_size_t) :: Nphotons, Nwedges, Nfaces, Nverts
    real(ireals) :: opt_photons

    integer(c_size_t) :: outer_id
    logical :: lcyclic, lflg
    integer(c_int) :: icyclic

    opt_photons = 100000
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-rayli_photons", opt_photons, lflg,ierr) ; call CHKERR(ierr)
    Nphotons = int(opt_photons, c_size_t)

    call PetscObjectGetComm(plex%dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(numnodes.gt.1) call CHKERR(numnodes, "Rayli currently only in single rank computations available")

    lsolar = solution%lsolar_rad
    lthermal = .not.solution%lsolar_rad

    if(present(plck).or.lthermal) &
      call CHKERR(1_mpiint, "You provided planck stuff, I guess you want to use thermal radiation computations."// &
                            "However, Rayli currently only supports solar radiation.")

    diffuse_point_origin = 0; idof=3
    call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-rayli_diff_flx_origin',diffuse_point_origin, idof, lflg, ierr); call CHKERR(ierr)
    if(lflg) call CHKERR(int(idof-i3, mpiint), 'must provide exactly 3 values for rayli_diff_flx_origin')

    lcyclic=.False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-rayli_cylic_bc', lcyclic, lflg, ierr); call CHKERR(ierr)
    if(lcyclic) then
      icyclic = 1
    else
      icyclic = 0
    endif

    if(solution%lsolar_rad) then
      call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    endif
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad .and. norm2(sundir).le.zero) then
      return
    endif

    call DMPlexGetDepthStratum(plex%dm, i3, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetDepthStratum(plex%dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
    call DMPlexGetDepthStratum(plex%dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
    Nwedges = cEnd - cStart
    Nfaces = fEnd - fStart
    Nverts = vEnd - vStart
    !print *,'Rayli Nwedges', Nwedges, 'Nfaces', Nfaces, 'Nverts', Nverts

    outer_id = -1

    allocate(verts_of_face(4, Nfaces), &
             faces_of_wedges(5, Nwedges), &
             vert_coords(3, Nverts), &
             rkabs(Nwedges), &
             rksca(Nwedges), &
             rg   (Nwedges), &
             ralbedo_on_faces(Nfaces), &
             flx_through_faces_edir(Nfaces), &
             flx_through_faces_ediff(2*Nfaces), &
             abso_in_cells(Nwedges) )

    do icell = cStart, cEnd-1
      call DMPlexGetCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      faces_of_wedges(:,icell-cStart+1) = faces_of_cell - fStart
      call DMPlexRestoreCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    do iface = fStart, fEnd-1
      call DMPlexGetTransitiveClosure(plex%dm, iface, PETSC_TRUE, trans_closure, ierr); call CHKERR(ierr)
      select case (size(trans_closure))
      case(14) ! 3 edges
        verts_of_face(1:3, iface-fStart+1) = trans_closure(9:size(trans_closure):2) - vStart
        verts_of_face(4  , iface-fStart+1) = outer_id
      case(18) ! 4 edges
        verts_of_face(:, iface-fStart+1) = trans_closure(11:size(trans_closure):2) - vStart
      end select
      call DMPlexRestoreTransitiveClosure(plex%dm, iface, PETSC_TRUE, trans_closure, ierr); call CHKERR(ierr)
    enddo

    call DMGetCoordinateSection(plex%dm, coord_section, ierr); call CHKERR(ierr)
    call DMGetCoordinatesLocal(plex%dm, coordinates, ierr); call CHKERR(ierr)
    call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

    do ivert = vStart, vEnd-1
      call PetscSectionGetOffset(coord_section, ivert, voff, ierr); call CHKERR(ierr)
      vert_coords(:, ivert-vStart+1) = coords(voff+1:voff+3)
    enddo
    call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    rkabs(:) = real(xkabs, kind(rkabs))
    rksca(:) = real(xksca, kind(rksca))
    rg   (:) = real(xg,    kind(rg   ))
    !call delta_scale(rkabs, rksca, rg, max_g=0._c_double)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call fill_albedo(plex, albedo, ralbedo_on_faces)

    rsundir = real(-sundir, kind(rsundir))

    call take_snap(comm, &
      Nphotons, Nwedges, Nfaces, Nverts, &
      verts_of_face, faces_of_wedges, vert_coords, &
      rkabs, rksca, rg, &
      ralbedo_on_faces, &
      rsundir, &
      solution, ierr)
    if(ierr.eq.1) return

    ierr = rfft_wedgeF90(Nphotons, Nwedges, Nfaces, Nverts, icyclic, &
      verts_of_face, faces_of_wedges, vert_coords, &
      rkabs, rksca, rg, &
      ralbedo_on_faces, rsundir, real(diffuse_point_origin, c_float), &
      flx_through_faces_edir, flx_through_faces_ediff, abso_in_cells ); call CHKERR(ierr)

    call rayli_get_result(plex, &
          flx_through_faces_edir, &
          flx_through_faces_ediff, &
          abso_in_cells, &
          solution)

    !Twostream solver returns fluxes as [W/m^2]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is up to date, otherwise abso would be overwritten
    solution%lchanged  = .False.

    end subroutine

    subroutine rayli_get_result(plex, &
        flx_through_faces_edir, &
        flx_through_faces_ediff, &
        abso_in_cells, &
        solution)
      type(t_plexgrid), intent(inout) :: plex
      real(c_float), intent(in) :: flx_through_faces_edir(:)
      real(c_float), intent(in) :: flx_through_faces_ediff(:)
      real(c_float), intent(in) :: abso_in_cells(:)
      type(t_state_container), intent(inout) :: solution

      type(tIS) :: toa_ids
      integer(iintegers) :: i, k, ke1, Ncol, ridx, numDof, geom_offset
      integer(iintegers) :: icell, iface, idof, voff
      integer(iintegers) :: ofStart, ofEnd, cStart, cEnd
      integer(iintegers), pointer :: xtoa_faces(:)
      real(ireals), pointer :: xedir(:), xediff(:), xabso(:), geoms(:)
      real(ireals) :: area
      type(tPetscSection) :: edir_section, ediff_section, abso_section, geomSection
      integer(mpiint) :: ierr

      call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      ! --------------------------- Direct Radiation ------------------
      call DMPlexGetDepthStratum(plex%dm, i2, ofStart, ofEnd, ierr); call CHKERR(ierr) ! vertices
      if(solution%lsolar_rad) then
        call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
        call VecGetArrayF90(solution%edir, xedir, ierr); call CHKERR(ierr)
        do iface = ofStart, ofEnd-1
          call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
          area = geoms(i1+geom_offset)

          call PetscSectionGetOffset(edir_section, iface, voff, ierr); call CHKERR(ierr)
          call PetscSectionGetDof(edir_section, iface, numDof, ierr); call CHKERR(ierr)
          do idof = 1, numDof
            xedir(voff+idof) = real(abs(flx_through_faces_edir(iface-ofStart+1)), ireals)! / area
          enddo
        enddo
        call VecRestoreArrayF90(solution%edir, xedir, ierr); call CHKERR(ierr)
      endif

      ! --------------------------- Diffuse Radiation -----------------
      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)

      call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', &
        TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

      ke1 = plex%Nlay+1

      call ISGetIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
      ridx=1
      do i = 1, size(xtoa_faces)
        iface = xtoa_faces(i)
        do k = 0, ke1-2
          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
          call PetscSectionGetFieldOffset(geomSection, iface+k, i2, geom_offset, ierr); call CHKERR(ierr)
          area = geoms(i1+geom_offset)

          xediff(i1+voff) = real(abs( flx_through_faces_ediff(ridx  ) ), ireals)! / area
          xediff(i2+voff) = real(abs( flx_through_faces_ediff(ridx+1) ), ireals)! / area

          ridx = ridx+2
        enddo

        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
        call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(geomSection, iface+ke1-1, i2, geom_offset, ierr); call CHKERR(ierr)
        area = geoms(i1+geom_offset)

        xediff(i2+voff) = real(abs( flx_through_faces_ediff(ridx  ) ), ireals)! / area
        xediff(i1+voff) = real(abs( flx_through_faces_ediff(ridx+1) ), ireals)! / area
        ridx = ridx+2
      enddo

      call ISRestoreIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)

      ! --------------------------- Absorption ------------------------
      call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
      call VecGetArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      do icell = cStart, cEnd-1
        call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr) ! cell volume
        call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
        xabso(i1+voff) = real(abso_in_cells(i1+icell-cStart), ireals) / geoms(i1+geom_offset)
      enddo
      call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    end subroutine

    subroutine fill_albedo(plex, albedo, ralbedo_on_faces)
      type(t_plexgrid), intent(inout) :: plex
      type(tVec), intent(in) :: albedo
      real(c_float), intent(out) :: ralbedo_on_faces(:)

      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xtoa_faces(:)
      real(ireals), pointer :: xalbedo(:)
      integer(iintegers) :: srfc_face, i, fStart, fEnd
      integer(mpiint) :: ierr

      call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      call ISGetIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(plex%dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces

      call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)

      ralbedo_on_faces(:) = -one
      do i=1,size(xtoa_faces)
        srfc_face = xtoa_faces(i) + plex%Nlay
        ralbedo_on_faces(srfc_face-fStart+1) = real(xalbedo(i), kind(ralbedo_on_faces))
      enddo

      call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
    end subroutine

    subroutine take_snap(comm, &
        Nphotons, Nwedges, Nfaces, Nverts, &
        verts_of_face, faces_of_wedges, vert_coords, &
        rkabs, rksca, rg, &
        ralbedo_on_faces, &
        rsundir, &
        solution, ierr)
      integer(mpiint),   intent(in) :: comm
      integer(c_size_t), intent(in) :: Nphotons, Nwedges, Nfaces, Nverts
      integer(c_size_t), intent(in) :: verts_of_face(:,:)
      integer(c_size_t), intent(in) :: faces_of_wedges(:,:)
      real(c_double),    intent(in) :: vert_coords(:,:)
      real(c_float),     intent(in) :: rkabs(:), rksca(:), rg(:)
      real(c_float),     intent(in) :: ralbedo_on_faces(:)
      real(c_float),     intent(in) :: rsundir(3)
      type(t_state_container), intent(in) :: solution
      integer(mpiint), intent(out) :: ierr

      character(len=default_str_len) :: snap_path, groups(2)
      logical :: lflg
      integer(c_size_t) :: Nx=400, Ny=300
      real(c_float), allocatable :: img(:,:)
      real(c_float) :: cam_loc(3), cam_viewing_dir(3), cam_up_vec(3)
      real(c_float) :: fov_width, fov_height
      real(ireals), dimension(3) :: visit_focus, visit_view_normal, visit_view_up
      real(ireals) :: visit_view_angle, visit_image_zoom, visit_parallel_scale
      integer(iintegers) :: narg

      integer(mpiint) :: myid

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-rayli_snapshot", &
        snap_path, lflg,ierr) ; call CHKERR(ierr)
      if(lflg) then
        if(len_trim(snap_path).eq.0) snap_path = 'rayli_snaphots.nc'
        if(myid.eq.0) print *,'Capturing scene to file: '//trim(snap_path), len_trim(snap_path)
        Nx = 400
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-rayli_snap_Nx', narg, lflg, ierr); call CHKERR(ierr)
        if(lflg) Nx = narg

        Ny = 300
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-rayli_snap_Ny', narg, lflg, ierr); call CHKERR(ierr)
        if(lflg) Ny = narg

        allocate(img(Nx, Ny))

        visit_view_angle = 30._ireals
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_view_angle', visit_view_angle, lflg, ierr); call CHKERR(ierr)

        visit_image_zoom = 1._ireals
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_image_zoom', visit_image_zoom, lflg, ierr); call CHKERR(ierr)

        visit_parallel_scale = 1e5_ireals
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_parallel_scale', visit_parallel_scale, lflg, ierr); call CHKERR(ierr)

        visit_focus = 0._ireals
        narg=3
        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_focus', visit_focus, narg, lflg, ierr); call CHKERR(ierr)
        if(lflg) &
          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')

        visit_view_normal = -rsundir
        narg=3
        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_view_normal', visit_view_normal, narg, lflg, ierr); call CHKERR(ierr)
        if(lflg) &
          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')

        visit_view_up = [0,0,1]
        narg=3
        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
          '-visit_view_up', visit_view_up, narg, lflg, ierr); call CHKERR(ierr)
        if(lflg) &
          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')

        cam_viewing_dir = -real(visit_view_normal, kind(cam_viewing_dir))
        cam_up_vec = real(visit_view_up, kind(cam_up_vec))
        cam_loc = real(visit_focus + visit_view_normal * visit_parallel_scale &
          / tan(deg2rad(visit_view_angle)/2), kind(cam_loc))

        fov_width = 2 * real(tan(deg2rad(visit_view_angle)/2) / visit_image_zoom, kind(fov_width))
        fov_height = fov_width * real(Ny, c_float) / real(Nx, c_float)

        ierr = rpt_img_wedgeF90( Nx, Ny, &
          Nphotons, Nwedges, Nfaces, Nverts, &
          verts_of_face, faces_of_wedges, vert_coords, &
          rkabs, rksca, rg, &
          ralbedo_on_faces, &
          rsundir, & ! DEBUG note the kabs/ksca/
          cam_loc, cam_viewing_dir, cam_up_vec, &
          fov_width, fov_height, &
          img); call CHKERR(ierr)

        groups(1) = trim(snap_path)
        groups(2) = "rpt_img_"//itoa(solution%uid)
        call ncwrite(groups, img, ierr); call CHKERR(ierr)
        ierr = 1
        return
      endif
    end subroutine
end module
