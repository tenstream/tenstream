module m_icon_plex_utils

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, i0, i1, i2, i3, i4, i5, zero, one

  use m_helper_functions, only: chkerr, itoa, norm

  use m_plex_grid, only: print_dmplex, create_plex_section

  implicit none

  private
  public :: dmplex_2D_to_3D, create_2d_fish_plex, dump_ownership

  logical, parameter :: ldebug=.True.

  contains

    subroutine dmplex_2D_to_3D(dm2d, hhl, dm3d)
      type(tDM), intent(in) :: dm2d
      real(ireals), intent(in) :: hhl(:) ! height levels of interfaces, those will be added to base height of 2D elements
      type(tDM), intent(out) :: dm3d

      integer(iintegers) :: p2dStart, p2dEnd
      integer(iintegers) :: f2dStart, f2dEnd
      integer(iintegers) :: e2dStart, e2dEnd
      integer(iintegers) :: v2dStart, v2dEnd
      integer(iintegers) :: ke, ke1, chartsize
      integer(iintegers) :: Nfaces2d, Nedges2d, Nverts2d
      integer(iintegers) :: Ncells, Nfaces, Nedges, Nverts
      integer(mpiint) :: comm, ierr
      logical, parameter :: ldebug=.False.
      type(tDMLabel) :: boundarylabel

      ke1 = size(hhl)
      ke = ke1-1

      call PetscObjectGetComm(dm2d, comm, ierr); call CHKERR(ierr)

      call print_dmplex(comm, dm2d)

      call DMPlexGetChart(dm2d, p2dStart, p2dEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm2d, i0, f2dStart, f2dEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetHeightStratum(dm2d, i1, e2dStart, e2dEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetHeightStratum(dm2d, i2, v2dStart, v2dEnd, ierr); call CHKERR(ierr) ! vertices

      Nfaces2d = f2dEnd - f2dStart
      Nedges2d = e2dEnd - e2dStart
      Nverts2d = v2dEnd - v2dStart
      print *,'Nfaces2d', Nfaces2d
      print *,'Nedges2d', Nedges2d
      print *,'Nverts2d', Nverts2d

      Ncells = Nfaces2d * ke
      Nfaces = Nfaces2d * ke1 + Nedges2d * ke
      Nedges = Nedges2d * ke1 + Nverts2d * ke
      Nverts = Nverts2d * ke1

      chartsize = Ncells + Nfaces + Nedges + Nverts

      print *,'Ncells', Ncells
      print *,'Nfaces', Nfaces
      print *,'Nedges', Nedges
      print *,'Nverts', Nverts
      print *,'Chartsize', chartsize

      call DMPlexCreate(comm, dm3d, ierr); call CHKERR(ierr)
      call DMSetDimension(dm3d, i3, ierr); call CHKERR(ierr)
      call DMPlexSetChart(dm3d, i0, chartsize, ierr); call CHKERR(ierr)

      call set_connectivity(dm2d, dm3d)

      call set_sf_graph(dm2d, dm3d)

      call set_coords(dm2d, dm3d)

    contains
      subroutine set_sf_graph(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        integer(mpiint) :: myid, ierr

        type(tDM) :: dmsf2d
        type(tPetscSF) :: sf2d, sf3d

        type(tPetscSection) :: section_2d_to_3d
        type(tVec) :: lVec, gVec
        real(ireals), pointer :: xv(:)

        integer(iintegers), pointer :: myidx(:) ! list of my indices that we do not own
        type(PetscSFNode), pointer  :: remote(:) ! rank and remote idx of those points
        integer(iintegers) :: nroots2d   ! number of root vertices on the current process (possible targets for leaves)
        integer(iintegers) :: nleaves2d  ! number of leaf vertices on the current process, references roots on any process
        integer(iintegers) :: nroots3d, nleaves3d
        integer(iintegers),allocatable :: ilocal_elements(:)
        type(PetscSFNode),allocatable :: iremote_elements(:)
        type(PetscCopyMode),parameter :: localmode=PETSC_COPY_VALUES, remotemode=PETSC_COPY_VALUES

        integer(iintegers) :: i, k, voff, ileaf, owner
        integer(iintegers), parameter :: idx_offset=0

        call DMClone(dm2d, dmsf2d, ierr); call CHKERR(ierr)

        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

        call DMGetPointSF(dmsf2d, sf2d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(sf2d, PETSC_NULL_SF, "-show_plex_sf2d", ierr); call CHKERR(ierr)

        call PetscSFGetGraph(sf2d, nroots2d, nleaves2d, myidx, remote, ierr); call CHKERR(ierr)

        call create_plex_section(comm, dmsf2d, 'plex_2d_to_3d_sf_graph_info', i1, &
          [i0], [i0], [ke1+ke], [ke1+ke], section_2d_to_3d)

        call DMSetDefaultSection(dmsf2d, section_2d_to_3d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(section_2d_to_3d, PETSC_NULL_SECTION, &
          '-show_dm2d_section_2d_to_3d', ierr); call CHKERR(ierr)

        call DMGetLocalVector(dmsf2d, lVec, ierr); call CHKERR(ierr)
        call DMGetGlobalVector(dmsf2d, gVec, ierr); call CHKERR(ierr)
        call VecSet(lVec, zero, ierr); call CHKERR(ierr)

        if(ldebug) then
          print *, myid, 'nroots', nroots2d, 'shape myidx', shape(myidx)
          print *, myid, 'nleaves', nleaves2d, 'shape remote', shape(remote)
          print *, myid, 'myidx', myidx
          print *, myid, 'remote', remote
        endif

        ! Distribute Info for 2D edges,
        call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        do i = e2dStart, e2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.lt.zero) then ! only add my local idx number if it belongs to me
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              xv(i1+voff+k) = iedge_top_icon_2_plex(i, k)
            enddo
            do k = 0, ke-1
              xv(i1+ke1+voff+k) = iface_side_icon_2_plex(i, k)
            enddo
          endif
        enddo
        do i = v2dStart, v2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.lt.zero) then
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              xv(i1+voff+k) = ivertex_icon_2_plex(i, k)
            enddo
            do k = 0, ke-1
              xv(i1+ke1+voff+k) = iedge_side_icon_2_plex(i, k)
            enddo
          endif
        enddo
        call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)

        call VecSet(gVec, zero, ierr); call CHKERR(ierr)
        call DMLocalToGlobalBegin(dmsf2d, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd(dmsf2d, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(dmsf2d, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(dmsf2d, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)

        nroots3d = chartsize
        nleaves3d = nleaves2d * (ke+ke1)
        allocate(ilocal_elements(nleaves3d))
        allocate(iremote_elements(nleaves3d))

        ileaf = 1
        call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        do i = e2dStart, e2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.ge.zero) then ! this is owned by someone else
            owner = remote(i1+voff)%rank
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              ilocal_elements(ileaf) = iedge_top_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = xv(i1+voff+k)
              if(ldebug) print *,myid,' 2dEdge top', i,'::', k,' local edge index', iedge_top_icon_2_plex(i, k), &
                'remote idx', xv(i1+voff+k)
              ileaf = ileaf+1
            enddo
            do k = 0, ke-1
              ilocal_elements(ileaf) = iface_side_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = xv(i1+ke1+voff+k)
              if(ldebug) print *,myid,' 2dEdge fac', i,'::', k,' local face index', ilocal_elements(ileaf), &
                'remote idx', xv(i1+ke1+voff+k)
              ileaf = ileaf+1
            enddo
          endif
        enddo
        do i = v2dStart, v2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.ge.zero) then ! this is owned by someone else
            owner = remote(i1+voff)%rank
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              ilocal_elements(ileaf) = ivertex_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = xv(i1+voff+k)
              if(ldebug) print *,myid,' 2dVert ver', i,'::', k,' local vert index', ilocal_elements(ileaf), &
                'remote idx', xv(i1+voff+k)
              ileaf = ileaf+1
            enddo
            do k = 0, ke-1
              ilocal_elements(ileaf) = iedge_side_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = xv(i1+ke1+voff+k)
              if(ldebug) print *,myid,' 2dVert edg', i,'::', k,' local face index', ilocal_elements(ileaf), &
                'remote idx', xv(i1+ke1+voff+k)
              ileaf = ileaf+1
            enddo
          endif
        enddo
        call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        if(ldebug) call CHKERR(int(ileaf-1-nleaves3d, mpiint), 'Does not add up... something is wrong')

        call PetscObjectViewFromOptions(lVec, PETSC_NULL_VEC, "-show_dm2d_sf_vec", ierr); call CHKERR(ierr)
        call DMRestoreLocalVector(dmsf2d, lVec, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(dmsf2d, gVec, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(section_2d_to_3d, ierr); call CHKERR(ierr)
        call DMDestroy(dmsf2d, ierr); call CHKERR(ierr)

        call DMGetPointSF(dm3d, sf3d, ierr); call CHKERR(ierr)
        call PetscSFSetGraph(sf3d, nroots3d, nleaves3d, ilocal_elements, localmode, &
          iremote_elements, remotemode, ierr); call CHKERR(ierr)
        call PetscSFSetUp(sf3d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(sf3d, PETSC_NULL_SF, "-show_plex_sf3d", ierr); call CHKERR(ierr)
      end subroutine

      subroutine set_connectivity(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        integer(iintegers) :: i, j, k
        integer(iintegers) :: icell, iface, iedge, ivert
        integer(iintegers), pointer :: cone(:)

        integer(iintegers) :: edge3(3), edge4(4), faces(5), vert2(2)

        ! Preallocation
        k = 0 ! k is running index for elements in DAG

        ! Every cell has 5 faces
        do i = 1, Ncells
          call DMPlexSetConeSize(dm3d, k, i5, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! top/bottom faces have 3 edges
        do i = 1, Nfaces2d * ke1
          call DMPlexSetConeSize(dm3d, k, i3, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! side faces have 4 edges
        do i = 1, Nedges2d * ke
          call DMPlexSetConeSize(dm3d, k, i4, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! Edges have 2 vertices
        do i = 1, Nedges
          call DMPlexSetConeSize(dm3d, k, i2, ierr); call CHKERR(ierr)
          k = k+1
        enddo
        call CHKERR(int((chartsize-Nverts)-k, mpiint), 'This does not add up, we forgot something?')

        call DMSetUp(dm3d, ierr); call CHKERR(ierr) ! Allocate space for cones

        ! Setup Connections
        ! First set five faces of cell
        do k = 0, ke-1
          do i = 0, Nfaces2d-1
            icell = icell_icon_2_plex(i, k)
            faces(1) = iface_top_icon_2_plex(i, k)
            faces(2) = iface_top_icon_2_plex(i, k+1)

            call DMPlexGetCone(dm2d, i, cone, ierr); call CHKERR(ierr) ! edges of face
            if(ldebug) print *,'iface2d', i, 'has edges', cone
            do j=1,size(cone)
              faces(2+j) = iface_side_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, i, cone, ierr); call CHKERR(ierr)

            call DMPlexSetCone(dm3d, icell, faces, ierr); call CHKERR(ierr)
            if(ldebug) print *,'icell',icell,'has faces:',faces
          enddo
        enddo

        ! set edges of top/bot faces
        do k = 0, ke1-1
          do i = 0, Nfaces2d-1
            iface = iface_top_icon_2_plex(i, k)

            call DMPlexGetCone(dm2d, i, cone, ierr); call CHKERR(ierr) ! edges of face
            do j=1,size(cone)
              edge3(j) = iedge_top_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, i, cone, ierr); call CHKERR(ierr)

            call DMPlexSetCone(dm3d, iface, edge3, ierr); call CHKERR(ierr)
            if(ldebug) print *,'iface', iface, 'gets edges', edge3
          enddo
        enddo

        ! set edges of vertical faces
        do k = 0, ke-1
          do i = 0, Nedges2d-1
            iedge = e2dStart+i
            iface = iface_side_icon_2_plex(iedge, k)
            edge4(1) = iedge_top_icon_2_plex(iedge, k)
            edge4(2) = iedge_top_icon_2_plex(iedge, k+1)

            call DMPlexGetCone(dm2d, iedge, cone, ierr); call CHKERR(ierr) ! vertices of 2d edge
            do j=1,size(cone)
              edge4(2+j) = iedge_side_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, iedge, cone, ierr); call CHKERR(ierr)
            if(ldebug) print *,'iface', iface, 'gets edges', edge4
            call DMPlexSetCone(dm3d, iface, edge4, ierr); call CHKERR(ierr)
          enddo
        enddo

        ! and then set the two vertices of edges in each level, i.e. vertices for edges in horizontal plane
        do k = 0, ke1-1
          do i = 0, Nedges2d-1
            iedge = e2dStart+i
            call DMPlexGetCone(dm2d, iedge, cone, ierr); call CHKERR(ierr) ! vertices of 2d edge
            do j=1,size(cone)
              vert2(j) = ivertex_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, iedge, cone, ierr); call CHKERR(ierr)
            if(ldebug) print *,'iedge', iedge_top_icon_2_plex(iedge, k), 'gets verts', vert2
            call DMPlexSetCone(dm3d, iedge_top_icon_2_plex(iedge, k), vert2, ierr); call CHKERR(ierr)
          enddo
        enddo

        ! and then set the two vertices of edges in each layer, i.e. vertices at the end of vertical edges
        do k = 0, ke-1
          do i = 0, Nverts2d-1
            ivert = v2dStart+i
            iedge = iedge_side_icon_2_plex(ivert, k)
            vert2(1) = ivertex_icon_2_plex(ivert, k)
            vert2(2) = ivertex_icon_2_plex(ivert, k+1)

            if(ldebug) print *,'edge', iedge, 'gets verts', vert2
            call DMPlexSetCone(dm3d, iedge, vert2, ierr); call CHKERR(ierr)
          enddo
        enddo

        call DMPlexSymmetrize(dm3d, ierr); call CHKERR(ierr)
        call DMPlexStratify(dm3d, ierr); call CHKERR(ierr)
      end subroutine

      subroutine set_coords(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        real(ireals), pointer :: coords2d(:), coords3d(:)
        type(tVec)            :: vec_coord2d, vec_coord3d
        type(tPetscSection)   :: coordSection2d, coordSection3d
        integer(iintegers)    :: coordsize, voff2d, voff3d
        integer(iintegers)    :: v2dStart, v2dEnd
        integer(iintegers)    :: v3dStart, v3dEnd

        real(ireals) :: distance, inv_distance
        integer(mpiint) :: ierr
        integer(iintegers) :: i, k, ivertex

        call DMPlexGetDepthStratum (dm3d, i0, v3dStart, v3dEnd, ierr); call CHKERR(ierr) ! 3D vertices
        call DMPlexGetDepthStratum (dm2d, i0, v2dStart, v2dEnd, ierr); call CHKERR(ierr) ! 2D vertices

        ! Create Coordinate stuff for 3D DM
        call DMGetCoordinateSection(dm3d, coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionSetNumFields(coordSection3d, i1, ierr); call CHKERR(ierr) ! just 1 field for spherical coords
        call PetscSectionSetUp(coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionSetFieldComponents(coordSection3d, i0, i3, ierr); call CHKERR(ierr)

        call PetscSectionSetChart(coordSection3d, v3dStart, v3dEnd, ierr);call CHKERR(ierr)
        do i = v3dStart, v3dEnd-1
          call PetscSectionSetDof(coordSection3d, i, i3, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldDof(coordSection3d, i, i0, i3, ierr); call CHKERR(ierr)
        enddo
        call PetscSectionSetUp(coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionGetStorageSize(coordSection3d, coordSize, ierr); call CHKERR(ierr)

        ! Create New Vec to hold 3D coordinates
        call VecCreate(PETSC_COMM_SELF, vec_coord3d, ierr); call CHKERR(ierr)
        call VecSetSizes(vec_coord3d, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
        call VecSetBlockSize(vec_coord3d, i3, ierr);call CHKERR(ierr)
        call VecSetType(vec_coord3d, VECSTANDARD, ierr);call CHKERR(ierr)

        call PetscObjectSetName(vec_coord3d, "coordinates", ierr); call CHKERR(ierr)

        ! Fill 3D Coord Vec
        call VecGetArrayF90(vec_coord3d, coords3d, ierr); call CHKERR(ierr)

        ! Get Coordinates from 2D DMPLEX
        call DMGetCoordinateSection(dm2d, coordSection2d, ierr); call CHKERR(ierr)
        call DMGetCoordinatesLocal(dm2d, vec_coord2d, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(vec_coord2d, coords2d, ierr); call CHKERR(ierr)

        do i = v2dStart, v2dEnd-1
          call PetscSectionGetOffset(coordSection2d, i, voff2d, ierr); call CHKERR(ierr)
          distance = norm(coords2d(voff2d+i1 : voff2d+i3))
          inv_distance = one / distance
          do k = 0, ke1-1
            ivertex = ivertex_icon_2_plex(i, k)
            call PetscSectionGetOffset(coordSection3d, ivertex, voff3d, ierr); call CHKERR(ierr)
            coords3d(voff3d+1:voff3d+3) = coords2d(voff2d+1:voff2d+3) * (distance + hhl(k+1)) * inv_distance
            if(ldebug) print *,'Setting coord for 2d', i, '3d vert', ivertex,':', &
              coords2d(voff2d+1:voff2d+3), '=>', &
              coords3d(voff3d+1:voff3d+3), &
              '(',(distance + hhl(k+1)) * inv_distance,')'
          enddo
        enddo

        call VecRestoreArrayReadF90(vec_coord2d, coords2d, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(vec_coord3d, coords3d, ierr); call CHKERR(ierr)

        call DMSetCoordinatesLocal(dm3d, vec_coord3d, ierr);call CHKERR(ierr)
        call PetscObjectViewFromOptions(vec_coord2d, PETSC_NULL_VEC, "-show_plex_coordinates2d", ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(vec_coord3d, PETSC_NULL_VEC, "-show_plex_coordinates3d", ierr); call CHKERR(ierr)
        call VecDestroy(vec_coord3d, ierr);call CHKERR(ierr)
      end subroutine

      !> @brief return the dmplex cell index for an  2d grid face index
      function icell_icon_2_plex(iface, k)
        integer(iintegers),intent(in) :: iface    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        icell_icon_2_plex = k + iface*ke
      end function

      !> @brief return the dmplex face index for an 2d face index situated at the top of a cell
      function iface_top_icon_2_plex(iface, k)
        integer(iintegers),intent(in) :: iface    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iface_top_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells
        iface_top_icon_2_plex = offset + k + iface*ke1
      end function

      !> @brief return the dmplex face index for an 2d edge index situated at the side of a cell, i.e. below a certain edge
      function iface_side_icon_2_plex(iedge, k)
        integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iface_side_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces2d*ke1
        iface_side_icon_2_plex = offset + k + (iedge-e2dStart)*ke
      end function

      !> @brief return the dmplex edge index for a given 2d edge index, i.e. the edges on the top/bot faces of cells
      function iedge_top_icon_2_plex(iedge, k)
        integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iedge_top_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces
        iedge_top_icon_2_plex = offset + k + (iedge-e2dStart)*ke1
      end function

      !> @brief return the dmplex edge index for a given 2d vertex index, i.e. the edges on the side faces of cells
      function iedge_side_icon_2_plex(ivertex, k)
        integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iedge_side_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces + Nedges2d*ke1
        iedge_side_icon_2_plex = offset + k + (ivertex-v2dStart)*ke
      end function

      !> @brief return the dmplex vertex index for a given 2d vertex index
      function ivertex_icon_2_plex(ivertex, k)
        integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: ivertex_icon_2_plex !< @param[out] icell_icon_2_plex, the vertex index in the dmplex, starts from plex%vStart and goes to plex%vEnd
        integer(iintegers) :: offset
        offset    = Ncells + Nfaces + Nedges
        ivertex_icon_2_plex = offset + k + (ivertex-v2dStart)*ke1
      end function
    end subroutine

    subroutine dump_ownership(dm, cmd_string)
      type(tDM), intent(in) :: dm
      character(len=*), intent(in) :: cmd_string

      type(tPetscSF) :: sf
      integer(iintegers) :: nroots, nleaves
      integer(iintegers), pointer :: myidx(:) ! list of my indices that we do not own
      type(PetscSFNode), pointer  :: remote(:) ! rank and remote idx of those points

      type(tDM) :: owner_dm
      type(tPetscSection) :: sec
      type(tVec) :: gVec
      real(ireals), pointer :: xv(:)
      integer(iintegers), pointer :: faces_of_cell(:)

      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: i, icell, idx, voff

      integer(mpiint) :: comm, myid, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call DMClone(dm, owner_dm, ierr); call CHKERR(ierr)

      call DMGetPointSF(owner_dm, sf, ierr); call CHKERR(ierr)
      call PetscSFGetGraph(sf, nroots, nleaves, myidx, remote, ierr); call CHKERR(ierr)

      call create_plex_section(comm, owner_dm, 'dmplex_ownership info', i1, &
        [i1], [i0], [i0], [i0], sec)
      call DMSetDefaultSection(owner_dm, sec, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(owner_dm, gVec, ierr); call CHKERR(ierr)

      ! Dump Cell Ownership
      call PetscObjectSetName(gVec, 'ownership_cells', ierr);call CHKERR(ierr)
      call VecGetArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      xv(:) = myid*one
      call VecRestoreArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gVec, PETSC_NULL_VEC, cmd_string, ierr); call CHKERR(ierr)

      ! Dump Face Ownership
      call PetscObjectSetName(gVec, 'ownership_non_local_faces', ierr);call CHKERR(ierr)

      call DMPlexGetHeightStratum(owner_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call VecGetArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      do icell = cStart, cEnd-1
        call PetscSectionGetOffset(sec, icell, voff, ierr); call CHKERR(ierr)
        xv(i1+voff) = zero
        call DMPlexGetCone(owner_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
        do i=1,size(faces_of_cell)
          call PetscFindInt(faces_of_cell(i), nleaves, myidx, idx, ierr); call CHKERR(ierr)
          if(idx.ge.i0) then ! not a local element
            xv(i1+voff) = xv(i1+voff) + i1
          endif
        enddo
      enddo
      call VecRestoreArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gVec, PETSC_NULL_VEC, cmd_string, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(owner_dm, gVec, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(sec, ierr); call CHKERR(ierr)
      call DMDestroy(owner_dm, ierr); call CHKERR(ierr)
    end subroutine

    ! Create a 2D Torus grid with Nx vertices horizontally and Ny rows of Vertices vertically
    subroutine create_2d_fish_plex(dm, Nx, Ny, lcyclic)
      type(tDM) :: dm, dmdist
      integer(iintegers), intent(in) :: Nx, Ny
      logical, intent(in) :: lcyclic

      integer(iintegers) :: chartsize, Nfaces, Nedges, Nvertices

      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      integer(mpiint) :: myid, numnodes, ierr

      logical, parameter :: ldebug=.False.

      call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(PETSC_COMM_WORLD, numnodes, ierr); call CHKERR(ierr)

      call DMPlexCreate(PETSC_COMM_WORLD, dm, ierr); call CHKERR(ierr)
      call PetscObjectSetName(dm, 'testplex Nx'//itoa(Nx)//'_Ny'//itoa(Ny), ierr); call CHKERR(ierr)
      call DMSetDimension(dm, i2, ierr); call CHKERR(ierr)

      if(modulo(Nx,i2).ne.0) call CHKERR(1_mpiint, 'Nx has to be even, e.g. 2,4,6...')
      if(modulo(Ny,i2).eq.0) call CHKERR(1_mpiint, 'Nx has to be uneven, e.g. 3,5...')

      if(myid.eq.0) then
        if(lcyclic) then
          Nfaces = (Nx-1) * (Ny-1) * 2        ! Number of faces
          Nedges = (Nx-1) * (Ny-1)            ! Number of edges on full height (y = 0, ...)
          Nedges = Nedges + (Nx*2-2) * (Ny-1) ! Number of edges on half height (y = {} + .5)
          Nvertices = (Nx-1) * (Ny-1)         ! Number of vertices
        else
          Nfaces = (Nx-1) * (Ny-1) * 2        ! Number of faces
          Nedges = (Nx-1) * Ny                ! Number of edges on full height (y = 0, ...)
          Nedges = Nedges + (Nx*2-1) * (Ny-1) ! Number of edges on half height (y = {} + .5)
          Nvertices = Nx * Ny                 ! Number of vertices
        endif
      else
        Nfaces = 0
        Nedges = 0
        Nvertices = 0
      endif
      if(ldebug) then
        print *, myid, 'Nfaces', Nfaces
        print *, myid, 'Nedges', Nedges
        print *, myid, 'Nverts', Nvertices
      endif

      chartsize = Nfaces + Nedges + Nvertices
      call DMPlexSetChart(dm, i0, chartsize, ierr); call CHKERR(ierr)

      call set_wedge_connectivity(dm, Nx, Ny, Nfaces, Nedges, lcyclic)

      call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
      call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

      call set_coords_serial(dm, Nx, lcyclic)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      if(myid.eq.0.and.ldebug) then
        print *,'pStart,End :: ',pStart, pEnd
        print *,'fStart,End :: ',fStart, fEnd
        print *,'eStart,End :: ',eStart, eEnd
        print *,'vStart,End :: ',vStart, vEnd
      endif

      ! True True is edges and vertices
      ! True False
      call DMPlexSetAdjacencyUseCone(dm, PETSC_TRUE, ierr); call CHKERR(ierr)
      call DMPlexSetAdjacencyUseClosure(dm, PETSC_FALSE, ierr); call CHKERR(ierr)
      !call DMPlexSetAdjacencyUseClosure(dm, PETSC_TRUE, ierr); call CHKERR(ierr)

      call DMPlexDistribute(dm, i0, PETSC_NULL_SF, dmdist, ierr); call CHKERR(ierr)
      if(dmdist.ne.PETSC_NULL_DM) then
        call DMDestroy(dm, ierr); call CHKERR(ierr)
        dm   = dmdist
      endif

      contains

        subroutine set_wedge_connectivity(dm, Nx, Ny, Nfaces, Nedges, lcyclic)
          type(tDM) :: dm
          integer(iintegers) :: Nx, Ny, Nfaces, Nedges
          logical, intent(in) :: lcyclic

          integer(iintegers) :: k, i, j, ioff, cone3(3), cone2(2)
          ! Preallocation
          ! Faces have 3 edges
          ioff = 0
          do k = 1, Nfaces
            call DMPlexSetConeSize(dm, ioff, i3, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
          ! Edges have 2 vertices
          do k = 1, Nedges
            call DMPlexSetConeSize(dm, ioff, i2, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo

          call DMSetUp(dm, ierr); call CHKERR(ierr) ! Allocate space for cones

          ioff = 0
          do k = 1, Nfaces
            j = (k-1) / (2*(Nx-1)) ! row of faces
            i = k-1 - j*(2*(Nx-1)) ! col of faces
            !print *,'faces',k,': i,j', i, j

            ! determine edges of a face
            if(modulo(i+modulo(j,i2),i2).eq.0) then ! this has a bot edge
              if(lcyclic) then
                cone3(1) = Nfaces + j*(Nx-1) + j*(Nx*2-2) + i/2  ! Nfaces offset + number of edges of full height + number of edges of half heights + i offset
                cone3(2) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2) + i  ! left edge  ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2) + i +1 ! right edge ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                if(i.eq.(Nx-1)*2-1) cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2)
              else
                cone3(1) = Nfaces + j*(Nx-1) + j*(Nx*2-1) + i/2  ! Nfaces offset + number of edges of full height + number of edges of half heights + i offset
                cone3(2) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-1) + i  ! left edge  ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-1) + i +1 ! right edge ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
              endif
              if(ldebug) print *,'upward edge of face', k-1, ':', cone3
            else
              if(lcyclic) then
                cone3(1) = Nfaces + (j+1)*(Nx-1) + (j+1)*(Nx*2-2) + i/2
                cone3(2) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2) + i
                cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2) + i +1
                if(i.eq.(Nx-1)*2-1) cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2)
                if(j.eq.Ny-2) cone3(1) = Nfaces + i/2
              else
                cone3(1) = Nfaces + (j+1)*(Nx-1) + (j+1)*(Nx*2-1) + i/2
                cone3(2) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-1) + i
                cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-1) + i +1
              endif
              if(ldebug) print *,'                                                                            downward edge of face', k-1, ':', i, j,':', cone3
            endif

            call DMPlexSetCone(dm,  ioff, cone3, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
          ! Edges have 2 vertices
          do k = 1, Nedges
            if(lcyclic) then
              j = (k-1) / (Nx-1 + Nx*2-2) ! row of edge ! goes up to Ny
              i = k-1 - j*(Nx-1 + Nx*2-2) ! col of edge
              if(i.lt.Nx-1) then ! bottom edge
                cone2(1) = Nfaces + Nedges + j*(Nx-1) + i
                cone2(2) = Nfaces + Nedges + j*(Nx-1) + i +1
                if(i.eq.Nx-2) cone2(2) = Nfaces + Nedges + j*(Nx-1)
                if(ldebug) print *,'bot edge',k,ioff,': i,j', i, j, ':', cone2
              else ! sideward edge
                if(modulo(i+j,i2).ne.0) then ! slash
                  cone2(1) = Nfaces + Nedges + j*(Nx-1) + (i-Nx+1)/2
                  cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1) + (i-Nx+1)/2 + modulo(j,i2)
                  if(i.eq.(Nx-1+Nx*2-2-1)) cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1)
                  if(j.eq.Ny-2) then
                    cone2(2) = Nfaces + Nedges + (i-Nx+1)/2 + modulo(j,i2)
                    if(i.eq.(Nx-1+Nx*2-2-1)) cone2(2) = Nfaces + Nedges
                  endif
                  if(ldebug) print *,'                                                                            slash side edge',k,ioff,': i,j', i, j, ':', cone2
                else ! backslash
                  cone2(1) = Nfaces + Nedges + j*(Nx-1) + (i-Nx+1)/2 + modulo(j+1,i2)
                  cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1) + (i-Nx+1)/2
                  if(i.eq.(Nx-1+Nx*2-2-1)) cone2(1) = Nfaces + Nedges + j*(Nx-1)
                  if(j.eq.Ny-2) cone2(2) = Nfaces + Nedges + (i-Nx+1)/2
                  if(ldebug) print *,'                                                                            backs side edge',k,ioff,': i,j', i, j, ':', cone2
                endif
              endif
            else
              j = (k-1) / (Nx-1 + Nx*2-1) ! row of edge ! goes up to Ny
              i = k-1 - j*(Nx-1 + Nx*2-1) ! col of edge
              if(i.lt.Nx-1) then ! bottom edge
                cone2(1) = Nfaces + Nedges + j*Nx + i
                cone2(2) = Nfaces + Nedges + j*Nx + i +1
                if(ldebug) print *,'bot edge',k,ioff,': i,j', i, j, ':', cone2
              else ! sideward edge
                if(modulo(i+j,i2).ne.0) then ! slash
                  cone2(1) = Nfaces + Nedges + j*Nx + (i-Nx+1)/2
                  cone2(2) = Nfaces + Nedges + (j+1)*Nx + (i-Nx+1)/2 + modulo(j,i2)
                  !print *,'                                                                            slash side edge',k,ioff,': i,j', i, j, ':', cone2
                else ! backslash
                  cone2(1) = Nfaces + Nedges + j*Nx + (i-Nx+1)/2 + modulo(j+1,i2)
                  cone2(2) = Nfaces + Nedges + (j+1)*Nx + (i-Nx+1)/2
                  !print *,'                                                                            backs side edge',k,ioff,': i,j', i, j, ':', cone2
                endif
              endif
            endif
            call DMPlexSetCone(dm,  ioff, cone2, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
        end subroutine

        subroutine set_coords_serial(dm, Nx, lcyclic)
          type(tDM) :: dm
          integer(iintegers), intent(in) :: Nx
          logical, intent(in) :: lcyclic
          real(ireals), pointer:: coords(:)
          type(tVec)           :: coordinates
          integer(iintegers)   :: dimEmbed, coordSize, vStart, vEnd, pStart, pEnd, eStart, eEnd, voff, cStart, cEnd
          type(tPetscSection)  :: coordSection
          integer(iintegers)   :: iv, i, j

          real(ireals), parameter :: dx=1, dy=1
          real(ireals), parameter :: ds=sqrt(dy**2 - (dx/2)**2)
          real(ireals) :: x, y, z

          call DMGetCoordinateDim(dm, dimEmbed, ierr); CHKERRQ(ierr)
          dimEmbed = 3

          call DMGetCoordinateSection(dm, coordSection, ierr); CHKERRQ(ierr)

          call PetscSectionSetNumFields(coordSection, i1, ierr); CHKERRQ(ierr)
          call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
          call PetscSectionSetFieldComponents(coordSection, i0, dimEmbed, ierr); CHKERRQ(ierr)

          call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
          call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
          call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
          call DMPlexGetHeightStratum (dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! faces

          call PetscSectionSetChart(coordSection, vStart, vEnd, ierr);CHKERRQ(ierr)

          do i = vStart, vEnd-1
            call PetscSectionSetDof(coordSection, i, dimEmbed, ierr); CHKERRQ(ierr)
            call PetscSectionSetFieldDof(coordSection, i, i0, dimEmbed, ierr); CHKERRQ(ierr)
          enddo

          call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
          call PetscSectionGetStorageSize(coordSection, coordSize, ierr); CHKERRQ(ierr)

          call VecCreate(PETSC_COMM_SELF, coordinates, ierr); CHKERRQ(ierr)
          call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);CHKERRQ(ierr)
          call VecSetBlockSize(coordinates, dimEmbed, ierr);CHKERRQ(ierr)
          call VecSetType(coordinates, VECSTANDARD, ierr);CHKERRQ(ierr)

          call PetscObjectSetName(coordinates, "coordinates", ierr); CHKERRQ(ierr)

          call VecGetArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)

          do iv = vStart, vEnd-1
            if(lcyclic) then
              j = (iv-vStart) / (Nx-1)
              i = (iv-vStart) - j*(Nx-1)
              z = sqrt(x**2 + y**2)
            else
              j = (iv-vStart) / Nx
              i = (iv-vStart) - j*Nx
              z = 100
            endif
            x = (real(modulo(j,i2),ireals)*.5_ireals + real(i, ireals))*dx
            y = j*ds
            if(ldebug) print *,'iv',iv,':', i, j,'=>', x, y, z

            call PetscSectionGetOffset(coordSection, iv, voff, ierr); coords(voff+1:voff+3) = [x, y, z]
          enddo

          call VecRestoreArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)
          call DMSetCoordinatesLocal(dm, coordinates, ierr);CHKERRQ(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); CHKERRQ(ierr)
          call VecDestroy(coordinates, ierr);CHKERRQ(ierr)
        end subroutine
      end subroutine
end module
