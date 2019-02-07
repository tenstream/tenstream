module m_plex2rayli

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i8

  use m_helper_functions, only: CHKERR

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
      integer(iintegers), pointer :: faces_of_cell(:), cell_support(:), edges_of_face(:), verts_of_edge(:)
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
        !call DMPlexGetSupport(dm, i, cell_support, ierr); call CHKERR(ierr)

        call DMPlexGetConeSize(dm, i, Nedges, ierr); call CHKERR(ierr)
        call DMPlexGetCone(dm, i, edges_of_face, ierr); call CHKERR(ierr)

        if(Nedges.eq.3) then
          !call DMPlexSetSupport(dmrayli, k, cell_support, ierr); call CHKERR(ierr)
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

        call DMPlexRestoreSupport(dm, iface, cell_support, ierr); call CHKERR(ierr)
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
end module
