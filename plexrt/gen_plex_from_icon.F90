module m_gen_plex_from_icon

#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload
  use m_helper_functions, only: CHKERR
  use m_icon_grid, only: t_icongrid, read_icon_grid_file, decompose_icon_grid
  use m_plex_grid, only: t_plexgrid, load_plex_from_file, update_plex_indices, &
    TOP_BOT_FACE, SIDE_FACE,  icell_icon_2_plex, icell_plex_2_icon
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    default_str_len, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one,       &
    init_mpi_data_parameters

  implicit none

  logical, parameter :: ldebug=.True.

  integer(mpiint) :: ierr

  contains

    subroutine create_dmplex_2d(icongrid, dmname, plex)
      type(t_icongrid), intent(in) :: icongrid
      type(t_plexgrid), intent(out) :: plex
      character(len=*), intent(in) :: dmname
      integer(iintegers) :: chartsize, depth
      integer(iintegers) :: i, icell, iedge
      integer(iintegers) :: offset_edges, offset_vertices, vert2(2), edge3(3)

      allocate(plex%dm)
      call DMPlexCreate(PETSC_COMM_WORLD, plex%dm, ierr);call CHKERR(ierr)

      call PetscObjectSetName(plex%dm, trim(dmname), ierr);call CHKERR(ierr)
      call DMSetDimension(plex%dm, i2, ierr);call CHKERR(ierr)

      chartsize = size(icongrid%cell_index) + size(icongrid%edge_index) + size(icongrid%vertex_index)
      call DMPlexSetChart(plex%dm, i0, chartsize, ierr); call CHKERR(ierr)

      offset_edges = size(icongrid%cell_index)
      offset_vertices = size(icongrid%cell_index) + size(icongrid%edge_index)
      print *,'Chartsize:', chartsize, 'offsets: ', offset_edges, offset_vertices

      ! Preallocation
      ! Every cell has three edges
      do i = 0, offset_edges-1
        call DMPlexSetConeSize(plex%dm, i, 3, ierr); call CHKERR(ierr)
      enddo

      ! Edges have 2 vertices
      do i = offset_edges, offset_vertices-1
        call DMPlexSetConeSize(plex%dm, i, 2, ierr); call CHKERR(ierr)
      enddo

      call DMSetUp(plex%dm, ierr); call CHKERR(ierr) ! Allocate space for cones

      ! Setup Connections
      ! First set three edges of cell
      do i = 1, size(icongrid%cell_index)
        icell = icongrid%cell_index(i)
        edge3 = icongrid%edge_of_cell(icell,:)

        call DMPlexSetCone(plex%dm, icell-i1, edge3+offset_edges-i1, ierr); call CHKERR(ierr)
      enddo

      !! and then set the two vertices of edge
      do i = 1, size(icongrid%edge_index)
        iedge = icongrid%edge_index(i)
        vert2 = icongrid%edge_vertices(iedge,:)

        call DMPlexSetCone(plex%dm, iedge+offset_edges-i1, vert2+offset_vertices-i1, ierr); call CHKERR(ierr)
      enddo

      ! Set cone orientations -- could traverse the DAG other way round, with values being -1
      !do i = 0, size(icongrid%cell_index)-1
      !  call DMPlexSetConeOrientation(plex%dm, i, [i0,i0,i0], ierr); call CHKERR(ierr)
      !enddo
      !do i = offset_edges, offset_vertices-1
      !  call DMPlexSetConeOrientation(plex%dm, i, [i0,i0], ierr); call CHKERR(ierr)
      !enddo

      call DMPlexSymmetrize(plex%dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug) then
        call DMPlexGetDepth(plex%dm, depth, ierr); CHKERRQ(ierr)
        print *,'Depth of Stratum:', depth
      endif

      call DMPlexGetChart(plex%dm, plex%pStart, plex%pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(plex%dm, i0, plex%cStart, plex%cEnd, ierr); call CHKERR(ierr) ! cells
      !call DMPlexGetHeightStratum(plex%dm, i1, plex%fStart, plex%fEnd, ierr); call CHKERR(ierr) ! faces / edges
      call DMPlexGetDepthStratum (plex%dm, i1, plex%eStart, plex%eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetDepthStratum (plex%dm, i0, plex%vStart, plex%vEnd, ierr); call CHKERR(ierr) ! vertices

      if(ldebug) then
        print *, 'pstart', plex%pstart, 'pEnd', plex%pEnd
        print *, 'cStart', plex%cStart, 'cEnd', plex%cEnd
        print *, 'fStart', plex%fStart, 'fEnd', plex%fEnd
        print *, 'eStart', plex%eStart, 'eEnd', plex%eEnd
        print *, 'vStart', plex%vStart, 'vEnd', plex%vEnd
      endif

      call set_coords()
      call label_domain()

      contains
        subroutine label_domain()
          integer(iintegers) :: e, numcells
          integer(iintegers), pointer :: icells(:)
          !DMLabel :: l

          call DMCreateLabel(plex%dm, plex%boundary_label, ierr); call CHKERR(ierr)

          do e = plex%eStart, plex%eEnd-i1
            call DMPlexGetSupportSize(plex%dm, e, numcells, ierr); call CHKERR(ierr)
            if(numcells.eq.1) then
              call DMPlexGetSupport(plex%dm, e, icells, ierr); call CHKERR(ierr)
              call DMSetLabelValue(plex%dm, plex%boundary_label, icells(1), 1, ierr); call CHKERR(ierr)
            endif
          enddo

          !call dmgetlabel(plex%dm, plex%boundary_label, l, ierr); call CHKERR(ierr)
          !call dmlabelview(l, petsc_viewer_stdout_world, ierr); call CHKERR(ierr)

        end subroutine
        subroutine set_coords()
          PetscScalar, pointer :: coords(:)
          Vec                  :: coordinates
          PetscInt             :: dimEmbed, coordSize, voff, ind
          PetscSection         :: coordSection

          real(ireals) :: cart_coord(2)

          call DMGetCoordinateDim(plex%dm, dimEmbed, ierr); call CHKERR(ierr)
          print *,'dimEmbed = ', dimEmbed

          call DMGetCoordinateSection(plex%dm, coordSection, ierr); call CHKERR(ierr)

          call PetscSectionSetNumFields(coordSection, i1, ierr); call CHKERR(ierr)
          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldComponents(coordSection, i0, dimEmbed, ierr); call CHKERR(ierr)

          call PetscSectionSetChart(coordSection, plex%vStart, plex%vEnd, ierr);call CHKERR(ierr)

          do i = plex%vStart, plex%vEnd-i1
            call PetscSectionSetDof(coordSection, i, dimEmbed, ierr); call CHKERR(ierr)
            call PetscSectionSetFieldDof(coordSection, i, i0, dimEmbed, ierr); call CHKERR(ierr)
          enddo

          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscSectionGetStorageSize(coordSection, coordSize, ierr); call CHKERR(ierr)
          print *,'Coord Section has size:', coordSize

          call VecCreate(PETSC_COMM_SELF, coordinates, ierr); call CHKERR(ierr)
          call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
          call VecSetBlockSize(coordinates, dimEmbed, ierr);call CHKERR(ierr)
          call VecSetType(coordinates, VECSTANDARD, ierr);call CHKERR(ierr)

          call PetscObjectSetName(coordinates, "coordinates", ierr); call CHKERR(ierr)

          call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)
          print *,'bounds coords:', lbound(coords), ubound(coords)

          ! set vertices as coordinates
          do i = i1, size(icongrid%vertex_index)
            ind = plex%vStart + i - i1
            call PetscSectionGetOffset(coordSection, ind, voff, ierr); call CHKERR(ierr)

            !cart_coord = [icongrid%cartesian_x_vertices(i), icongrid%cartesian_y_vertices(i), icongrid%cartesian_z_vertices(i)]
            cart_coord = [icongrid%cartesian_x_vertices(i), icongrid%cartesian_y_vertices(i)]
            coords(voff+i1 : voff+dimEmbed) = cart_coord(i1:dimEmbed)
            !print *,'setting coords',cart_coord,'to',voff+i1 , voff+dimEmbed
          enddo

          !print *,'coords', shape(coords), '::', coords
          call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

          call DMSetCoordinatesLocal(plex%dm, coordinates, ierr);call CHKERR(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); call CHKERR(ierr)

          call VecDestroy(coordinates, ierr);call CHKERR(ierr)
        end subroutine

    end subroutine

   ! Nz is number of layers
    subroutine create_dmplex_3d(comm, icongrid, dmname, Nz, plex)
      integer(mpiint),intent(in) :: comm
      type(t_icongrid),intent(in) :: icongrid
      type(t_plexgrid),intent(out) :: plex
      character(len=*), intent(in) :: dmname
      integer(iintegers), intent(in) :: Nz
      integer(mpiint) :: myid
      integer(iintegers) :: chartsize, depth
      integer(iintegers) :: i, k, icell, ivertex, iedge
      integer(iintegers) :: Ncells, Nfaces, Nedges, Nvertices
      integer(iintegers) :: offset_faces, offset_edges, offset_vertices
      integer(iintegers) :: offset_faces_sides, offset_edges_vertical
      integer(iintegers) :: vert2(2), edge3(3), edge4(4), faces(5)

      type(tDMLabel) :: faceposlabel, zindexlabel, TOAlabel

      ! Create Plex
      allocate(plex%dm)
      plex%comm = comm
      call DMPlexCreate(plex%comm, plex%dm, ierr);call CHKERR(ierr)

      call PetscObjectSetName(plex%dm, trim(dmname), ierr);call CHKERR(ierr)
      call DMSetDimension(plex%dm, i3, ierr);call CHKERR(ierr)

      call mpi_comm_rank(plex%comm, myid, ierr); CHKERRQ(ierr)
      if(myid.ne.0) return

      plex%Nz = Nz

      Ncells    = icongrid%Nfaces * plex%Nz
      Nfaces    = icongrid%Nfaces * (plex%Nz+i1) + icongrid%Nedges * plex%Nz
      Nedges    = icongrid%Nedges * (plex%Nz+i1) + icongrid%Nvertices * plex%Nz
      Nvertices = icongrid%Nvertices * (plex%Nz+i1)

      chartsize = Ncells + Nfaces + Nedges + Nvertices
      call DMPlexSetChart(plex%dm, i0, chartsize, ierr); call CHKERR(ierr)

      offset_faces = Ncells
      offset_faces_sides = offset_faces + (plex%Nz+1)*icongrid%Nfaces
      offset_edges = Ncells + Nfaces
      offset_edges_vertical = offset_edges + icongrid%Nedges * (plex%Nz+i1)
      offset_vertices = Ncells + Nfaces + Nedges
      print *,'offsets faces:', offset_faces, offset_faces_sides
      print *,'offsets edges:', offset_edges, offset_edges_vertical
      print *,'offsets verte:', offset_vertices
      print *,'Chartsize:', chartsize

      ! Create some labels ... those are handy later when setting up matrices etc
      call DMCreateLabel(plex%dm, "Face Position", ierr); CHKERRQ(ierr)
      call DMCreateLabel(plex%dm, "Vertical Index", ierr); CHKERRQ(ierr)
      call DMCreateLabel(plex%dm, "TOA", ierr); CHKERRQ(ierr)

      call DMGetLabel(plex%dm, "Face Position", faceposlabel, ierr); CHKERRQ(ierr)
      call DMGetLabel(plex%dm, "Vertical Index", zindexlabel, ierr); CHKERRQ(ierr)
      call DMGetLabel(plex%dm, "TOA", TOAlabel, ierr); CHKERRQ(ierr)


      ! Preallocation
      ! Every cell has 5 faces
      do i = i0, offset_faces-i1
        call DMPlexSetConeSize(plex%dm, i, i5, ierr); call CHKERR(ierr)
      enddo

      ! top/bottom faces have 3 edges
      do i = offset_faces, offset_faces_sides-i1
        call DMPlexSetConeSize(plex%dm, i, i3, ierr); call CHKERR(ierr)
      enddo

      ! side faces have 4 edges
      do i = offset_faces_sides, offset_edges-i1
        call DMPlexSetConeSize(plex%dm, i, i4, ierr); call CHKERR(ierr)
      enddo

      ! Edges have 2 vertices
      do i = offset_edges, offset_vertices-i1
        call DMPlexSetConeSize(plex%dm, i, i2, ierr); call CHKERR(ierr)
      enddo

      call DMSetUp(plex%dm, ierr); call CHKERR(ierr) ! Allocate space for cones

      ! Setup Connections
      ! First set five faces of cell
      do k = 1, plex%Nz
        do i = 1, icongrid%Nfaces
          icell = icongrid%cell_index(i)
          edge3 = icongrid%edge_of_cell(icell,:)

          faces(1) = offset_faces + icongrid%Nfaces*(k-1) + icell-i1  ! top face
          faces(2) = offset_faces + icongrid%Nfaces*k + icell-i1      ! bot face

          faces(3:5) = offset_faces_sides + (edge3-i1) + (k-1)*icongrid%Nedges

          call DMPlexSetCone(plex%dm, icell_icon_2_plex(icongrid, plex, icell, k), faces, ierr); call CHKERR(ierr)

          call DMLabelSetValue(zindexlabel, icell, k, ierr); call CHKERR(ierr)
        enddo
      enddo

      ! set edges of top/bot faces
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nfaces
          icell = icongrid%cell_index(i)
          edge3 = offset_edges + icongrid%edge_of_cell(icell,:)-i1 + icongrid%Nedges*(k-i1)

          call DMPlexSetCone(plex%dm, offset_faces + icongrid%Nfaces*(k-1) + icell -i1, edge3, ierr); call CHKERR(ierr)
          !print *,'edges @ horizontal faces',icell,':',edge3,'petsc:', offset_faces + Nfaces*(k-1) + icell -i1, '::', edge3
          call DMLabelSetValue(faceposlabel, offset_faces + icongrid%Nfaces*(k-1) + icell -i1, TOP_BOT_FACE, ierr); call CHKERR(ierr)
          if (k.eq.i1) then
            call DMLabelSetValue(TOAlabel, offset_faces + icongrid%Nfaces*(k-1) + icell -i1, i1, ierr); call CHKERR(ierr)
          endif
        enddo
      enddo

      ! set edges of vertical faces
      do k = 1, plex%Nz ! layers
        do i = 1, icongrid%Nedges
          iedge = icongrid%edge_index(i)
          edge4(1) = offset_edges + iedge-i1 + icongrid%Nedges*(k-i1)
          edge4(2) = offset_edges + iedge-i1 + icongrid%Nedges*(k)

          vert2 = icongrid%edge_vertices(iedge,:)-i1 + icongrid%Nvertices*(k-i1)
          edge4(3:4) = offset_edges_vertical + vert2

          call DMPlexSetCone(plex%dm, offset_faces_sides + iedge-i1 + icongrid%Nedges*(k-i1),edge4, ierr); call CHKERR(ierr)
          !print *,'edges @ vertical faces',iedge,':',edge4,'petsc:', offset_faces_sides + iedge-i1 + Nedges*(k-i1), '::', edge4
          call DMLabelSetValue(faceposlabel, offset_faces_sides + iedge-i1 + icongrid%Nedges*(k-i1), SIDE_FACE, ierr); call CHKERR(ierr)
        enddo
      enddo


      ! and then set the two vertices of edges in each level
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nedges
          iedge = icongrid%edge_index(i)
          vert2 = offset_vertices + icongrid%edge_vertices(iedge,:) + icongrid%Nvertices*(k-i1)

          call DMPlexSetCone(plex%dm, offset_edges + iedge-i1 + icongrid%Nedges*(k-i1), vert2-i1, ierr); call CHKERR(ierr)
          !print *,'vertices @ edge2d',iedge,':',vert2,'petsc:',offset_edges + iedge-i1 + icongrid%Nedges*(k-i1),'::', vert2-i1
        enddo
      enddo

      ! and then set the two vertices of edges in each layer
      do k = 1, plex%Nz ! layer
        do i = 1, icongrid%Nvertices
          ivertex = icongrid%vertex_index(i)
          vert2(1) = offset_vertices + ivertex + icongrid%Nvertices*(k-i1)
          vert2(2) = offset_vertices + ivertex + icongrid%Nvertices*(k)

          call DMPlexSetCone(plex%dm, offset_edges_vertical + ivertex-i1 + icongrid%Nvertices*(k-i1), vert2 - i1, ierr); call CHKERR(ierr)
          !print *,'vertices @ edge_vertical',ivertex,':',vert2,'petsc:',offset_edges_vertical + ivertex-i1 + Nvertices*(k-i1),':', vert2 - i1
        enddo
      enddo

      ! Set cone orientations -- could traverse the DAG other way round, with values being -1
      !do i = 0, size(icongrid%cell_index)-1
      !  call DMPlexSetConeOrientation(plex%dm, i, [i0,i0,i0], ierr); call CHKERR(ierr)
      !enddo
      !do i = offset_edges, offset_vertices-1
      !  call DMPlexSetConeOrientation(plex%dm, i, [i0,i0], ierr); call CHKERR(ierr)
      !enddo

      print *,'Symmetrize'
      call DMPlexSymmetrize(plex%dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug) then
        call DMPlexGetDepth(plex%dm, depth, ierr); CHKERRQ(ierr)
        print *,'Depth of Stratum:', depth
      endif

      call update_plex_indices(plex)

      call set_coords()
      !call label_domain()

      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
      contains
        subroutine label_domain()
          integer(iintegers) :: e, numcells
          integer(iintegers), pointer :: icells(:)
          !DMLabel :: l

          call DMCreateLabel(plex%dm, plex%boundary_label, ierr); call CHKERR(ierr)

          do e = plex%eStart, plex%eEnd-i1
            call DMPlexGetSupportSize(plex%dm, e, numcells, ierr); call CHKERR(ierr)
            if(numcells.eq.1) then
              call DMPlexGetSupport(plex%dm, e, icells, ierr); call CHKERR(ierr)
              call DMSetLabelValue(plex%dm, plex%boundary_label, icells(1), 1, ierr); call CHKERR(ierr)
            endif
          enddo

          !call dmgetlabel(plex%dm, plex%boundary_label, l, ierr); call CHKERR(ierr)
          !call dmlabelview(l, petsc_viewer_stdout_world, ierr); call CHKERR(ierr)

        end subroutine
        subroutine set_coords()
          PetscScalar, pointer :: coords(:)
          Vec                  :: coordinates
          PetscInt             :: dimEmbed, coordSize, voff, ind
          PetscSection         :: coordSection

          real(ireals) :: cart_coord(3)

          real(ireals), parameter :: sphere_radius=6371229  ! [m]
          logical :: l_is_spherical_coords

          l_is_spherical_coords = any(icongrid%cartesian_z_vertices.ne.0)

          call DMGetCoordinateDim(plex%dm, dimEmbed, ierr); call CHKERR(ierr)
          print *,'dimEmbed = ', dimEmbed

          call DMGetCoordinateSection(plex%dm, coordSection, ierr); call CHKERR(ierr)

          call PetscSectionSetNumFields(coordSection, i1, ierr); call CHKERR(ierr)
          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldComponents(coordSection, i0, dimEmbed, ierr); call CHKERR(ierr)

          call PetscSectionSetChart(coordSection, plex%vStart, plex%vEnd, ierr);call CHKERR(ierr)

          do i = plex%vStart, plex%vEnd-i1
            call PetscSectionSetDof(coordSection, i, dimEmbed, ierr); call CHKERR(ierr)
            call PetscSectionSetFieldDof(coordSection, i, i0, dimEmbed, ierr); call CHKERR(ierr)
          enddo

          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(coordSection, PETSC_NULL_SECTION, "-show_coordinates_section", ierr); call CHKERR(ierr)
          call PetscSectionGetStorageSize(coordSection, coordSize, ierr); call CHKERR(ierr)
          print *,'Coord Section has size:', coordSize

          call VecCreate(PETSC_COMM_SELF, coordinates, ierr); call CHKERR(ierr)
          call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
          call VecSetBlockSize(coordinates, dimEmbed, ierr);call CHKERR(ierr)
          call VecSetType(coordinates, VECSTANDARD, ierr);call CHKERR(ierr)

          call PetscObjectSetName(coordinates, "coordinates", ierr); call CHKERR(ierr)

          call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)
          print *,'bounds coords:', lbound(coords), ubound(coords)

          ! set vertices as coordinates
          do k = 1, plex%Nz+1
            do i = 1, icongrid%Nvertices
              ind = plex%vStart + i - i1 + icongrid%Nvertices*(k-i1)
              call PetscSectionGetOffset(coordSection, ind, voff, ierr); call CHKERR(ierr)

              cart_coord = [icongrid%cartesian_x_vertices(i), icongrid%cartesian_y_vertices(i), &
                            icongrid%cartesian_z_vertices(i)]

              if(l_is_spherical_coords) then
                cart_coord = cart_coord * (sphere_radius + (plex%Nz-k)*200)
              else
                cart_coord(3) = (plex%Nz+1-k)*200
              endif
              coords(voff+i1 : voff+dimEmbed) = cart_coord(i1:dimEmbed)
            enddo
          enddo

          print *,'coords', shape(coords), '::', coords
          call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

          call DMSetCoordinatesLocal(plex%dm, coordinates, ierr);call CHKERR(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); call CHKERR(ierr)

          call VecDestroy(coordinates, ierr);call CHKERR(ierr)
        end subroutine

    end subroutine

    subroutine create_mass_vec(icongrid, plex)
      type(t_icongrid), intent(in) :: icongrid
      type(t_plexgrid), intent(in) :: plex
      PetscInt    :: i, depth, section_size, vecsize, voff
      PetscSection :: s

      Vec :: globalVec
      PetscScalar, pointer :: xv(:)

      integer(iintegers) :: icell_k(2)
      integer(iintegers) :: cell_ownership(plex%cEnd-plex%cStart)
      integer(iintegers) :: edge_ownership(plex%eEnd-plex%eStart)
      integer(iintegers) :: vertex_ownership(plex%vEnd-plex%vStart)

      call DMPlexGetDepth(plex%dm, depth, ierr); CHKERRQ(ierr)
      print *,'Depth of Stratum:', depth

      ! Create Default Section
      call PetscSectionCreate(PETSC_COMM_WORLD, s, ierr); CHKERRQ(ierr)
      call PetscSectionSetNumFields(s, 1, ierr); CHKERRQ(ierr)
      call PetscSectionSetChart(s, plex%cStart, plex%cEnd, ierr); CHKERRQ(ierr)

      do i = plex%cStart, plex%cEnd-1
        call PetscSectionSetDof(s, i, 1, ierr); CHKERRQ(ierr)
        call PetscSectionSetFieldDof(s, i, 0, 1, ierr); CHKERRQ(ierr)
      enddo

      call PetscSectionSetUp(s, ierr); CHKERRQ(ierr)

      call PetscSectionGetStorageSize(s, section_size, ierr);
      print *,'Section has size:', section_size

      call DMSetDefaultSection(plex%dm, s, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(s, 'CellSection', ierr);CHKERRQ(ierr)

      ! Now lets get vectors!
      call DMGetGlobalVector(plex%dm, globalVec,ierr); CHKERRQ(ierr)

      call VecGetSize(globalVec,vecsize, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(globalVec, 'massVec', ierr);CHKERRQ(ierr)

      call decompose_icon_grid(icongrid, 2, cell_ownership, edge_ownership, vertex_ownership)

      call VecGetArrayF90(globalVec, xv, ierr); CHKERRQ(ierr)
      do i = plex%cStart, plex%cEnd-1
        call PetscSectionGetOffset(s, i, voff, ierr); call CHKERR(ierr)
        icell_k = icell_plex_2_icon(icongrid, plex, i)
        !xv(voff+i1) = icell_k(2)*plex%Nfaces + icell_k(1)
        if(cell_ownership(icell_k(1)).eq.0) then
          xv(voff+i1) = icell_k(1)
        else
          xv(voff+i1) = -icell_k(1) !+ cell_ownership(icell_k(1))*1000
        endif
      enddo
      call VecRestoreArrayF90(globalVec, xv, ierr); CHKERRQ(ierr)

      call PetscObjectViewFromOptions(globalVec, PETSC_NULL_VEC, '-show_massVec', ierr); CHKERRQ(ierr)

      call DMRestoreGlobalVector(plex%dm, globalVec, ierr); CHKERRQ(ierr)

      !do i = plex%cStart, plex%cEnd-1
      !  call DMPlexGetCone(plex%dm, i, cone, ierr); CHKERRQ(ierr)
      !  print *,'Cone',i,'::',cone
      !  call DMPlexRestoreCone(plex%dm, i, cone, ierr); CHKERRQ(ierr)
      !enddo

      call PetscSectionDestroy(s, ierr); CHKERRQ(ierr)
    end subroutine

    !> @brief distribute plex mesh and info about the  icon base grid
    subroutine distribute_dmplex(comm, plex)
      integer(mpiint) :: comm                      !< @param[in] Global MPI Communicator
      type(t_plexgrid), intent(inout) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(mpiint) :: numnodes, myid, mpierr
      PetscSF         :: pointSF
      DM              :: dmdist

      call mpi_comm_rank( comm, myid, mpierr)
      call mpi_comm_size( comm, numnodes, mpierr)

      if(.not.(allocated(plex%dm))) then
        if(myid.eq.0) stop 'distribute_dmplex :: rank 0 has to have an already instantiated dmplex!'
        allocate(plex%dm)
        call DMPlexCreate(comm, plex%dm, ierr);call CHKERR(ierr)
      endif

      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
      if (numnodes.gt.1) then
        if(ldebug) then
          print *,'Distributing DMPlex to ',numnodes,' ranks'
        endif
        call DMPlexDistribute(plex%dm, i0, pointSF, dmdist, ierr); call CHKERR(ierr)
        call DMDestroy(plex%dm, ierr); call CHKERR(ierr)
        plex%dm   = dmdist
        plex%comm = comm
      endif

      call update_plex_indices(plex)
    end subroutine
end module


program main
  use m_gen_plex_from_icon
  implicit none

  logical :: lflg_grid2d=.False., lflg_grid3d=.False., lflg_plex=.False., lflg

  character(len=default_str_len) :: gridfile

  type(t_plexgrid) :: plex
  type(t_icongrid),allocatable :: icongrid
  PetscInt :: petscint, Nz=3
  PetscScalar :: petscreal
  integer(mpiint) :: numnodes, myid, mpierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)

  call init_mpi_data_parameters(PETSC_COMM_WORLD)
  call mpi_comm_rank(PETSC_COMM_WORLD, myid, mpierr)
  call mpi_comm_size(PETSC_COMM_WORLD, numnodes, mpierr)

  print *,'Kind ints',kind(petscint), kind(i0)
  print *,'Kind reals',kind(petscreal), kind(one)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid2d', gridfile, lflg_grid2d, ierr); call CHKERR(ierr)
  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid3d', gridfile, lflg_grid3d, ierr); call CHKERR(ierr)
  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-plex', gridfile, lflg_plex, ierr); call CHKERR(ierr)
  if((.not.lflg_grid2d) .and. (.not.lflg_grid3d) .and. (.not.lflg_plex)) then
    print *,'Please specify a grid file with either option:'
    print *,'-grid2d <path_to_icon_grid_file.nc>'
    print *,'-grid3d <path_to_icon_grid_file.nc>'
    print *,'-plex <path_to_preprocessed_plex.h5>'
    stop 'Required Option missing'
  endif

  if(myid.eq.0) then
    print *,lflg_grid2d,' ', lflg_grid3d, ' ', lflg_plex
    if (lflg_grid2d) then
      call read_icon_grid_file(gridfile, icongrid)
      call create_dmplex_2d(icongrid, trim('plex'//gridfile), plex)
    else if (lflg_grid3d) then
      call read_icon_grid_file(gridfile, icongrid)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-Nz', Nz, lflg, ierr); call CHKERR(ierr)
      call create_dmplex_3d(PETSC_COMM_WORLD, icongrid, trim('plex3d'//gridfile), Nz, plex)
    else if (lflg_plex) then
      call load_plex_from_file(PETSC_COMM_WORLD, gridfile, plex)
    endif
  endif


  if(myid.eq.0) then
    call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
    call create_mass_vec(icongrid, plex)
  endif ! rank0

  call distribute_dmplex(PETSC_COMM_WORLD, plex)
  call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex_dist", ierr); call CHKERR(ierr)

  call DMDestroy(plex%dm, ierr); call CHKERR(ierr)
  call PetscFinalize(ierr)
  print *,'m_gen_plex_from_icon... done'
end program
