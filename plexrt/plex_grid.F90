module m_plex_grid
#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload
  use m_helper_functions, only: CHKERR, compute_normal_3d, spherical_2_cartesian, norm, cross_3d, &
    determine_normal_direction, angle_between_two_vec, rad2deg, deg2rad, hit_plane, &
    rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world, vec_proj_on_plane, &
    get_arg, imp_bcast, itoa
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, zero, one, pi, &
    default_str_len
  use m_icon_grid, only : t_icongrid

  implicit none

  private
  public :: t_plexgrid, load_plex_from_file, create_plex_from_icongrid, &
    icell_icon_2_plex, update_plex_indices, &
    compute_face_geometry, setup_edir_dmplex, print_dmplex,       &
    setup_abso_dmplex, compute_edir_absorption, create_edir_mat,  &
    TOP_BOT_FACE, SIDE_FACE

  logical, parameter :: ldebug=.True.

  integer(iintegers), parameter :: TOP_BOT_FACE=1, SIDE_FACE=2

  type :: t_plexgrid
    integer(mpiint) :: comm

    DM, allocatable :: dm

    DM, allocatable :: abso_dm

    DM, allocatable :: edir_dm

    DM, allocatable :: geom_dm
    type(tVec) :: geomVec

    character(len=8) :: boundary_label='boundary'

    real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

    ! Index counters on plex:
    integer(iintegers) :: pStart, pEnd ! points
    integer(iintegers) :: cStart, cEnd ! cells
    integer(iintegers) :: fStart, fEnd ! faces
    integer(iintegers) :: eStart, eEnd ! edges
    integer(iintegers) :: vStart, vEnd ! vertices

    integer(iintegers) :: Nz = 1      ! Number of layers
    integer(iintegers) :: Ncells      ! Number of cells in DMPlex (wedge form)
    integer(iintegers) :: Nfaces      ! Number of faces (for each cell 2 for top/bottom and 3 sideward-faces)
    integer(iintegers) :: Nedges      ! Number of edges (each cell has 3 at top/bottom faces plus 3 vertically)
    integer(iintegers) :: Nvertices   ! Number of vertices (each cell has 3 at top/bottom faces)

    ! In consecutive numbering in Plex chart we first have the cells until offset_faces-1
    integer(iintegers) :: offset_faces           ! then the top/bot faces
    integer(iintegers) :: offset_faces_sides     ! and next the side faces
    integer(iintegers) :: offset_edges
    integer(iintegers) :: offset_edges_vertical
    integer(iintegers) :: offset_vertices

    type(tDMLabel) :: faceposlabel    ! TOP_BOT_FACE=1, SIDE_FACE=2
    type(tDMLabel) :: iconindexlabel  ! local index of face, edge, vertex on icongrid
    type(tDMLabel) :: zindexlabel     ! vertical layer / level
    type(tDMLabel) :: TOAlabel        ! 1 if top level, 0 otherwise
    type(tDMLabel) :: ownerlabel      ! rank that posses this element
  end type

  contains

    subroutine create_plex_from_icongrid(comm, Nz, icongrid, plex)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nz
      type(t_icongrid), allocatable, intent(in) :: icongrid
      type(t_plexgrid), allocatable, intent(inout) :: plex

      integer(iintegers) :: chartsize
      !integer(iintegers) :: offset_faces, offset_edges, offset_vertices
      !integer(iintegers) :: offset_faces_sides, offset_edges_vertical
      integer(iintegers) :: edge3(3), edge4(4), faces(5), vert2(2)

      integer(iintegers) :: i, k, j, icell, iedge, iface, ivertex
      integer(iintegers) :: depth
      integer(mpiint) :: myid, numnodes, ierr

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(allocated(plex)) stop 'create_plex_from_icongrid :: plex should not be allocated!'
      if(.not.allocated(icongrid)) stop 'create_plex_from_icongrid :: icongrid should be allocated!'

      allocate(plex)
      plex%comm = comm

      plex%Nz = Nz

      plex%Ncells    = icongrid%Nfaces * plex%Nz
      plex%Nfaces    = icongrid%Nfaces * (plex%Nz+1) + icongrid%Nedges * plex%Nz
      plex%Nedges    = icongrid%Nedges * (plex%Nz+1) + icongrid%Nvertices * plex%Nz
      plex%Nvertices = icongrid%Nvertices * (plex%Nz+i1)

      chartsize = plex%Ncells + plex%Nfaces + plex%Nedges + plex%Nvertices

      allocate(plex%dm)
      call DMPlexCreate(plex%comm, plex%dm, ierr); call CHKERR(ierr)
      call DMSetDimension(plex%dm, i3, ierr); call CHKERR(ierr)

      call DMPlexSetChart(plex%dm, i0, chartsize, ierr); call CHKERR(ierr)

      plex%offset_faces = plex%Ncells
      plex%offset_faces_sides = plex%offset_faces + (plex%Nz+1)*icongrid%Nfaces
      plex%offset_edges = plex%Ncells + plex%Nfaces
      plex%offset_edges_vertical = plex%offset_edges + icongrid%Nedges * (plex%Nz+i1)
      plex%offset_vertices = plex%Ncells + plex%Nfaces + plex%Nedges
      print *,myid,'offsets faces:', plex%offset_faces, plex%offset_faces_sides
      print *,myid,'offsets edges:', plex%offset_edges, plex%offset_edges_vertical
      print *,myid,'offsets verte:', plex%offset_vertices
      print *,myid,'Chartsize:', chartsize

      ! Create some labels ... those are handy later when setting up matrices etc
      call DMCreateLabel(plex%dm, "Face Position", ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "Vertical Index", ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "Icon Index", ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "TOA", ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "Owner", ierr); call CHKERR(ierr)

      call DMGetLabel(plex%dm, "Face Position", plex%faceposlabel, ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "Vertical Index", plex%zindexlabel, ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "Icon Index", plex%iconindexlabel, ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "TOA", plex%TOAlabel, ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "Owner", plex%ownerlabel, ierr); call CHKERR(ierr)

      ! Preallocation
      ! Every cell has 5 faces
      do i = i0, plex%offset_faces-i1
        call DMPlexSetConeSize(plex%dm, i, i5, ierr); call CHKERR(ierr)
      enddo

      ! top/bottom faces have 3 edges
      do i = plex%offset_faces, plex%offset_faces_sides-i1
        call DMPlexSetConeSize(plex%dm, i, i3, ierr); call CHKERR(ierr)
      enddo

      ! side faces have 4 edges
      do i = plex%offset_faces_sides, plex%offset_edges-i1
        call DMPlexSetConeSize(plex%dm, i, i4, ierr); call CHKERR(ierr)
      enddo

      ! Edges have 2 vertices
      do i = plex%offset_edges, plex%offset_vertices-i1
        call DMPlexSetConeSize(plex%dm, i, i2, ierr); call CHKERR(ierr)
      enddo

      call DMSetUp(plex%dm, ierr); call CHKERR(ierr) ! Allocate space for cones

      ! Setup Connections
      ! First set five faces of cell
      do k = 1, plex%Nz
        do i = 1, icongrid%Nfaces
          edge3 = icongrid%edge_of_cell(i,:)

          faces(1) = iface_top_icon_2_plex(icongrid, plex, i, k) ! top face
          faces(2) = iface_top_icon_2_plex(icongrid, plex, i, k+1) ! bot face

          do j=1,3
            faces(2+j) = iface_side_icon_2_plex(icongrid, plex, edge3(j), k)
          enddo

          icell = icell_icon_2_plex(icongrid, plex, i, k)
          call DMPlexSetCone(plex%dm, icell, faces, ierr); call CHKERR(ierr)

          call DMLabelSetValue(plex%iconindexlabel, icell, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel, icell, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel, icell, icongrid%cellowner(icongrid%cell_index(i)), ierr); call CHKERR(ierr)
          !print *,myid,'Setting cell indices', i, icongrid%cell_index(i), 'edge3', edge3, 'faces', faces
        enddo
      enddo

      ! set edges of top/bot faces
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nfaces
          do j=1,3
            iedge = icongrid%edge_of_cell(i,j)
            edge3(j) = iedge_top_icon_2_plex(icongrid, plex, iedge, k)
          enddo
          iface = iface_top_icon_2_plex(icongrid, plex, i, k)

          call DMPlexSetCone(plex%dm, iface, edge3, ierr); call CHKERR(ierr)
          !print *,'edges @ horizontal faces',icell,':',edge3,'petsc:', offset_faces + Nfaces*(k-1) + i-i1, '::', edge3
          call DMLabelSetValue(plex%iconindexlabel, iface, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel, iface, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%faceposlabel, iface, TOP_BOT_FACE, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel, iface, icongrid%cellowner(icongrid%cell_index(i)), ierr); call CHKERR(ierr)
          if (k.eq.i1) then
            call DMLabelSetValue(plex%TOAlabel, iface, i1, ierr); call CHKERR(ierr)
          endif
        enddo
      enddo

      ! set edges of vertical faces
      do k = 1, plex%Nz ! layers
        do i = 1, icongrid%Nedges
          iface = iface_side_icon_2_plex(icongrid, plex, i, k)
          edge4(1) = iedge_top_icon_2_plex(icongrid, plex, i, k)
          edge4(2) = iedge_top_icon_2_plex(icongrid, plex, i, k+1)
          do j=1,2
            ivertex = icongrid%edge_vertices(i,j)
            edge4(2+j) = iedge_side_icon_2_plex(icongrid, plex, ivertex, k)
          enddo

          call DMPlexSetCone(plex%dm, iface, edge4, ierr); call CHKERR(ierr)
          !print *,myid,'edges @ vertical faces',i,icongrid%edge_index(i),':',edge4,'petsc:', iface
          call DMLabelSetValue(plex%iconindexlabel, iface, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel, iface, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%faceposlabel, iface, SIDE_FACE, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel, iface, icongrid%edgeowner(icongrid%edge_index(i)), ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      ! and then set the two vertices of edges in each level, i.e. vertices for edges in horizontal plane
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nedges
          iedge = iedge_top_icon_2_plex(icongrid, plex, i, k)
          do j=1,2
            ivertex = icongrid%edge_vertices(i,j)
            vert2(j) = ivertex_icon_2_plex(icongrid, plex, ivertex, k)
          enddo
          call DMPlexSetCone(plex%dm, iedge, vert2, ierr); call CHKERR(ierr)

          call DMLabelSetValue(plex%iconindexlabel, iedge, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel, iedge, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel, iedge, icongrid%edgeowner(icongrid%edge_index(i)), ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%faceposlabel, iedge, TOP_BOT_FACE, ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      ! and then set the two vertices of edges in each layer, i.e. vertices at the end of vertical edges
      do k = 1, plex%Nz ! layer
        do i = 1, icongrid%Nvertices
          iedge = iedge_side_icon_2_plex(icongrid, plex, i, k)
          vert2(1) = ivertex_icon_2_plex(icongrid, plex, i, k)
          vert2(2) = ivertex_icon_2_plex(icongrid, plex, i, k+1)
          call DMPlexSetCone(plex%dm, iedge, vert2, ierr); call CHKERR(ierr)

          call DMLabelSetValue(plex%iconindexlabel, iedge, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel, iedge, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel, iedge, icongrid%vertexowner(icongrid%vertex_index(i)), ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%faceposlabel, iedge, SIDE_FACE, ierr); call CHKERR(ierr)
        enddo
      enddo

      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nvertices
          ivertex = ivertex_icon_2_plex(icongrid, plex, i, k)

          call DMLabelSetValue(plex%iconindexlabel, ivertex, i, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%zindexlabel   , ivertex, k, ierr); call CHKERR(ierr)
          call DMLabelSetValue(plex%ownerlabel    , ivertex, icongrid%vertexowner(icongrid%vertex_index(i)), ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      print *,'Symmetrize'
      call DMPlexSymmetrize(plex%dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug) then
        call DMPlexGetDepth(plex%dm, depth, ierr); call CHKERR(ierr)
        print *,'Depth of Stratum:', depth
      endif

      call update_plex_indices(plex)
      !if(numnodes.gt.1) &
        call set_sf_graph(icongrid, plex)

      call set_coords()

      call mpi_barrier(plex%comm, ierr)
      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
      call mpi_barrier(plex%comm, ierr)

      contains

        subroutine set_coords()
          real(ireals), pointer :: coords(:)
          type(tVec)            :: coordinates
          integer(iintegers)    :: dimEmbed, coordSize, voff, ind
          type(tPetscSection)   :: coordSection

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
    subroutine set_sf_graph(icongrid, plexgrid)
      type(t_icongrid), intent(in) :: icongrid
      type(t_plexgrid), intent(inout) :: plexgrid

      type(tPetscSF) :: sf

      type(PetscInt) :: nroots, nleaves
      type(PetscInt),allocatable :: ilocal_elements(:)
      type(PetscSFNode),allocatable :: iremote_elements(:)
      type(PetscCopyMode),parameter :: localmode=PETSC_COPY_VALUES, remotemode=PETSC_COPY_VALUES

      integer(mpiint) :: myid, numnodes, ierr
      integer(iintegers) :: N_remote_cells, N_remote_edges, N_remote_vertices
      integer(iintegers) :: icell, iface, iedge, ivertex
      integer(iintegers) :: facepos, k, ilocal, iparent, iremote, owner, ileaf

      call mpi_comm_rank(plexgrid%comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(plexgrid%comm, numnodes, ierr); call CHKERR(ierr)

      N_remote_cells = 0
      do ilocal = 1, icongrid%Nfaces
        iparent = icongrid%cell_index(ilocal)
        if(icongrid%cellowner(iparent).ne.myid) N_remote_cells = N_remote_cells + 1
      enddo

      N_remote_edges = 0
      do ilocal = 1, icongrid%Nedges
        iparent = icongrid%edge_index(ilocal)
        if(icongrid%edgeowner(iparent).ne.myid) N_remote_edges = N_remote_edges + 1
      enddo

      N_remote_vertices = 0
      do ilocal = 1, icongrid%Nvertices
        iparent = icongrid%vertex_index(ilocal)
        if(icongrid%vertexowner(iparent).ne.myid) N_remote_vertices = N_remote_vertices + 1
      enddo

      nleaves = 0
      do k=plexgrid%pStart, plexgrid%pEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel, k, owner, ierr); call CHKERR(ierr)
        print *,myid,'ownership of plex element',k,owner
        if(owner.ne.myid) nleaves = nleaves + 1
      enddo

      !nleaves = N_remote_cells + N_remote_edges + N_remote_vertices ! non local elements in icongrid

      print *,myid,'remote elements:',nleaves, '(',N_remote_cells, N_remote_edges, N_remote_vertices, ')'

      allocate(ilocal_elements(nleaves))  ! local indices of elements
      allocate(iremote_elements(nleaves)) ! remote indices of elements

      ileaf = 0
      do icell = plexgrid%cStart, plexgrid%cEnd-1
        call DMLabelGetValue(plexgrid%zindexlabel, icell, k, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%ownerlabel, icell, owner, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%iconindexlabel, icell, ilocal, ierr); call CHKERR(ierr)

        if(owner.ne.myid) then
          iremote = icell_icon_2_plex(icongrid, plexgrid, ilocal, k, owner) ! plex index of cell at neighbor
          ileaf = ileaf + 1
          ilocal_elements(ileaf) = icell
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote

          print *,myid, 'local cell', icell, 'parent', iparent, '@', owner, '->', iremote
        endif
      enddo

      do iface = plexgrid%fStart, plexgrid%fEnd-1
        call DMLabelGetValue(plexgrid%zindexlabel, iface, k, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%iconindexlabel, iface, ilocal, ierr); call CHKERR(ierr) ! either a local index for cells or edges ...
        call DMLabelGetValue(plexgrid%faceposlabel, iface, facepos, ierr); call CHKERR(ierr) ! ... depending on faceposition

        if(facepos.eq.TOP_BOT_FACE) then
          iparent = icongrid%cell_index(ilocal) ! global cell index in parent grid
          iparent = icongrid%local_cell_index(iparent) ! local icon index @ neighbour
          iremote = iface_top_icon_2_plex(icongrid, plexgrid, iparent, k, owner) ! plex index at neighbor
          !print *,'top bot face: plex owner vs icongrid owner ::',iface, ilocal, iparent, '::', owner, icongrid%cellowner(iparent)
        elseif(facepos.eq.SIDE_FACE) then
          iparent = icongrid%edge_index(ilocal) ! global edge index in parent grid
          iparent = icongrid%local_edge_index(iparent) ! local icon index @ neighbour
          iremote = iface_side_icon_2_plex(icongrid, plexgrid, iparent, k, owner) ! plex index at neighbor

          !print *,'side face: plex owner vs icongrid owner ::',iface, ilocal, iparent, '::', owner, icongrid%edgeowner(iparent)
        else
          stop 'wrong or no side face label'
        endif

        if(owner.ne.myid) then
          ileaf = ileaf + 1
          ilocal_elements(ileaf) = iface
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          print *,myid, 'local cell', iface, 'parent', iparent, '@', owner, '->', iremote
        endif
      enddo

      do iedge = plexgrid%eStart, plexgrid%eEnd-1
        call DMLabelGetValue(plexgrid%zindexlabel, iedge, k, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%ownerlabel, iedge, owner, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%iconindexlabel, iedge, ilocal, ierr); call CHKERR(ierr) ! either a local index for edges or vertices ...
        call DMLabelGetValue(plexgrid%faceposlabel, iedge, facepos, ierr); call CHKERR(ierr) ! ... depending on faceposition

        if(facepos.eq.TOP_BOT_FACE) then
          iparent = icongrid%edge_index(ilocal) ! global cell index in parent grid
          iparent = icongrid%local_edge_index(iparent) ! local icon index @ neighbour
          iremote = iedge_top_icon_2_plex(icongrid, plexgrid, iparent, k, owner) ! plex index at neighbor
        elseif(facepos.eq.SIDE_FACE) then
          iparent = icongrid%vertex_index(ilocal) ! global edge index in parent grid
          iparent = icongrid%local_vertex_index(iparent) ! local icon index @ neighbour
          iremote = iedge_side_icon_2_plex(icongrid, plexgrid, iparent, k, owner) ! plex index at neighbor
        else
          stop 'wrong or no side face label'
        endif

        if(owner.ne.myid) then
          ileaf = ileaf + 1
          ilocal_elements(ileaf) = iedge
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          print *,myid, 'local edge', iedge, 'parent(edges/vertices)', iparent, '@', owner, '->', iremote
        endif
      enddo

      call mpi_barrier(plexgrid%comm, ierr)
      do ivertex = plexgrid%vStart, plexgrid%vEnd-1
        call DMLabelGetValue(plexgrid%zindexlabel   , ivertex, k, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%ownerlabel    , ivertex, owner, ierr); call CHKERR(ierr)
        call DMLabelGetValue(plexgrid%iconindexlabel, ivertex, ilocal, ierr); call CHKERR(ierr) ! local index for vertices ...

        iparent = icongrid%vertex_index(ilocal) ! global edge index in parent grid
        iparent = icongrid%local_vertex_index(iparent) ! local icon index @ neighbour
        iremote = ivertex_icon_2_plex(icongrid, plexgrid, iparent, k, owner) ! plex index at neighbor

        if(owner.ne.myid) then
          ileaf = ileaf + 1
          ilocal_elements(ileaf) = ivertex
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          print *,myid, 'local vertex', ivertex, 'parent', iparent, '@', owner, '->', iremote
        endif
      enddo

      if(ileaf.ne.nleaves) call abort('Seems like we forgot some remote elements?'//itoa(ileaf)//'/'//itoa(nleaves))

      nroots = plexgrid%pEnd

      call mpi_barrier(plexgrid%comm, ierr)
      print *,'calling sfsetgraph:',nroots,nleaves
      call DMGetPointSF(plexgrid%dm, sf, ierr); call CHKERR(ierr)
      call PetscSFSetGraph(sf, nroots, nleaves, ilocal, localmode, iremote, remotemode, ierr); call CHKERR(ierr)

      call mpi_barrier(plexgrid%comm, ierr)

      call PetscObjectViewFromOptions(sf, PETSC_NULL_SF, "-show_pointsf", ierr); call CHKERR(ierr)

      call PetscSFSetUp(sf, ierr); call CHKERR(ierr)

      call mpi_barrier(plexgrid%comm, ierr)
      stop 'debug'

      call DMSetDefaultSF(plexgrid%dm, sf, ierr); call CHKERR(ierr)
      call DMGetDefaultSF(plexgrid%dm, sf, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(sf, PETSC_NULL_SF, "-show_sf", ierr); call CHKERR(ierr)
      call PetscSFView(sf, PETSC_VIEWER_STDOUT_WORLD,ierr);call CHKERR(ierr);

    end subroutine

    subroutine load_plex_from_file(comm, gridfile, plex)
      integer(mpiint), intent(in) :: comm
      character(len=default_str_len),intent(in) :: gridfile
      type(t_plexgrid), intent(inout) :: plex
      integer(mpiint) :: numnodes, ierr

      type(tPetscSF)  :: sf
      DM       :: dmdist

      if(ldebug) print *,'Loading Plex Grid from File:', trim(gridfile)
      plex%comm = comm
      allocate(plex%dm)
      call DMPlexCreateFromFile(comm, trim(gridfile), PETSC_TRUE, plex%dm, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
      if (numnodes.gt.1_mpiint) then
        call DMPlexSetAdjacencyUseCone(plex%dm, PETSC_TRUE, ierr); call CHKERR(ierr)
        call DMPlexSetAdjacencyUseClosure(plex%dm, PETSC_TRUE, ierr); call CHKERR(ierr)
        call DMPlexDistribute(plex%dm, 0_mpiint, sf, dmdist, ierr); call CHKERR(ierr)
        call DMDestroy(plex%dm, ierr); call CHKERR(ierr)
        plex%dm   = dmdist
      endif
      call PetscObjectViewFromOptions(sf, PETSC_NULL_SF, "-show_plex_sf", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
    end subroutine

  subroutine compute_face_geometry(plex)
    type(t_plexgrid), intent(inout) :: plex

    type(tVec) :: coordinates
    integer(iintegers) :: cStart, cEnd, icell
    integer(iintegers) :: fStart, fEnd, iface
    integer(iintegers) :: vStart, vEnd, ivert
    integer(iintegers) :: Nedges, Nvertices

    integer(iintegers), pointer :: transclosure(:)
    real(ireals), pointer :: coords(:) ! pointer to coordinates vec

    integer(iintegers), target :: vertices3(3), vertices4(4), vertices6(6)
    integer(iintegers), pointer :: vertices(:)

    type(tPetscSection) :: coordSection, geomSection
    integer(iintegers) :: Ndim, voff

    real(ireals), allocatable :: vertex_coord(:,:) ! shape (Nvertices, dims)
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec

    integer(mpiint) :: ierr

    if(allocated(plex%geom_dm)) stop 'called compute_face_geometry on an already allocated geom_dm'
    allocate(plex%geom_dm)
    call DMClone(plex%dm, plex%geom_dm, ierr); ; call CHKERR(ierr)
    call DMGetCoordinateDim(plex%geom_dm, Ndim, ierr); call CHKERR(ierr)
    call DMGetCoordinateSection(plex%geom_dm, coordSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(coordSection, PETSC_NULL_SECTION, "-show_dm_coord_section", ierr); call CHKERR(ierr)

    call create_plex_section(plex%comm, plex%geom_dm, 'Geometry Section', i1*3, i1*6, i0, i0, geomSection)  ! Contains 3 dof for centroid on cells and faces plus 3 for normal vecs on faces
    call DMSetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call PetscSectionDestroy(geomSection, ierr); call CHKERR(ierr)
    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(geomSection, PETSC_NULL_SECTION, "-show_dm_geom_section", ierr); call CHKERR(ierr)


    call DMGetCoordinatesLocal(plex%geom_dm, coordinates, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_dm_coord", ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%geom_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr)  ! cells
    call DMPlexGetHeightStratum(plex%geom_dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (plex%geom_dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    call DMGetLocalVector(plex%geom_dm, plex%geomVec,ierr); call CHKERR(ierr)
    call PetscObjectSetName(plex%geomVec, 'geomVec', ierr);call CHKERR(ierr)
    call VecGetArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1
      call DMPlexGetTransitiveClosure(plex%geom_dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      !print *,'cell transclosure:',icell, transclosure
      if(size(transclosure).ne.21*2) then
        print *,'len of transclosure with size', size(transclosure), 'not supported -- is this a wedge?'
        stop('geometry not supported -- is this a wedge?')
      endif

      vertices => vertices6

      vertices = transclosure(size(transclosure)-11:size(transclosure):2)
      !print *,'cell vertices',vertices

      allocate(vertex_coord(Ndim, size(vertices)))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:Ndim+voff)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      ! print *,'centroid of cell:',icell,'::', sum(vertex_coord,dim=2)/size(vertices)
      call PetscSectionGetOffset(geomSection, icell, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2)/size(vertices)

      call DMPlexRestoreTransitiveClosure(plex%geom_dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      deallocate(vertex_coord)
    enddo

    do iface = fStart, fEnd-1
      call DMPlexGetConeSize(plex%geom_dm, iface, Nedges, ierr); call CHKERR(ierr)

      select case(Nedges)
      case (3)
        vertices => vertices3
        Nvertices = 3

      case (4)
        vertices => vertices4
        Nvertices = 4

      case default
        print *,'Computation of face normal not implemented for faces with ',Nedges,'edges/vertices'
        stop 'Invalid number of edges for face normal computation'
      end select

      call DMPlexGetTransitiveClosure(plex%geom_dm, iface, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      !print *,'transclosure', iface,'::',Nedges,'::',transclosure

      vertices = transclosure(3+2*Nedges:size(transclosure):2) ! indices come from: 1 for faceid, Nedge, and then vertices, each with two entries, one for index, and one for orientation
      !print *,'iface',iface,'vertices', vertices

      call DMPlexRestoreTransitiveClosure(plex%geom_dm, iface, PETSC_True, transclosure, ierr); call CHKERR(ierr)

      ! Get the coordinates of vertices
      allocate(vertex_coord(Ndim, Nvertices))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:voff+Ndim)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      !print *,'centroid of face:',iface,'::', sum(vertex_coord,dim=1)/Nvertices
      call PetscSectionGetOffset(geomSection, iface, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2)/Nvertices

      ! and use 3 coordinates to compute normal
      ! print *,'normal of face', iface,'::', compute_normal_3d(vertex_coord(1,:),vertex_coord(2,:),vertex_coord(3,:))
      geoms(voff+Ndim+i1: voff+Ndim*2) = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,2),vertex_coord(:,3))

      deallocate(vertex_coord)
    enddo

    call VecRestoreArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(plex%geomVec, PETSC_NULL_VEC, "-show_dm_geom_vec", ierr); call CHKERR(ierr)
  end subroutine

  subroutine create_plex_section(comm, dm, sectionname, cdof, fdof, edof, vdof, section)
    MPI_Comm, intent(in) :: comm
    type(tDM), intent(in) :: dm
    character(len=*) :: sectionname
    integer(iintegers), intent(in) :: cdof, fdof, edof, vdof
    type(tPetscSection), intent(out) :: section
    integer(iintegers)    :: i, section_size

    integer(iintegers) :: pStart, pEnd
    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(iintegers) :: eStart, eEnd
    integer(iintegers) :: vStart, vEnd

    integer(mpiint) :: ierr

    ! This is code to manually create a section, e.g. as in DMPlexCreateSection

    call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    ! Create Default Section
    call PetscSectionCreate(comm, section, ierr); call CHKERR(ierr)
    call PetscSectionSetNumFields(section, 1, ierr); call CHKERR(ierr)
    call PetscSectionSetChart(section, pStart, pEnd, ierr); call CHKERR(ierr)

    do i = cStart, cEnd-1
      call PetscSectionSetDof(section, i, cdof, ierr); call CHKERR(ierr)
      call PetscSectionSetFieldDof(section, i, 0, cdof, ierr); call CHKERR(ierr)
    enddo
    do i = fStart, fEnd-1
      call PetscSectionSetDof(section, i, fdof, ierr); call CHKERR(ierr)
      call PetscSectionSetFieldDof(section, i, 0, fdof, ierr); call CHKERR(ierr)
    enddo
    do i = eStart, eEnd-1
      call PetscSectionSetDof(section, i, edof, ierr); call CHKERR(ierr)
      call PetscSectionSetFieldDof(section, i, 0, edof, ierr); call CHKERR(ierr)
    enddo
    do i = vStart, vEnd-1
      call PetscSectionSetDof(section, i, vdof, ierr); call CHKERR(ierr)
      call PetscSectionSetFieldDof(section, i, 0, vdof, ierr); call CHKERR(ierr)
    enddo

    call PetscSectionSetUp(section, ierr); call CHKERR(ierr)

    call PetscSectionGetStorageSize(section, section_size, ierr);

    call PetscObjectSetName(section, trim(sectionname), ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, '-show_'//trim(sectionname), ierr); call CHKERR(ierr)
  end subroutine

  subroutine setup_edir_dmplex(plex, dm)
    type(t_plexgrid), intent(inout) :: plex
    type(tDM), allocatable :: dm
    type(tPetscSection) :: edirSection
    integer(mpiint) :: ierr

    if(allocated(dm)) stop 'called setup_edir_dmplex on an already allocated DM'
    allocate(dm)

    call DMClone(plex%dm, dm, ierr); call CHKERR(ierr)

    call PetscObjectSetName(dm, 'plex_edir', ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex_edir", ierr); call CHKERR(ierr)

    call create_plex_section(plex%comm, dm, 'Face Section', i0, i1*1, i0, i0, edirSection)  ! Contains 1 dof on each side for direct radiation
    call DMSetDefaultSection(dm, edirSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(edirSection, PETSC_NULL_SECTION, '-show_edir_section', ierr); call CHKERR(ierr)
    call PetscSectionDestroy(edirSection, ierr); call CHKERR(ierr)
  end subroutine

  subroutine setup_abso_dmplex(plex, dm)
    type(t_plexgrid), intent(inout) :: plex
    type(tDM), allocatable :: dm
    type(tPetscSection) :: s
    integer(mpiint) :: ierr


    if(allocated(dm)) stop 'called setup_abso_dmplex on an already allocated DM'
    allocate(dm)

    call DMClone(plex%dm, dm, ierr); call CHKERR(ierr)

    call PetscObjectSetName(dm, 'plex_abso', ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex_abso", ierr); call CHKERR(ierr)

    call create_plex_section(plex%comm, dm, 'Absorption Section', i1, i0, i0, i0, s)  ! Contains 1 dof on each cell
    call DMSetDefaultSection(dm, s, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(s, PETSC_NULL_SECTION, '-show_abso_section', ierr); call CHKERR(ierr)
    call PetscSectionDestroy(s, ierr); call CHKERR(ierr)
  end subroutine

  subroutine print_dmplex(comm, dm)
    integer(mpiint), intent(in) :: comm
    type(tDM),intent(in) :: dm

    integer(mpiint) :: myid, ierr
    integer(iintegers) :: pStart, pEnd
    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(iintegers) :: eStart, eEnd
    integer(iintegers) :: vStart, vEnd

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    print *,myid,'pStart,End :: ',pStart, pEnd
    print *,myid,'cStart,End :: ',cStart, cEnd
    print *,myid,'fStart,End :: ',fStart, fEnd
    print *,myid,'eStart,End :: ',eStart, eEnd
    print *,myid,'vStart,End :: ',vStart, vEnd

  end subroutine

  subroutine compute_edir_absorption(plex, edir, abso)
    type(t_plexgrid), intent(inout) :: plex
    type(tVec),intent(in) :: edir
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
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    real(ireals) :: angle_to_sun

    if(.not.allocated(plex%edir_dm) .or. .not.allocated(plex%abso_dm)) stop 'called compute_edir_absorption with a dm which is not allocated?'

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

      call DMPlexGetCone(plex%abso_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(edir_section, faces_of_cell(iface), edir_offset, ierr); call CHKERR(ierr)

        call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
        cell_center = geoms(1+geom_offset:3+geom_offset)

        call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); call CHKERR(ierr)
        face_center = geoms(1+geom_offset: 3+geom_offset)
        face_normal = geoms(4+geom_offset: 6+geom_offset)

        ! Determine the inward normal vec for the face
        face_normal = face_normal * determine_normal_direction(face_normal, face_center, cell_center)

        ! Then determine if the face is src(in line with the sun vec) or if it is destination(contra sun direction)
        angle_to_sun = angle_between_two_vec(face_normal, plex%sundir)

        if(angle_to_sun.lt.pi/2) then
          lsrc(iface) = .True.
          xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1)
        else
          lsrc(iface) = .False.
          xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1)
        endif
      enddo
      call DMPlexRestoreCone(plex%abso_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)
    enddo

    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
  end subroutine

  subroutine create_edir_mat(plex, A)
    type(t_plexgrid), intent(inout) :: plex
    type(tMat), intent(out) :: A

    type(tPetscSection) :: sec

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(mpiint) :: ierr

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: iface, irow, icol, icell, idst

    type(tDMLabel) :: faceposlabel, zindexlabel

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec

    real(ireals) :: zenith, azimuth

    integer(iintegers) :: base_face   ! index of face which is closest to sun angle, index regarding faces_of_cell
    integer(iintegers) :: left_face   ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: right_face  ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: upper_face  ! index of face which is top/bot of base face, index regarding faces_of_cell
    integer(iintegers) :: bottom_face ! index of face which is top/bot of base face, index regarding faces_of_cell

    real(ireals) :: dir2dir(5,5), S(5), T(5) ! Nface**2
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, sec, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(plex%edir_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(plex%edir_dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMGetLabel(plex%edir_dm, "Face Position", faceposlabel, ierr); call CHKERR(ierr)
    call DMGetLabel(plex%edir_dm, "Vertical Index", zindexlabel, ierr); call CHKERR(ierr)

    call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      call compute_local_wedge_ordering(icell, faces_of_cell, geomSection, geoms, faceposlabel, zindexlabel, plex%sundir, &
        zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)

      do iface = 1, size(faces_of_cell)

        if(lsrc(iface)) then ! this is really a source

          !retrieve coeffs:
          if (iface.eq.upper_face) then
            call compute_dir2dir_coeff(i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.bottom_face) then
            call compute_dir2dir_coeff(5*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.base_face) then
            call compute_dir2dir_coeff(2*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.left_face) then
            call compute_dir2dir_coeff(3*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.right_face) then
            call compute_dir2dir_coeff(4*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else
            stop 'iface is not in local wedge face numbering... something must have gone terribly wrong in compute_local_wedge_ordering'
          endif

          dir2dir([upper_face,base_face,left_face,right_face,bottom_face], iface) = T

          do idst=1,size(faces_of_cell)
            call PetscSectionGetOffset(sec, faces_of_cell(idst), irow, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
            call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces

            call MatSetValue(A, irow, icol, -dir2dir(idst, iface), INSERT_VALUES, ierr); call CHKERR(ierr)
          enddo

        endif

      enddo

      ! Set diagonal entries
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(sec, faces_of_cell(iface), irow, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
        call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
        call MatSetValue(A, irow, icol, one, INSERT_VALUES, ierr); call CHKERR(ierr)
      enddo

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_VEC, '-show_A', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

  end subroutine

  subroutine compute_local_wedge_ordering(icell, faces_of_cell, geomSection, geoms, faceposlabel, zindexlabel, sundir, &
      zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)
    integer(iintegers), intent(in) :: icell
    integer(iintegers), intent(in), pointer :: faces_of_cell(:)
    type(tPetscSection) :: geomSection
    real(ireals), intent(in), pointer :: geoms(:) ! pointer to coordinates vec
    type(tDMLabel), intent(in) :: faceposlabel, zindexlabel
    real(ireals), intent(in) :: sundir(3)
    real(ireals), intent(out) :: zenith, azimuth
    integer(iintegers), intent(out) :: upper_face, bottom_face, base_face, left_face, right_face
    logical, intent(out) :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: iface, geom_offset, ifacepos, izindex(2)

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
      face_normals(:,iface) = face_normals(:,iface) * determine_normal_direction(face_normals(:,iface), face_centers(:, iface), cell_center)
      ! to find the face whose normal is closest to sun vector:
      side_faces_angles_to_sun(iface) = angle_between_two_vec(face_normals(:,iface), sundir) ! in [rad]
      lsrc(iface) = side_faces_angles_to_sun(iface).lt.pi/2-deg2rad(.1_ireals) ! dont propagate energy along edges where sun is only .1 degrees off
    enddo

    ! get numbering for 2 top/bot faces and 3 side faces which indicate the position in the faces_of_cell vec
    iside_faces = 1; itop_faces = 1
    do iface = 1, size(faces_of_cell)
      call DMLabelGetValue(faceposlabel, faces_of_cell(iface), ifacepos, ierr); call CHKERR(ierr)
      select case(ifacepos)
      case(TOP_BOT_FACE)
        top_faces(itop_faces) = iface
        itop_faces = itop_faces+1
      case(SIDE_FACE)
        side_faces(iside_faces) = iface
        iside_faces = iside_faces+1
      case default
        stop 'wrong or no side face label'
      end select
    enddo

    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(1)), izindex(1), ierr); call CHKERR(ierr)
    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(2)), izindex(2), ierr); call CHKERR(ierr)
    if(izindex(1).gt.izindex(2)) then
      upper_face = top_faces(2)
      bottom_face = top_faces(1)
    else
      upper_face = top_faces(1)
      bottom_face = top_faces(2)
    endif

    do iface=1,size(side_faces)
      proj_normal = vec_proj_on_plane(sundir, face_normals(:,upper_face))
      proj_angles_to_sun(iface) = angle_between_two_vec(proj_normal, face_normals(:,side_faces(iface)))
    enddo

    ibase_face = minloc(proj_angles_to_sun,dim=1)
    base_face = side_faces(ibase_face)

    e_y = face_normals(:, base_face)   ! inward facing normal -> in local wedgemc coordinates
    e_z = -face_normals(:, upper_face) ! outward facing normal with respect to the top plate
    e_x = cross_3d(e_y, e_z)           ! in local wedge_coords, this is y=0 coordinate

    MrotWorld2Local = rotation_matrix_world_to_local_basis(e_x, e_y, e_z)

    left_face  = side_faces(modulo(ibase_face,size(side_faces))+1)
    right_face = side_faces(modulo(ibase_face+1,size(side_faces))+1)

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

    subroutine compute_dir2dir_coeff(src, phi, theta, S,T)
      use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_5
      integer(iintegers), intent(in) :: src
      real(ireals), intent(in) :: phi, theta
      real(ireals), intent(out) :: S(5),T(5)

      type(t_boxmc_wedge_5_5) :: bmc_wedge_5_5
      real(ireals) :: bg(3), dx,dy,dz
      real(ireals) :: S_tol(5),T_tol(5)

      call bmc_wedge_5_5%init(PETSC_COMM_SELF)
      !print *,'computing coeffs for src/phi/theta',src,phi,theta

      bg  = [1e-3_ireals, zero, one/2 ]

      !phi   =  0
      !theta = 45

      dx = 100
      dy = dx
      dz = 200

      call bmc_wedge_5_5%get_coeff(PETSC_COMM_SELF, bg, src, .True., &
        phi, theta, dx, dy, dz, S, T, S_tol, T_tol, inp_atol=1e-2_ireals, inp_rtol=1e-1_ireals)
    end subroutine

    !> @brief return the dmplex cell index for an icon base grid cell index
    function icell_icon_2_plex(icon, plex, icell, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: Nfaces
      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
      else
        Nfaces = icon%Nfaces
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz) then
          print *,'icell_icon_2_plex :: inp',icell, k, present(owner),'::',Nfaces
          stop 'icell_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1 .or. icell.gt.Nfaces) then
          print *,'icell_icon_2_plex :: inp',icell, k, present(owner),'::',Nfaces
          stop 'icell_icon_2_plex :: icon cell index out of range'
        endif
      endif
      icell_icon_2_plex = (k-i1)*Nfaces + icell - i1
    end function

    !> @brief return the dmplex face index for an icongrid index situated at the top of a cell
    function iface_top_icon_2_plex(icon, plex, icell, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: iface_top_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: Nfaces, offset
      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
        offset = Nfaces * plex%Nz
      else
        Nfaces = icon%Nfaces
        offset = plex%offset_faces
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz+1) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner), '::', Nfaces
          stop 'iface_top_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1 .or. icell.gt.Nfaces) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner), '::', Nfaces
          stop 'iface_top_icon_2_plex :: icon cell index out of range'
        endif
      endif
      iface_top_icon_2_plex = offset + Nfaces*(k-1) + icell-i1
    end function

    !> @brief return the dmplex face index for an icongrid index situated at the top of a cell
    function iface_side_icon_2_plex(icon, plex, iedge, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: iface_side_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: Nfaces, Nedges, offset

      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
        Nedges = icon%parNedges(owner)
        offset = Nfaces * plex%Nz + (plex%Nz + 1) * Nfaces
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        offset = plex%offset_faces_sides
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz) then
          print *,'iface_side_icon_2_plex :: inp', iedge, k, present(owner), '::', Nfaces, Nedges, offset
          call abort('iface_side_icon_2_plex :: vertical index k out of range '//itoa(k))
        endif
        if(iedge.lt.i1 .or. iedge.gt.Nedges) then
          print *,'iface_side_icon_2_plex :: inp', iedge, k, present(owner), '::', Nfaces, Nedges, offset
          call abort('iface_side_icon_2_plex :: icon edge index out of range '//itoa(iedge))
        endif
      endif
      iface_side_icon_2_plex = offset + (k-1)*Nedges + (iedge-1)
    end function

    !> @brief return the dmplex edge index for a given icon edge index, i.e. the edges on the top/bot faces of cells
    function iedge_top_icon_2_plex(icon, plex, iedge, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: iedge_top_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: Nfaces, Nedges, offset
      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
        Nedges = icon%parNedges(owner)
        offset = Nfaces * plex%Nz + Nfaces * (plex%Nz+1) + Nedges * plex%Nz
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        offset = plex%offset_edges
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz+1) stop 'iedge_top_icon_2_plex :: vertical index k out of range'
        if(iedge.lt.i1 .or. iedge.gt.Nedges) stop 'iedge_top_icon_2_plex :: icon cell index out of range'
      endif

      iedge_top_icon_2_plex = offset + Nedges*(k-1) + iedge-1
    end function


    !> @brief return the dmplex edge index for a given icon vertex index, i.e. the edges on the side faces of cells
    function iedge_side_icon_2_plex(icon, plex, ivertex, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: iedge_side_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: Nfaces, Nedges, Nvertices, offset
      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
        Nedges = icon%parNedges(owner)
        Nvertices = icon%parNvertices(owner)
        offset = Nfaces * plex%Nz + Nfaces * (plex%Nz+1) + Nedges * plex%Nz + Nedges * (plex%Nz+i1)
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        Nvertices = icon%Nvertices
        offset = plex%offset_edges_vertical
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz) stop 'iedge_side_icon_2_plex :: vertical index k out of range'
        if(ivertex.lt.i1 .or. ivertex.gt.Nvertices) stop 'iedge_side_icon_2_plex :: icon vertex index out of range'
      endif

      iedge_side_icon_2_plex = offset + Nvertices*(k-1) + ivertex-1
    end function

    !> @brief return the dmplex vertex index for a given icon vertex index
    function ivertex_icon_2_plex(icon, plex, ivertex, k, owner)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: ivertex_icon_2_plex !< @param[out] icell_icon_2_plex, the vertex index in the dmplex, starts from plex%vStart and goes to plex%vEnd
      integer(iintegers) :: Nfaces, Nedges, Nvertices, offset
      if(present(owner)) then
        Nfaces = icon%parNfaces(owner)
        Nedges = icon%parNedges(owner)
        Nvertices = icon%parNvertices(owner)
        offset = Nfaces * plex%Nz + Nfaces * (plex%Nz+1) + Nedges * plex%Nz + Nedges * (plex%Nz+1) + Nvertices * plex%Nz
      else
        Nfaces    = icon%Nfaces
        Nedges    = icon%Nedges
        Nvertices = icon%Nvertices
        offset    = plex%offset_vertices
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz+1) then
          print *,'ivertex_icon_2_plex error! input was',ivertex, k
          stop 'ivertex_side_icon_2_plex :: vertical index k out of range'
        endif
        if(ivertex.lt.i1 .or. ivertex.gt.Nvertices) stop 'ivertex_side_icon_2_plex :: icon vertex index out of range'
      endif

      ivertex_icon_2_plex = offset + Nvertices*(k-1) + ivertex-1
    end function

    !!> @brief return the vertical and icon base grid cell index for a given dmplex cell index.
    !function icell_plex_2_icon(icon, plex, ielement)
    !  type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
    !  type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
    !  integer(iintegers),intent(in) :: icell    !< @param[in] ielement starts with 0 up to cEnd-1 (size of DMPlex cells)
    !  integer(iintegers) :: icell_plex_2_icon(2)     !< @param[out] return icon index for base grid(starting with 1, to Nfaces) and vertical index(1 at top of domain)
    !  if(ldebug) then
    !    if(icell.lt.0 .or. icell.ge.plex%Ncells) then
    !      print *,'icell_plex_2_icon :: ', icell, ' not in valid range', plex%Ncells
    !      stop 'icell_plex_2_icon, dmplex cell index out of range'
    !    endif
    !  endif
    !  icell_plex_2_icon(1) = modulo(icell, icon%Nfaces) + i1
    !  icell_plex_2_icon(2) = icell / icon%Nfaces + i1
    !  if(ldebug) then
    !    if(icell_plex_2_icon(1).lt.i1 .or. icell_plex_2_icon(1).gt.icon%Nfaces) then
    !      print *,'icell_plex_2_icon error! input was',icell, icon%Nfaces, 'result:', icell_plex_2_icon
    !      stop 'icell_icon_2_plex :: icon cell index out of range'
    !    endif
    !    if(icell_plex_2_icon(2).lt.i1 .or. icell_plex_2_icon(2).gt.plex%Nz) then
    !      print *,'icell_plex_2_icon error! input was',icell, icon%Nfaces, 'result:', icell_plex_2_icon
    !      stop 'icell_icon_2_plex :: vertical index k out of range'
    !    endif
    !  endif
    !end function

    subroutine update_plex_indices(plex)
      type(t_plexgrid), intent(inout) :: plex
      integer(mpiint) :: myid, ierr

      call DMPlexGetChart(plex%dm, plex%pStart, plex%pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(plex%dm, i0, plex%cStart, plex%cEnd, ierr); call CHKERR(ierr) ! cells
      call DMPlexGetHeightStratum(plex%dm, i1, plex%fStart, plex%fEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetDepthStratum (plex%dm, i1, plex%eStart, plex%eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetDepthStratum (plex%dm, i0, plex%vStart, plex%vEnd, ierr); call CHKERR(ierr) ! vertices

      if(ldebug) then
        call mpi_comm_rank( plex%comm, myid, ierr)
        print *,myid, 'pstart', plex%pstart, 'pEnd', plex%pEnd
        print *,myid, 'cStart', plex%cStart, 'cEnd', plex%cEnd
        print *,myid, 'fStart', plex%fStart, 'fEnd', plex%fEnd
        print *,myid, 'eStart', plex%eStart, 'eEnd', plex%eEnd
        print *,myid, 'vStart', plex%vStart, 'vEnd', plex%vEnd
      endif
    end subroutine
end module
