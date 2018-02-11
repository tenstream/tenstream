module m_plex_grid
#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload
  use m_helper_functions, only: CHKERR, itoa, compute_normal_3d, approx, strF2C, distance, &
    triangle_area_by_vertices, swap, norm
  use m_data_parameters, only : ireals, iintegers, mpiint, zero, &
    i0, i1, i2, i3, i4, i5, i6, i7, default_str_len
  use m_icon_grid, only : t_icongrid, ICONULL

  implicit none

  private
  public :: t_plexgrid, load_plex_from_file, create_plex_from_icongrid, &
    icell_icon_2_plex, iface_top_icon_2_plex, update_plex_indices, distribute_plexgrid_dm, &
    compute_face_geometry, setup_edir_dmplex, print_dmplex,       &
    setup_abso_dmplex, ncvar2d_to_globalvec, facevec2cellvec

  logical, parameter :: ldebug=.True.

  type :: t_plexgrid
    integer(mpiint) :: comm

    type(tDM), allocatable :: dm
    type(tDM), allocatable :: abso_dm
    type(tDM), allocatable :: edir_dm
    type(tDM), allocatable :: geom_dm
    type(tVec) :: geomVec
    type(tPetscSF) :: default_sf

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

    logical,allocatable :: ltopfacepos(:)                ! TOP_BOT_FACE or SIDE_FACE of faces and edges, fStart..eEnd-1
    integer(iintegers),allocatable :: zindex(:)          ! vertical layer / level of cells/faces/edges/vertices , pStart..pEnd-1
    integer(iintegers),allocatable :: localiconindex(:)  ! local index of face, edge, vertex on icongrid, pStart, pEnd-1

    !type(tDMLabel) :: faceposlabel         ! TOP_BOT_FACE=1, SIDE_FACE=2
    !type(tDMLabel) :: iconindexlabel       ! index of face, edge, vertex on icongrid
    !type(tDMLabel) :: localiconindexlabel  ! local index of face, edge, vertex on icongrid
    !type(tDMLabel) :: zindexlabel          ! vertical layer / level
    type(tDMLabel) :: TOAlabel             ! 1 if top level, 0 otherwise
    type(tDMLabel) :: boundarylabel        ! 1 if top level, 2 if side face, 0 otherwise
    type(tDMLabel) :: domainboundarylabel  ! 1 if top level, 2 if side face, 0 otherwise
    type(tDMLabel) :: ownerlabel           ! rank that posses this element

    AO :: cell_ao

    real(ireals),allocatable :: hhl(:) ! vertical height of horizontal faces, i.e. height of levels, starting at TOA
  end type

  contains

    subroutine create_plex_from_icongrid(comm, Nz, hhl, cell_ao, icongrid, plex)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nz
      real(ireals), intent(in) :: hhl(:)
      AO, intent(in) :: cell_ao
      type(t_icongrid), allocatable, intent(in) :: icongrid
      type(t_plexgrid), allocatable, intent(inout) :: plex

      integer(iintegers) :: chartsize
      integer(iintegers) :: edge3(3), edge4(4), faces(5), vert2(2)

      integer(iintegers) :: i, j, k, icell, iedge, iface, ivertex
      integer(iintegers) :: iparent, owner
      integer(mpiint)    :: myid, numnodes, ierr

      type(tPetscSection) :: cell_section

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: starting'

      if(allocated(plex)) stop 'create_plex_from_icongrid :: plex should not be allocated!'
      if(.not.allocated(icongrid)) stop 'create_plex_from_icongrid :: icongrid should be allocated!'

      allocate(plex)
      plex%comm = comm
      plex%cell_ao = cell_ao

      plex%Nz = Nz
      if(size(hhl).ne.Nz+1) stop 'plex_grid::create_plex_from_icongrid -> hhl does not fit the size of Nz+1'
      allocate(plex%hhl(Nz+1), source=hhl)

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

      plex%pStart = 0                    ; plex%pEnd = plex%offset_vertices + plex%Nvertices
      plex%cStart = 0                    ; plex%cEnd = plex%Ncells
      plex%fStart = plex%Ncells          ; plex%fEnd = plex%offset_edges
      plex%eStart = plex%offset_edges    ; plex%eEnd = plex%offset_vertices
      plex%vStart = plex%offset_vertices ; plex%vEnd = plex%offset_vertices + plex%Nvertices

      print *,myid, 'pstart', plex%pstart, 'pEnd', plex%pEnd
      print *,myid, 'cStart', plex%cStart, 'cEnd', plex%cEnd
      print *,myid, 'fStart', plex%fStart, 'fEnd', plex%fEnd
      print *,myid, 'eStart', plex%eStart, 'eEnd', plex%eEnd
      print *,myid, 'vStart', plex%vStart, 'vEnd', plex%vEnd


      ! Create some labels ... those are handy later when setting up matrices etc
      call DMCreateLabel(plex%dm, "TOA"     , ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "Boundary", ierr); call CHKERR(ierr)
      call DMCreateLabel(plex%dm, "Owner"   , ierr); call CHKERR(ierr)

      call DMGetLabel(plex%dm, "TOA"     , plex%TOAlabel           , ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "Boundary", plex%boundarylabel, ierr); call CHKERR(ierr)
      call DMGetLabel(plex%dm, "Owner"   , plex%ownerlabel         , ierr); call CHKERR(ierr)
      call DMLabelSetDefaultValue(plex%ownerlabel, int(myid, kind=iintegers), ierr); call CHKERR(ierr)

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

      allocate(plex%ltopfacepos   (plex%fStart:plex%eEnd-1))
      allocate(plex%zindex        (plex%pStart:plex%pEnd-1))
      allocate(plex%localiconindex(plex%pStart:plex%pEnd-1))

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections'
      ! Setup Connections
      ! First set five faces of cell
      do k = 1, plex%Nz
        do i = 1, icongrid%Nfaces
          edge3 = icongrid%edge_of_cell(i,:)

          faces(1) = iface_top_icon_2_plex(plex, i, k) ! top face
          faces(2) = iface_top_icon_2_plex(plex, i, k+1) ! bot face

          do j=1,3
            faces(2+j) = iface_side_icon_2_plex(icongrid, plex, edge3(j), k)
          enddo

          icell = icell_icon_2_plex(plex, i, k)
          iparent = icongrid%cell_index(i)
          owner = icongrid%cellowner(iparent)

          call DMPlexSetCone(plex%dm, icell, faces, ierr); call CHKERR(ierr)
          call DMPlexSetConeOrientation(plex%dm, icell, i0*faces, ierr); call CHKERR(ierr)

          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel         , icell, owner  , ierr); call CHKERR(ierr)
          plex%zindex(icell) = k
          plex%localiconindex(icell) = i
          !print *,myid,'cells :: icon', i, icongrid%cell_index(i),'plex:',icell, 'edge3', edge3, 'faces', faces
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : edges of top/bot faces'
      ! set edges of top/bot faces
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nfaces
          do j=1,3
            iedge = icongrid%edge_of_cell(i,j)
            edge3(j) = iedge_top_icon_2_plex(icongrid, plex, iedge, k)
          enddo
          iface = iface_top_icon_2_plex(plex, i, k)
          iparent = icongrid%cell_index(i)
          owner = icongrid%cellowner(iparent)

          call DMPlexSetCone(plex%dm, iface, edge3, ierr); call CHKERR(ierr)
          call DMPlexSetConeOrientation(plex%dm, iface, i0*edge3, ierr); call CHKERR(ierr)

          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel         , iface, owner       , ierr); call CHKERR(ierr)
          plex%ltopfacepos(iface) = .True.
          plex%zindex(iface) = k
          plex%localiconindex(iface) = i

          if (k.eq.i1) then
            call DMLabelSetValue(plex%TOAlabel, iface, i1, ierr); call CHKERR(ierr)
            call DMLabelSetValue(plex%boundarylabel, iface, i1, ierr); call CHKERR(ierr)
          endif
          !print *,myid,'top/bot faces :: icon face', i, iparent,'plexface:',iface, 'edge3', edge3, 'owner', owner
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : edges of vertical faces'
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
          iparent = icongrid%edge_index(i)
          owner = icongrid%edgeowner(iparent)

          call DMPlexSetCone(plex%dm, iface, edge4, ierr); call CHKERR(ierr)
          call DMPlexSetConeOrientation(plex%dm, iface, i0*edge4, ierr); call CHKERR(ierr)

          plex%ltopfacepos(iface) = .False.
          plex%zindex(iface) = k
          plex%localiconindex(iface) = i
          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
          if(any(icongrid%adj_cell_of_edge(i,:).eq.ICONULL)) then
            call DMLabelSetValue(plex%boundarylabel, iface, i2, ierr); call CHKERR(ierr)
          endif

          !print *,myid,'side faces :: icon edge', i, iparent, 'plexface:', iface, 'edge4', edge4, 'owner', owner, &
          !  'adj cells', icongrid%adj_cell_of_edge(i,:)
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : vertices of horizontal edges'
      ! and then set the two vertices of edges in each level, i.e. vertices for edges in horizontal plane
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nedges
          iedge = iedge_top_icon_2_plex(icongrid, plex, i, k)
          do j=1,2
            ivertex = icongrid%edge_vertices(i,j)
            vert2(j) = ivertex_icon_2_plex(icongrid, plex, ivertex, k)
          enddo

          iparent = icongrid%edge_index(i)
          owner = icongrid%edgeowner(iparent)

          call DMPlexSetCone(plex%dm, iedge, vert2, ierr); call CHKERR(ierr)
          call DMPlexSetConeOrientation(plex%dm, iedge, i0*vert2, ierr); call CHKERR(ierr)

          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel         , iedge, owner       , ierr); call CHKERR(ierr)
          plex%zindex(iedge) = k
          plex%localiconindex(iedge) = i
          plex%ltopfacepos(iedge) = .True.
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : vertices of vertical edges'
      ! and then set the two vertices of edges in each layer, i.e. vertices at the end of vertical edges
      do k = 1, plex%Nz ! layer
        do i = 1, icongrid%Nvertices
          iedge = iedge_side_icon_2_plex(icongrid, plex, i, k)
          vert2(1) = ivertex_icon_2_plex(icongrid, plex, i, k)
          vert2(2) = ivertex_icon_2_plex(icongrid, plex, i, k+1)

          iparent = icongrid%vertex_index(i)
          owner = icongrid%vertexowner(iparent)

          call DMPlexSetCone(plex%dm, iedge, vert2, ierr); call CHKERR(ierr)
          call DMPlexSetConeOrientation(plex%dm, iedge, i0*vert2, ierr); call CHKERR(ierr)

          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel         , iedge, owner    , ierr); call CHKERR(ierr)
          plex%zindex(iedge) = k
          plex%localiconindex(iedge) = i
          plex%ltopfacepos(iedge) = .False.
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : labels of vertices'
      do k = 1, plex%Nz+1 ! levels
        do i = 1, icongrid%Nvertices
          ivertex = ivertex_icon_2_plex(icongrid, plex, i, k)

          iparent = icongrid%vertex_index(i)
          owner = icongrid%vertexowner(iparent)

          if(owner.ne.myid) call DMLabelSetValue(plex%ownerlabel         , ivertex, owner  , ierr); call CHKERR(ierr)
          plex%zindex(ivertex) = k
          plex%localiconindex(ivertex) = i
        enddo
      enddo

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : symmetrize dm'
      call DMPlexSymmetrize(plex%dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : putting up sf_graph'
      call set_sf_graph(icongrid, plex)

      call update_plex_indices(plex)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : set coords'
      call set_coords(plex, icongrid)

      call PetscObjectSetName(plex%dm, 'Icon DMPLEX', ierr);call CHKERR(ierr)

      call DMSetFromOptions(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : create default section'
      call create_plex_section(plex%comm, plex%dm, 'cell section', i1, [i1], [i0], [i0], [i0], cell_section)  ! Contains 1 dof for centroid on cells
      call DMSetDefaultSection(plex%dm, cell_section, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(cell_section, ierr); call CHKERR(ierr)

      !call DMCreateLabel(plex%dm, "Boundary", ierr); call CHKERR(ierr)
      !call DMGetLabel(plex%dm, "Boundary", plex%boundarylabel, ierr); call CHKERR(ierr)
      !call DMPlexMarkBoundaryFaces(plex%dm, i1, plex%boundarylabel, ierr); call CHKERR(ierr)
      !call DMPlexLabelAddCells(plex%dm, plex%boundarylabel, ierr); call CHKERR(ierr)
      call label_domain_boundary(plex, plex%dm, plex%boundarylabel, plex%domainboundarylabel)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : show plex'
      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex_dump", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: finished'

      call dump_ownership(icongrid, plex)
      end subroutine
      subroutine label_domain_boundary(plex, dm, boundarylabel, domainboundarylabel)
        type(t_plexgrid), intent(in) :: plex
        type(tDM), intent(in) :: dm
        type(tDMLabel), intent(out) :: boundarylabel, domainboundarylabel

        type(tDM) :: facedm
        type(tPetscSection) :: facesection
        type(tVec) :: lVec, gVec
        real(ireals), pointer :: xv(:)
        type(tIS) :: boundary_ids
        integer(iintegers), pointer :: xbndry_iface(:)
        integer(iintegers) :: i, iface, voff, lv
        integer(mpiint) :: myid, ierr

        call DMCreateLabel(dm, "DomainBoundary", ierr); call CHKERR(ierr)
        call DMGetLabel(dm, "DomainBoundary", domainboundarylabel, ierr); call CHKERR(ierr)

        call DMCreateLabel(dm, "Boundary", ierr); call CHKERR(ierr)
        call DMGetLabel(dm, "Boundary", boundarylabel, ierr); call CHKERR(ierr)
        call DMPlexMarkBoundaryFaces(dm, i1, boundarylabel, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_Boundary_DM', ierr); call CHKERR(ierr)

        call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

        call DMClone(dm, facedm, ierr); call CHKERR(ierr)

        call create_plex_section(plex%comm, facedm, 'Face Section', i1, [i0], [i1], [i0], [i0], facesection)  ! Contains 1 dof on each side
        call DMSetDefaultSection(facedm, facesection, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(facesection, ierr); call CHKERR(ierr)
        call DMGetDefaultSection(facedm, facesection, ierr); call CHKERR(ierr)

        call DMGetGlobalVector(facedm, gVec, ierr); call CHKERR(ierr)
        call VecSet(gVec, zero, ierr); call CHKERR(ierr)

        call DMGetLocalVector(facedm, lVec, ierr); call CHKERR(ierr)
        call VecSet(lVec, zero, ierr); call CHKERR(ierr)

        call DMGetStratumIS(dm, "Boundary", i1, boundary_ids, ierr); call CHKERR(ierr)

        if (boundary_ids.ne.PETSC_NULL_IS) then
          call ISGetIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)

          call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
          do i = 1, size(xbndry_iface)
            iface = xbndry_iface(i)
            call DMLabelGetValue(boundarylabel, iface, lv, ierr); call CHKERR(ierr)
            call PetscSectionGetOffset(facesection, iface, voff, ierr); call CHKERR(ierr)
            xv(i1+voff) = lv
          enddo
          call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)
          call ISRestoreIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)
        endif

        call DMLocalToGlobalBegin(facedm, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd  (facedm, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(facedm, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd  (facedm, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)

        call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        do iface = plex%fStart, plex%fEnd-1
          call PetscSectionGetOffset(facesection, iface, voff, ierr); call CHKERR(ierr)
          if(int(xv(i1+voff), iintegers) .eq. i1) then ! if the additive val is not 2 it must be at the domain edge
            if(plex%ltopfacepos(iface)) then
              if(plex%zindex(iface).eq.1) then
                call DMLabelSetValue(domainboundarylabel, iface, i1, ierr); call CHKERR(ierr)
              else
                call DMLabelSetValue(domainboundarylabel, iface, i3, ierr); call CHKERR(ierr)
              endif
            else
              call DMLabelSetValue(domainboundarylabel, iface, i2, ierr); call CHKERR(ierr)
            endif
          endif
        enddo
        call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)

        call DMRestoreLocalVector(facedm, lVec, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(facedm, gVec, ierr); call CHKERR(ierr)
      end subroutine

      subroutine dump_ownership(icongrid, plex)
        type(t_icongrid), intent(in) :: icongrid
        type(t_plexgrid), intent(inout) :: plex
        type(tVec) :: globalVec
        real(ireals), pointer :: xv(:)
        type(tPetscSection) :: cellSection
        integer(iintegers) :: icell, cStart, cEnd, voff, labelval

        integer(mpiint) :: myid, ierr

        call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

        if(ldebug.and.myid.eq.0) print *,'dump_ownership :: starting'

        call DMGetDefaultSection(plex%dm, cellSection, ierr); call CHKERR(ierr)

        call DMPlexGetHeightStratum(plex%dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

        ! Now lets get vectors!
        call DMGetGlobalVector(plex%dm, globalVec,ierr); call CHKERR(ierr)

        call PetscObjectSetName(globalVec, 'ownershipvec', ierr);call CHKERR(ierr)
        call VecGetArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        do icell = cStart, cEnd-1
          call DMLabelGetValue(plex%ownerlabel, icell, labelval, ierr); call CHKERR(ierr)
          call PetscSectionGetOffset(cellSection, icell, voff, ierr); call CHKERR(ierr)
          xv(voff+1) = real(labelval, kind=ireals)
        enddo
        call VecRestoreArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(globalVec, PETSC_NULL_VEC, '-show_ownership', ierr); call CHKERR(ierr)

        call PetscObjectSetName(globalVec, 'iconindexvec', ierr);call CHKERR(ierr)
        call VecGetArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        do icell = cStart, cEnd-1
          labelval = icongrid%cell_index(plex%localiconindex(icell))
          call PetscSectionGetOffset(cellSection, icell, voff, ierr); call CHKERR(ierr)
          xv(voff+1) = real(labelval, kind=ireals)
        enddo
        call VecRestoreArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(globalVec, PETSC_NULL_VEC, '-show_iconindex', ierr); call CHKERR(ierr)

        call PetscObjectSetName(globalVec, 'zindexlabel', ierr);call CHKERR(ierr)
        call VecGetArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        do icell = cStart, cEnd-1
          labelval = plex%zindex(icell)
          call PetscSectionGetOffset(cellSection, icell, voff, ierr); call CHKERR(ierr)
          xv(voff+1) = real(labelval, kind=ireals)
        enddo
        call VecRestoreArrayF90(globalVec, xv, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(globalVec, PETSC_NULL_VEC, '-show_zindex', ierr); call CHKERR(ierr)

        call DMRestoreGlobalVector(plex%dm, globalVec, ierr); call CHKERR(ierr)
        if(ldebug.and.myid.eq.0) print *,'dump_ownership :: finished'
      end subroutine

      subroutine facevec2cellvec(plex, faceVec_dm, global_faceVec)
        type(t_plexgrid), intent(in) :: plex
        type(tDM), intent(in) :: faceVec_dm
        type(tVec),intent(in) :: global_faceVec

        type(tVec) :: cellVec, global_cellVec, faceVec
        real(ireals), pointer :: xcellVec(:)
        real(ireals), pointer :: xfaceVec(:)

        type(tDM) :: celldm
        type(tPetscSection) :: cellSection, faceVecSection
        integer(iintegers) :: icell, cStart, cEnd, cell_offset, iface, faceVec_offset
        integer(iintegers),pointer :: faces_of_cell(:)

        integer(mpiint) :: myid, ierr
        character(len=default_str_len) :: faceVecname, cellVecname

        call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

        call PetscObjectGetName(global_faceVec, faceVecname, ierr); call CHKERR(ierr)
        if(ldebug.and.myid.eq.0) print *,'facevec2cellvec :: starting..'//trim(faceVecname)

        call DMClone(plex%dm, celldm, ierr); ; call CHKERR(ierr)

        call create_plex_section(plex%comm, celldm, 'Faces to Cells Section', i1, [i5], [i0], [i0], [i0], cellSection)
        call DMSetDefaultSection(celldm, cellSection, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(cellSection, ierr); call CHKERR(ierr)

        call DMGetDefaultSection(celldm, cellSection, ierr); call CHKERR(ierr)
        call DMGetDefaultSection(faceVec_dm, faceVecSection, ierr); call CHKERR(ierr)

        call DMGetLocalVector(faceVec_dm, faceVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(faceVec_dm, global_faceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd  (faceVec_dm, global_faceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

        call DMGetLocalVector(celldm, cellVec, ierr); call CHKERR(ierr)
        call VecSet(cellVec, zero, ierr); call CHKERR(ierr)

        call VecGetArrayF90(cellVec, xcellVec, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(faceVec, xfaceVec, ierr); call CHKERR(ierr)

        call DMPlexGetHeightStratum(celldm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
        do icell = cStart, cEnd-1
          call PetscSectionGetOffset(cellSection, icell, cell_offset, ierr); call CHKERR(ierr)

          call DMPlexGetCone(celldm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
          do iface = 1, size(faces_of_cell)
            call PetscSectionGetOffset(faceVecSection, faces_of_cell(iface), faceVec_offset, ierr); call CHKERR(ierr)
            xcellVec(cell_offset+iface) = xfaceVec(i1+faceVec_offset)
          enddo
          call DMPlexRestoreCone(celldm, icell, faces_of_cell,ierr); call CHKERR(ierr)
        enddo

        call VecRestoreArrayReadF90(faceVec, xfaceVec, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(cellVec, xcellVec, ierr); call CHKERR(ierr)
        call DMRestoreLocalVector(faceVec_dm, faceVec, ierr); call CHKERR(ierr)

        call DMGetGlobalVector(celldm, global_cellVec, ierr); call CHKERR(ierr)
        call VecSet(global_cellVec, zero, ierr); call CHKERR(ierr)
        call DMLocalToGlobalBegin(celldm, cellVec, ADD_VALUES, global_cellVec, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd  (celldm, cellVec, ADD_VALUES, global_cellVec, ierr); call CHKERR(ierr)

        call DMRestoreLocalVector(celldm, cellVec, ierr); call CHKERR(ierr)

        cellVecname = 'fV2cV_'//trim(faceVecname)
        call PetscObjectSetName(global_cellVec, cellVecname, ierr); call CHKERR(ierr)

        call PetscObjectGetName(global_cellVec, cellVecname, ierr); call CHKERR(ierr)

        call PetscObjectViewFromOptions(global_cellVec, PETSC_NULL_VEC, '-show_'//trim(cellVecname), ierr); call CHKERR(ierr)

        call DMRestoreGlobalVector(celldm, global_cellVec, ierr); call CHKERR(ierr)
        call DMDestroy(celldm, ierr); call CHKERR(ierr)
        if(ldebug.and.myid.eq.0) print *,'facevec2cellvec :: end'//trim(faceVecname)
      end subroutine

      subroutine set_coords(plex, icongrid)
        type(t_icongrid), intent(in) :: icongrid
        type(t_plexgrid), intent(inout) :: plex

        real(ireals), pointer :: coords(:)
        type(tVec)            :: coordinates
        integer(iintegers)    :: dimEmbed, coordSize, voff
        type(tPetscSection)   :: coordSection

        real(ireals) :: cart_coord(3)

        real(ireals), parameter :: sphere_radius= 6371229  ! [m]
        logical :: l_is_spherical_coords

        integer(mpiint) :: ierr
        integer(iintegers) :: i, k, ivertex

        if(.not.allocated(plex%hhl)) stop 'plex_grid::set_coords -> plex%hhl is not allocated. Need to know the height of the grid before I can setup the coordinates!'

        l_is_spherical_coords = any(.not. approx(icongrid%cartesian_z_vertices, zero))

        call DMGetCoordinateDim(plex%dm, dimEmbed, ierr); call CHKERR(ierr)
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
        print *,'Coord Section has size:', coordSize, '(',coordSize/3,' vertices)'

        call VecCreate(PETSC_COMM_SELF, coordinates, ierr); call CHKERR(ierr)
        call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
        call VecSetBlockSize(coordinates, dimEmbed, ierr);call CHKERR(ierr)
        call VecSetType(coordinates, VECSTANDARD, ierr);call CHKERR(ierr)

        call PetscObjectSetName(coordinates, "coordinates", ierr); call CHKERR(ierr)

        call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

        ! set vertices as coordinates
        do k = 1, plex%Nz+1
          do i = 1, icongrid%Nvertices
            ivertex = ivertex_icon_2_plex(icongrid, plex, i, k)
            call PetscSectionGetOffset(coordSection, ivertex, voff, ierr); call CHKERR(ierr)

            cart_coord = [icongrid%cartesian_x_vertices(i), &
                          icongrid%cartesian_y_vertices(i), &
                          icongrid%cartesian_z_vertices(i)]

            if(l_is_spherical_coords) then
              cart_coord = cart_coord * (sphere_radius + plex%hhl(k))
            else
              cart_coord(3) = plex%hhl(k)
            endif
            coords(voff+i1 : voff+dimEmbed) = cart_coord(i1:dimEmbed)
          enddo
        enddo

        call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

        call DMSetCoordinatesLocal(plex%dm, coordinates, ierr);call CHKERR(ierr)
        call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); call CHKERR(ierr)

        call VecDestroy(coordinates, ierr);call CHKERR(ierr)
      end subroutine

    subroutine set_sf_graph(icongrid, plexgrid)
      type(t_icongrid), intent(in) :: icongrid
      type(t_plexgrid), intent(inout) :: plexgrid

      type(tPetscSF) :: sf

      integer(iintegers) :: nroots, nleaves
      integer(iintegers),allocatable :: ilocal_elements(:)
      type(PetscSFNode),allocatable :: iremote_elements(:)
      type(PetscCopyMode),parameter :: localmode=PETSC_COPY_VALUES, remotemode=PETSC_COPY_VALUES

      integer(mpiint) :: myid, numnodes, ierr
      integer(iintegers) :: icell, iface, iedge, ivertex
      integer(iintegers) :: k, ilocal, iparent, iremote, owner, ileaf, jparent

      call mpi_comm_rank(plexgrid%comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(plexgrid%comm, numnodes, ierr); call CHKERR(ierr)

      nleaves = 0
      do k=plexgrid%fStart, plexgrid%pEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel, k, owner, ierr); call CHKERR(ierr)
        !print *,myid,'ownership of plex element',k,owner
        if(owner.ne.myid) nleaves = nleaves + 1
      enddo

      if(ldebug) print *,myid,'remote elements:',nleaves

      allocate(ilocal_elements(nleaves))  ! local indices of elements
      allocate(iremote_elements(nleaves)) ! remote indices of elements

      ileaf = 0
      do icell = plexgrid%cStart, plexgrid%cEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel, icell, owner, ierr); call CHKERR(ierr)
        if(owner.ne.myid) then
          k = plexgrid%zindex(icell)
          ilocal = plexgrid%localiconindex(icell)

          iparent = icongrid%cell_index(ilocal) ! global cell index in parent grid
          jparent = icongrid%remote_cell_index(iparent) ! local icon index @ neighbour

          iremote = icell_icon_2_plex(plexgrid, jparent, k) ! plex index of cell at neighbor
          ileaf = ileaf + 1
          ilocal_elements(ileaf) = icell
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote

          !print *,myid, 'local cell', icell, 'parent', iparent, '@', owner, '->', iremote
        endif
      enddo

      do iface = plexgrid%fStart, plexgrid%fEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
        if(owner.ne.myid) then
          k = plexgrid%zindex(iface)
          ilocal = plexgrid%localiconindex(iface)

          if(plexgrid%ltopfacepos(iface)) then
            iparent = icongrid%cell_index(ilocal) ! global cell index in parent grid
            jparent = icongrid%remote_cell_index(iparent) ! local icon index @ neighbour
            iremote = iface_top_icon_2_plex(plexgrid, jparent, k, icongrid, owner) ! plex index at neighbor
            !if(myid.eq.0) print *,myid,'top bot face:: plex:',iface, 'icon:', ilocal, iparent, jparent, '::', iremote,'@', owner
          else
            iparent = icongrid%edge_index(ilocal) ! global edge index in parent grid
            jparent = icongrid%remote_edge_index(iparent) ! local icon index @ neighbour
            iremote = iface_side_icon_2_plex(icongrid, plexgrid, jparent, k, owner) ! plex index at neighbor

            !if(myid.eq.0) print *,myid,'side face:: plex',iface, 'icon:', ilocal, iparent, jparent, '::', iremote,'@', owner
          endif

          ileaf = ileaf + 1
          ilocal_elements(ileaf) = iface
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          !print *,myid, 'local iface',facepos, iface, 'parent', iparent, jparent, '@', owner, '->', iremote
        endif
      enddo

      do iedge = plexgrid%eStart, plexgrid%eEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel, iedge, owner, ierr); call CHKERR(ierr)
        if(owner.ne.myid) then
          k = plexgrid%zindex(iedge)
          ilocal = plexgrid%localiconindex(iedge)

          if(plexgrid%ltopfacepos(iedge)) then
            iparent = icongrid%edge_index(ilocal) ! global cell index in parent grid
            jparent = icongrid%remote_edge_index(iparent) ! local icon index @ neighbour
            iremote = iedge_top_icon_2_plex(icongrid, plexgrid, jparent, k, owner) ! plex index at neighbor
          else
            iparent = icongrid%vertex_index(ilocal) ! global edge index in parent grid
            jparent = icongrid%remote_vertex_index(iparent) ! local icon index @ neighbour
            iremote = iedge_side_icon_2_plex(icongrid, plexgrid, jparent, k, owner) ! plex index at neighbor
          endif

          ileaf = ileaf + 1
          ilocal_elements(ileaf) = iedge
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          !print *,myid, 'local edge', iedge, 'parent(edges/vertices)', iparent, '@', owner, '->', iremote
        endif
      enddo

      do ivertex = plexgrid%vStart, plexgrid%vEnd-1
        call DMLabelGetValue(plexgrid%ownerlabel    , ivertex, owner, ierr); call CHKERR(ierr)
        if(owner.ne.myid) then
          k = plexgrid%zindex(ivertex)
          ilocal = plexgrid%localiconindex(ivertex)

          iparent = icongrid%vertex_index(ilocal) ! global edge index in parent grid
          jparent = icongrid%remote_vertex_index(iparent) ! local icon index @ neighbour
          iremote = ivertex_icon_2_plex(icongrid, plexgrid, jparent, k, owner) ! plex index at neighbor

          ileaf = ileaf + 1
          ilocal_elements(ileaf) = ivertex
          iremote_elements(ileaf)%rank = owner
          iremote_elements(ileaf)%index = iremote
          !print *,myid, 'local vertex', ivertex, 'parent', iparent, '@', owner, '->', iremote
        endif
      enddo

      if(ileaf.ne.nleaves) stop 'Seems like we forgot some remote elements? ileaf .ne. nleaves'
      do ileaf=1,nleaves
        !print *,myid,'starforest topology',ileaf, ilocal_elements(ileaf),'-->',iremote_elements(ileaf)%index,'@',iremote_elements(ileaf)%rank
      enddo

      nroots = plexgrid%pEnd

      call DMGetPointSF(plexgrid%dm, sf, ierr); call CHKERR(ierr)
      call PetscSFSetGraph(sf, nroots, nleaves, ilocal_elements, localmode, iremote_elements, remotemode, ierr); call CHKERR(ierr)

      call PetscSFSetUp(sf, ierr); call CHKERR(ierr)
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

    subroutine distribute_plexgrid_dm(plex, dm)
      type(t_plexgrid), intent(inout) :: plex
      type(tDM), intent(inout) :: dm
      integer(mpiint) :: numnodes, ierr

      type(tDM)       :: dmdist

      stop 'distribute_plexgrid_dm :: this is not really useful for plexrt,  we setup the topology ourselves.'

      call mpi_comm_size(plex%comm, numnodes, ierr); call CHKERR(ierr)

      call DMPlexSetAdjacencyUseCone(dm, PETSC_TRUE, ierr); call CHKERR(ierr)
      call DMPlexSetAdjacencyUseClosure(dm, PETSC_FALSE, ierr); call CHKERR(ierr)

      call DMPlexDistribute(dm, i0, PETSC_NULL_SF, dmdist, ierr); call CHKERR(ierr)
      if(dmdist.ne.PETSC_NULL_DM) then
        call DMDestroy(dm, ierr); call CHKERR(ierr)
        dm   = dmdist
      endif

      call update_plex_indices(plex)
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_dist_plex", ierr); call CHKERR(ierr)
    end subroutine

  subroutine compute_face_geometry(plex, dm)
    type(t_plexgrid), intent(inout) :: plex
    type(tDM), intent(out), allocatable :: dm

    type(tVec) :: coordinates
    integer(iintegers) :: cStart, cEnd, icell
    integer(iintegers) :: fStart, fEnd, iface
    integer(iintegers) :: eStart, eEnd, iedge
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

    real(ireals) :: AB(3), CD(3), dotp, area
    integer(iintegers) :: iface_up, iface_dn

    integer(mpiint) :: ierr

    if(allocated(dm)) stop 'called compute_face_geometry on an already allocated geom_dm'
    allocate(dm)
    call DMClone(plex%dm, dm, ierr); ; call CHKERR(ierr)
    call DMGetCoordinateDim(dm, Ndim, ierr); call CHKERR(ierr)
    call DMGetCoordinateSection(dm, coordSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(coordSection, PETSC_NULL_SECTION, "-show_dm_coord_section", ierr); call CHKERR(ierr)
    ! Geometry Vec Contains 3 Fields:
    ! 1: 3 dof for centroid on cells and faces
    ! 2: 3 for normal vecs on faces
    ! 3: and 1 cell, face and edge entry for volume, area, length
    call create_plex_section(plex%comm, dm, 'Geometry Section', i3, &
      [i3, i0, i1], [i3, i3, i1], [i0, i0, i1], [i0, i0, i0], geomSection)
    call DMSetDefaultSection(dm, geomSection, ierr); call CHKERR(ierr)
    call PetscSectionDestroy(geomSection, ierr); call CHKERR(ierr)
    call DMGetDefaultSection(dm, geomSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(geomSection, PETSC_NULL_SECTION, "-show_dm_geom_section", ierr); call CHKERR(ierr)


    call DMGetCoordinatesLocal(dm, coordinates, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_dm_coord", ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(dm, i0, cStart, cEnd, ierr); call CHKERR(ierr)  ! cells
    call DMPlexGetHeightStratum(dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! vertices
    call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    call DMGetLocalVector(dm, plex%geomVec,ierr); call CHKERR(ierr)
    call PetscObjectSetName(plex%geomVec, 'geomVec', ierr);call CHKERR(ierr)
    call VecGetArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    ! First Set lengths of edges
    allocate(vertex_coord(Ndim, 2))
    do iedge = eStart, eEnd-1
      call DMPlexGetCone(dm, iedge, vertices, ierr); call CHKERR(ierr)

      ! Get the coordinates of vertices
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:voff+Ndim)
      enddo

      !print *,'edge length:',iedge,'::', distance(vertex_coord(:,1), vertex_coord(:,2))
      call PetscSectionGetOffset(geomSection, iedge, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff) = distance(vertex_coord(:,1), vertex_coord(:,2))

      call DMPlexRestoreCone(dm, iedge, vertices, ierr); call CHKERR(ierr)
    enddo
    deallocate(vertex_coord)

    ! Then define geom info for faces
    do iface = fStart, fEnd-1
      call DMPlexGetConeSize(dm, iface, Nedges, ierr); call CHKERR(ierr)

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

      call DMPlexGetTransitiveClosure(dm, iface, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      !print *,'transclosure', iface,'::',Nedges,'::',transclosure

      vertices = transclosure(3+2*Nedges:size(transclosure):2) ! indices come from: 1 for faceid, Nedge, and then vertices, each with two entries, one for index, and one for orientation
      !print *,'iface',iface,'vertices', vertices

      call DMPlexRestoreTransitiveClosure(dm, iface, PETSC_True, transclosure, ierr); call CHKERR(ierr)

      ! Get the coordinates of vertices
      allocate(vertex_coord(Ndim, Nvertices))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:voff+Ndim)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      !print *,'centroid of face:',iface,'::', sum(vertex_coord,dim=2)/real(Nvertices,ireals)
      call PetscSectionGetOffset(geomSection, iface, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2) / real(Nvertices, kind=ireals)

      ! and use 3 coordinates to compute normal
      ! print *,'normal of face', iface,'::', compute_normal_3d(vertex_coord(1,:),vertex_coord(2,:),vertex_coord(3,:))
      geoms(voff+Ndim+i1: voff+Ndim*2) = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,2),vertex_coord(:,3))

      if(Nvertices.eq.3) then
        geoms(i1+voff+Ndim*2) = triangle_area_by_vertices(vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3))
        !print *,'face triangle area:', geoms(i1+voff+Ndim*2)
      else
        AB = vertex_coord(:,2) - vertex_coord(:,1)
        CD = vertex_coord(:,4) - vertex_coord(:,3)
        dotp = dot_product(AB/norm(AB),CD/norm(CD))
        if(dotp.gt.zero) then
          !print *,'swapping vertex coordinates because dot_product>0', dotp
          call swap(vertex_coord(:,1),vertex_coord(:,2))
        endif
        geoms(i1+voff+Ndim*2) = triangle_area_by_vertices(vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3)) + &
                                triangle_area_by_vertices(vertex_coord(:,3), vertex_coord(:,4), vertex_coord(:,1))
        !print *,'face rectangle area:', geoms(i1+voff+Ndim*2)
      endif

      deallocate(vertex_coord)
    enddo

    ! Last but not least define geom info for cells
    do icell = cStart, cEnd-1
      call DMPlexGetTransitiveClosure(dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      !print *,'cell transclosure:',icell, transclosure
      if(size(transclosure).ne.21*2) then
        print *,'len of transclosure with size', size(transclosure), 'not supported -- is this a wedge?'
        stop('geometry not supported -- is this a wedge?')
      endif

      vertices => vertices6

      vertices = transclosure(size(transclosure)-11:size(transclosure):2)

      allocate(vertex_coord(Ndim, size(vertices)))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:Ndim+voff)
      enddo

      !print *,'centroid of cell:',icell,'::', sum(vertex_coord,dim=2)/size(vertices)
      call PetscSectionGetOffset(geomSection, icell, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2)/size(vertices)

      ! Compute volume of wedges
      iface_up = transclosure(3)
      iface_dn = transclosure(5)

      call PetscSectionGetOffset(geomSection, iface_dn, voff, ierr); call CHKERR(ierr)
      area = geoms(i1+voff+i7)
      call PetscSectionGetOffset(geomSection, iface_dn, voff, ierr); call CHKERR(ierr)
      !print *,'cell top area', area, 'down_area', geoms(i1+voff+i7)
      area = (area + geoms(i1+voff+i7))/2

      call PetscSectionGetOffset(geomSection, icell, voff, ierr); call CHKERR(ierr)
      geoms(i1+voff+Ndim+i1) = area * (plex%hhl(plex%zindex(icell))-plex%hhl(plex%zindex(icell)+i1))
      !print *,'cell volume', geoms(i1+voff+Ndim+i1)

      call DMPlexRestoreTransitiveClosure(dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      deallocate(vertex_coord)
    enddo

    call VecRestoreArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)
    call VecRestoreArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(plex%geomVec, PETSC_NULL_VEC, "-show_dm_geom_vec", ierr); call CHKERR(ierr)
    if(ldebug) call mpi_barrier(plex%comm, ierr); call CHKERR(ierr)
  end subroutine

  subroutine create_plex_section(comm, dm, sectionname, numfields, cdof, fdof, edof, vdof, section, fieldnames)
    MPI_Comm, intent(in) :: comm
    type(tDM), intent(in) :: dm
    character(len=*) :: sectionname
    integer(iintegers), intent(in) :: numfields
    integer(iintegers), intent(in) :: cdof(:), fdof(:), edof(:), vdof(:) ! dim=numfields
    type(tPetscSection), intent(out) :: section
    character(len=*), intent(in), optional :: fieldnames(:)

    integer(iintegers)    :: i, ifield, section_size, sum_cdof, sum_fdof, sum_edof, sum_vdof

    integer(iintegers) :: pStart, pEnd
    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(iintegers) :: eStart, eEnd
    integer(iintegers) :: vStart, vEnd

    integer(mpiint) :: ierr

    ! This is code to manually create a section, e.g. as in DMPlexCreateSection
    sum_cdof = sum(cdof)
    sum_fdof = sum(fdof)
    sum_edof = sum(edof)
    sum_vdof = sum(vdof)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    ! Create Default Section
    call PetscSectionCreate(comm, section, ierr); call CHKERR(ierr)
    call PetscSectionSetNumFields(section, numfields, ierr); call CHKERR(ierr)
    call PetscSectionSetChart(section, pStart, pEnd, ierr); call CHKERR(ierr)

    if(sum_cdof.gt.i0) then
      do i = cStart, cEnd-1
        call PetscSectionSetDof(section, i, sum_cdof, ierr); call CHKERR(ierr)
        do ifield = 1, numfields
          call PetscSectionSetFieldDof(section, i, ifield-1, cdof(ifield), ierr); call CHKERR(ierr)
        enddo
      enddo
    endif
    if(sum_fdof.gt.i0) then
      do i = fStart, fEnd-1
        call PetscSectionSetDof(section, i, sum_fdof, ierr); call CHKERR(ierr)
        do ifield = 1, numfields
          call PetscSectionSetFieldDof(section, i, ifield-1, fdof(ifield), ierr); call CHKERR(ierr)
        enddo
      enddo
    endif
    if(sum_edof.gt.i0) then
      do i = eStart, eEnd-1
        call PetscSectionSetDof(section, i, sum_edof, ierr); call CHKERR(ierr)
        do ifield = 1, numfields
          call PetscSectionSetFieldDof(section, i, ifield-1, edof(ifield), ierr); call CHKERR(ierr)
        enddo
      enddo
    endif
    if(sum_vdof.gt.i0) then
      do i = vStart, vEnd-1
        call PetscSectionSetDof(section, i, sum_vdof, ierr); call CHKERR(ierr)
        do ifield = 1, numfields
          call PetscSectionSetFieldDof(section, i, ifield-1, vdof(ifield), ierr); call CHKERR(ierr)
        enddo
      enddo
    endif

    call PetscSectionSetUp(section, ierr); call CHKERR(ierr)

    call PetscSectionGetStorageSize(section, section_size, ierr);

    call PetscObjectSetName(section, trim(sectionname), ierr);call CHKERR(ierr)
    if(present(fieldnames)) then
      if(size(fieldnames).ne.numfields) stop 'plex_grid::create_plex_section : dim of fieldnames.ne.numfields'
      do i = 1, numfields
        call PetscSectionSetFieldName(section, i-1, fieldnames(i), ierr); call CHKERR(ierr)
      enddo
    endif
    call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, '-show_'//trim(sectionname), ierr); call CHKERR(ierr)
  end subroutine

  subroutine setup_edir_dmplex(plex, dm)
    type(t_plexgrid), intent(inout) :: plex
    type(tDM), allocatable, intent(inout) :: dm
    type(tPetscSection) :: edirSection
    integer(mpiint) :: ierr

    if(allocated(dm)) stop 'called setup_edir_dmplex on an already allocated DM'
    allocate(dm)

    call DMClone(plex%dm, dm, ierr); call CHKERR(ierr)

    call PetscObjectSetName(dm, 'plex_edir', ierr);call CHKERR(ierr)
    call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex_edir", ierr); call CHKERR(ierr)

    call create_plex_section(plex%comm, dm, 'Face Section', i1, [i0], [i1], [i0], [i0], edirSection)  ! Contains 1 dof on each side for direct radiation
    call DMSetDefaultSection(dm, edirSection, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(edirSection, PETSC_NULL_SECTION, '-show_edir_section', ierr); call CHKERR(ierr)
    call PetscSectionDestroy(edirSection, ierr); call CHKERR(ierr)

    !call setup_edir_dmplex_test()
    contains
      subroutine setup_edir_dmplex_test()
        type(tVec) :: lvec, lvec2, gvec

        integer(mpiint) :: myid, i
        real(ireals), pointer :: xv(:), xv2(:)

        call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

        call DMCreateLocalVector(dm, lvec, ierr); call CHKERR(ierr)
        call DMCreateLocalVector(dm, lvec2, ierr); call CHKERR(ierr)
        call DMCreateGlobalVector(dm, gvec, ierr); call CHKERR(ierr)

        call VecSet(lvec, real(myid,kind=ireals), ierr); call CHKERR(ierr)
        call VecSet(lvec2, real(-1,kind=ireals), ierr); call CHKERR(ierr)
        call VecSet(gvec, real(-1,kind=ireals), ierr); call CHKERR(ierr)

        call DMLocalToGlobalBegin(dm,lvec,ADD_VALUES,gvec, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd  (dm,lvec,ADD_VALUES,gvec, ierr); call CHKERR(ierr)

        call DMGlobalToLocalBegin(dm,gvec,INSERT_VALUES,lvec2, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd  (dm,gvec,INSERT_VALUES,lvec2, ierr); call CHKERR(ierr)

        call VecGetArrayReadF90(lvec,xv,ierr)
        call VecGetArrayReadF90(lvec2,xv2,ierr)
        do i=lbound(xv,1), ubound(xv,1)
          print *,myid,'xv/xv2 ::',i, int(xv(i)), int(xv2(i))
        enddo
        call VecRestoreArrayReadF90(lvec2,xv2,ierr)
        call VecRestoreArrayReadF90(lvec,xv,ierr)

        call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_edir_gvec', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(lvec, PETSC_NULL_VEC, '-show_edir_lvec', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(lvec2, PETSC_NULL_VEC, '-show_edir_lvec2', ierr); call CHKERR(ierr)

        call DMRestoreLocalVector(dm, lvec,ierr); call CHKERR(ierr)
        call DMRestoreLocalVector(dm, lvec2,ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(dm, gvec,ierr); call CHKERR(ierr)

        call mpi_barrier(plex%comm, ierr)
        stop 'debug'

      end subroutine
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

    call create_plex_section(plex%comm, dm, 'Absorption Section', i1, [i1], [i0], [i0], [i0], s)  ! Contains 1 dof on each cell
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
    call DMPlexGetHeightStratum(dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    print *,myid,'pStart,End :: ',pStart, pEnd
    print *,myid,'cStart,End :: ',cStart, cEnd
    print *,myid,'fStart,End :: ',fStart, fEnd
    print *,myid,'eStart,End :: ',eStart, eEnd
    print *,myid,'vStart,End :: ',vStart, vEnd

  end subroutine

  subroutine ncvar2d_to_globalvec(plexgrid, filename, varname, vec)
    type(t_plexgrid), intent(in) :: plexgrid
    character(len=*), intent(in) :: filename, varname
    type(tVec), intent(inout) :: vec

    type(tVec) :: local
    type(tVecScatter) :: scatter_context
    real(ireals), pointer :: xloc(:)=>null()
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: ic, icell, icell_global, k, Ncells
    real(ireals), allocatable :: arr(:,:)
    character(len=default_str_len) :: ncgroups(2)

    call mpi_comm_rank(plexgrid%comm, myid, ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,'Loading from nc file: ',trim(filename), ' varname: ', trim(varname)

    call print_dmplex(plexgrid%comm, plexgrid%dm)

    ! Now lets get vectors!
    call DMCreateGlobalVector(plexgrid%dm, vec, ierr); call CHKERR(ierr)
    call PetscObjectSetName(vec, 'VecfromNC_'//trim(varname), ierr);call CHKERR(ierr)

    call VecScatterCreateToZero(vec, scatter_context, local, ierr); call CHKERR(ierr)
    !call AOView(plexgrid%cell_ao, PETSC_VIEWER_STDOUT_WORLD, ierr)

    if(myid.eq.0) then
      ncgroups(1) = trim(filename)
      ncgroups(2) = trim(varname)
      call ncload(ncgroups, arr, ierr); call CHKERR(ierr) ! arr probably has shape ncell, hhl
      Ncells = size(arr, dim=1)

      call VecGetArrayF90(local,xloc,ierr); call CHKERR(ierr)
      do icell=1,Ncells
        icell_global = icell-1
        call AOApplicationToPetsc(plexgrid%cell_ao, i1, icell_global, ierr); call CHKERR(ierr) ! icell_global is petsc index for cells from 0 to Ncells-1

        do k=1,size(arr, dim=2)
          ic = icell_global*size(arr, dim=2) + k
          xloc(ic) = arr(icell, k)
        enddo
      enddo
      call VecRestoreArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
    endif

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: scatter reverse....'
    call VecScatterBegin(scatter_context, local, vec, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
    call VecScatterEnd  (scatter_context, local, vec, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

    call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)
    call VecDestroy(local,ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(vec, PETSC_NULL_VEC, '-show_ncvar_to_globalvec', ierr); call CHKERR(ierr)
    if(ldebug) call mpi_barrier(plexgrid%comm, ierr); call CHKERR(ierr)
  end subroutine

    !> @brief return the dmplex cell index for an icon base grid cell index
    function icell_icon_2_plex(plex, icell, k)
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz) then
          print *,'icell_icon_2_plex :: inp',icell, k
          stop 'icell_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1) then
          print *,'icell_icon_2_plex :: inp',icell, k
          stop 'icell_icon_2_plex :: icon cell index out of range'
        endif
      endif
      icell_icon_2_plex = (k-i1) + (icell-i1)*plex%Nz
    end function

    !> @brief return the dmplex face index for an icongrid index situated at the top of a cell
    function iface_top_icon_2_plex(plex, icell, k, icon, owner)
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      type(t_icongrid), intent(in),optional  :: icon  !< @param[in] icon mesh object, holding info about number of grid cells
      integer(iintegers),intent(in),optional :: owner !< @param[in], optional the mpi rank on which to lookup the index
      integer(iintegers) :: iface_top_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      integer(iintegers) :: offset
      if(present(owner)) then
        offset = icon%parNfaces(owner) * plex%Nz
      else
        offset = plex%offset_faces
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz+1) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner)
          stop 'iface_top_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner)
          stop 'iface_top_icon_2_plex :: icon cell index out of range'
        endif
      endif
      iface_top_icon_2_plex = offset + (k-i1) + (icell-i1)*(plex%Nz+1)
    end function

    !> @brief return the dmplex face index for an icongrid index situated at the side of a cell, i.e. below a certain edge
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
          stop 'iface_side_icon_2_plex :: vertical index k out of range'
        endif
        if(iedge.lt.i1 .or. iedge.gt.Nedges) then
          print *,'iface_side_icon_2_plex :: inp', iedge, k, present(owner), '::', Nfaces, Nedges, offset
          stop 'iface_side_icon_2_plex :: icon edge index out of range'
        endif
      endif
      iface_side_icon_2_plex = offset + (k-i1) + (iedge-i1)*plex%Nz
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

      iedge_top_icon_2_plex = offset + (k-i1) + (iedge-i1)*(plex%Nz+i1)
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

      iedge_side_icon_2_plex = offset + (k-i1) + (ivertex-i1)*plex%Nz
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

      ivertex_icon_2_plex = offset + (k-i1) + (ivertex-i1)*(plex%Nz+i1)
    end function

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
