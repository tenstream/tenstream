module m_plex_grid
#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload
  use m_helper_functions, only: CHKERR, compute_normal_3d, spherical_2_cartesian, norm, cross_3d, &
    determine_normal_direction, angle_between_two_vec, rad2deg, deg2rad, hit_plane, &
    rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world, vec_proj_on_plane, &
    get_arg, imp_bcast
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, zero, one, pi, &
    default_str_len
  use m_icon_grid, only : t_icongrid

  implicit none

  private
  public :: t_plexgrid, load_plex_from_file, create_plex_from_local_icongrid, &
    icell_icon_2_plex, icell_plex_2_icon, update_plex_indices, &
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

    integer(iintegers) :: Nz = 1 ! Number of layers, in 2D set to one
  end type

  contains

    subroutine create_plex_from_local_icongrid(comm, Nz, icongrid, local_icongrid, plex)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nz
      type(t_icongrid), allocatable, intent(in) :: icongrid, local_icongrid
      type(t_plexgrid), allocatable, intent(inout) :: plex

      integer(iintegers) :: chartsize, Ncells, Nfaces, Nedges, Nvertices
      integer(iintegers) :: offset_faces, offset_edges, offset_vertices
      integer(iintegers) :: offset_faces_sides, offset_edges_vertical
      integer(iintegers) :: edge3(3), edge4(4), faces(5), vert2(2)

      type(tDMLabel) :: faceposlabel, zindexlabel, TOAlabel

      integer(iintegers) :: i, k, icell, iedge, ivertex
      integer(iintegers) :: depth
      integer(mpiint) :: myid, numnodes, ierr

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(allocated(plex)) stop 'create_plex_from_local_icongrid :: plex should not be allocated!'
      if(.not.allocated(icongrid)) stop 'create_plex_from_local_icongrid :: icongrid should be allocated!'
      if(.not.allocated(local_icongrid)) stop 'create_plex_from_local_icongrid :: local_icongrid should be allocated!'

      allocate(plex)
      plex%comm = comm

      plex%Nz = Nz

      Ncells    = local_icongrid%Nfaces * plex%Nz
      Nfaces    = local_icongrid%Nfaces * (plex%Nz+1) + local_icongrid%Nedges * plex%Nz
      Nedges    = local_icongrid%Nedges * (plex%Nz+1) + local_icongrid%Nvertices * plex%Nz
      Nvertices = local_icongrid%Nvertices * (plex%Nz+i1)

      chartsize = Ncells + Nfaces + Nedges + Nvertices

      allocate(plex%dm)
      call DMPlexCreate(plex%comm, plex%dm, ierr); call CHKERR(ierr)
      call DMSetDimension(plex%dm, i3, ierr); call CHKERR(ierr)

      call DMPlexSetChart(plex%dm, i0, chartsize, ierr); call CHKERR(ierr)

      offset_faces = Ncells
      offset_faces_sides = offset_faces + (plex%Nz+1)*local_icongrid%Nfaces
      offset_edges = Ncells + Nfaces
      offset_edges_vertical = offset_edges + local_icongrid%Nedges * (plex%Nz+i1)
      offset_vertices = Ncells + Nfaces + Nedges
      print *,myid,'offsets faces:', offset_faces, offset_faces_sides
      print *,myid,'offsets edges:', offset_edges, offset_edges_vertical
      print *,myid,'offsets verte:', offset_vertices
      print *,myid,'Chartsize:', chartsize

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
        do i = 1, local_icongrid%Nfaces
          icell = i !local_icongrid%cell_index(i)
          edge3 = local_icongrid%edge_of_cell(i,:)

          faces(1) = offset_faces + local_icongrid%Nfaces*(k-1) + icell-i1  ! top face
          faces(2) = offset_faces + local_icongrid%Nfaces*k + icell-i1      ! bot face

          faces(3:5) = offset_faces_sides + (edge3-i1) + (k-1)*local_icongrid%Nedges

          call DMPlexSetCone(plex%dm, icell_icon_2_plex(local_icongrid, plex, icell, k), faces, ierr); call CHKERR(ierr)

          call DMLabelSetValue(zindexlabel, icell, k, ierr); call CHKERR(ierr)
          print *,myid,'Setting cell indices', i, icell, 'edge3', edge3, 'faces', faces
        enddo
      enddo

      ! set edges of top/bot faces
      do k = 1, plex%Nz+1 ! levels
        do i = 1, local_icongrid%Nfaces
          icell = i !local_icongrid%cell_index(i)
          edge3 = offset_edges + local_icongrid%edge_of_cell(icell,:)-i1 + local_icongrid%Nedges*(k-i1)

          call DMPlexSetCone(plex%dm, offset_faces + local_icongrid%Nfaces*(k-1) + icell -i1, edge3, ierr); call CHKERR(ierr)
          !print *,'edges @ horizontal faces',icell,':',edge3,'petsc:', offset_faces + Nfaces*(k-1) + icell -i1, '::', edge3
          call DMLabelSetValue(faceposlabel, offset_faces + local_icongrid%Nfaces*(k-1) + icell -i1, TOP_BOT_FACE, ierr); call CHKERR(ierr)
          if (k.eq.i1) then
            call DMLabelSetValue(TOAlabel, offset_faces + local_icongrid%Nfaces*(k-1) + icell -i1, i1, ierr); call CHKERR(ierr)
          endif
        enddo
      enddo

      ! set edges of vertical faces
      do k = 1, plex%Nz ! layers
        do i = 1, local_icongrid%Nedges
          iedge = i !local_icongrid%edge_index(i)
          edge4(1) = offset_edges + iedge-i1 + local_icongrid%Nedges*(k-i1)
          edge4(2) = offset_edges + iedge-i1 + local_icongrid%Nedges*(k)

          vert2 = local_icongrid%edge_vertices(iedge,:)-i1 + local_icongrid%Nvertices*(k-i1)
          edge4(3:4) = offset_edges_vertical + vert2

          call DMPlexSetCone(plex%dm, offset_faces_sides + iedge-i1 + local_icongrid%Nedges*(k-i1),edge4, ierr); call CHKERR(ierr)
          print *,'edges @ vertical faces',iedge,':',edge4,'petsc:', offset_faces_sides + iedge-i1 + Nedges*(k-i1), '::', edge4
          call DMLabelSetValue(faceposlabel, offset_faces_sides + iedge-i1 + local_icongrid%Nedges*(k-i1), SIDE_FACE, ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      ! and then set the two vertices of edges in each level
      do k = 1, plex%Nz+1 ! levels
        do i = 1, local_icongrid%Nedges
          iedge = i !local_icongrid%edge_index(i)
          vert2 = offset_vertices + local_icongrid%edge_vertices(iedge,:) + local_icongrid%Nvertices*(k-i1)

          print *,'vertices @ edge2d',iedge,':',vert2,'petsc:',offset_edges + iedge-i1 + local_icongrid%Nedges*(k-i1),'::', vert2-i1
          call DMPlexSetCone(plex%dm, offset_edges + iedge-i1 + local_icongrid%Nedges*(k-i1), vert2-i1, ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      ! and then set the two vertices of edges in each layer
      do k = 1, plex%Nz ! layer
        do i = 1, local_icongrid%Nvertices
          ivertex = i !local_icongrid%vertex_index(i)
          vert2(1) = offset_vertices + ivertex + local_icongrid%Nvertices*(k-i1)
          vert2(2) = offset_vertices + ivertex + local_icongrid%Nvertices*(k)

          print *,'vertices @ edge_vertical',ivertex,':',vert2,'petsc:',offset_edges_vertical + ivertex-i1 + Nvertices*(k-i1),':', vert2 - i1
          call DMPlexSetCone(plex%dm, offset_edges_vertical + ivertex-i1 + local_icongrid%Nvertices*(k-i1), vert2 - i1, ierr); call CHKERR(ierr)
        enddo
      enddo

      if(ldebug) call mpi_barrier(comm, ierr)

      print *,'Symmetrize'
      call DMPlexSymmetrize(plex%dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(plex%dm, ierr); call CHKERR(ierr)

      if(ldebug) then
        call DMPlexGetDepth(plex%dm, depth, ierr); CHKERRQ(ierr)
        print *,'Depth of Stratum:', depth
      endif

      call update_plex_indices(plex)

      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
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
    call DMClone(plex%dm, plex%geom_dm, ierr); ; CHKERRQ(ierr)
    call DMGetCoordinateDim(plex%geom_dm, Ndim, ierr); CHKERRQ(ierr)
    call DMGetCoordinateSection(plex%geom_dm, coordSection, ierr); CHKERRQ(ierr)
    call PetscObjectViewFromOptions(coordSection, PETSC_NULL_SECTION, "-show_dm_coord_section", ierr); CHKERRQ(ierr)

    call create_plex_section(plex%comm, plex%geom_dm, 'Geometry Section', i1*3, i1*6, i0, i0, geomSection)  ! Contains 3 dof for centroid on cells and faces plus 3 for normal vecs on faces
    call DMSetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call PetscSectionDestroy(geomSection, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call PetscObjectViewFromOptions(geomSection, PETSC_NULL_SECTION, "-show_dm_geom_section", ierr); CHKERRQ(ierr)


    call DMGetCoordinatesLocal(plex%geom_dm, coordinates, ierr); CHKERRQ(ierr)
    call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_dm_coord", ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(coordinates, coords, ierr); CHKERRQ(ierr)

    call DMPlexGetHeightStratum(plex%geom_dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr)  ! cells
    call DMPlexGetHeightStratum(plex%geom_dm, i1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMPlexGetDepthStratum (plex%geom_dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices

    call DMGetLocalVector(plex%geom_dm, plex%geomVec,ierr); CHKERRQ(ierr)
    call PetscObjectSetName(plex%geomVec, 'geomVec', ierr);CHKERRQ(ierr)
    call VecGetArrayF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

    do icell = cStart, cEnd-1
      call DMPlexGetTransitiveClosure(plex%geom_dm, icell, PETSC_TRUE, transclosure, ierr); CHKERRQ(ierr)
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
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); CHKERRQ(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:Ndim+voff)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      ! print *,'centroid of cell:',icell,'::', sum(vertex_coord,dim=2)/size(vertices)
      call PetscSectionGetOffset(geomSection, icell, voff, ierr); CHKERRQ(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2)/size(vertices)

      call DMPlexRestoreTransitiveClosure(plex%geom_dm, icell, PETSC_TRUE, transclosure, ierr); CHKERRQ(ierr)
      deallocate(vertex_coord)
    enddo

    do iface = fStart, fEnd-1
      call DMPlexGetConeSize(plex%geom_dm, iface, Nedges, ierr); CHKERRQ(ierr)

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

      call DMPlexGetTransitiveClosure(plex%geom_dm, iface, PETSC_TRUE, transclosure, ierr); CHKERRQ(ierr)
      !print *,'transclosure', iface,'::',Nedges,'::',transclosure

      vertices = transclosure(3+2*Nedges:size(transclosure):2) ! indices come from: 1 for faceid, Nedge, and then vertices, each with two entries, one for index, and one for orientation
      !print *,'iface',iface,'vertices', vertices

      call DMPlexRestoreTransitiveClosure(plex%geom_dm, iface, PETSC_True, transclosure, ierr); CHKERRQ(ierr)

      ! Get the coordinates of vertices
      allocate(vertex_coord(Ndim, Nvertices))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); CHKERRQ(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:voff+Ndim)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      !print *,'centroid of face:',iface,'::', sum(vertex_coord,dim=1)/Nvertices
      call PetscSectionGetOffset(geomSection, iface, voff, ierr); CHKERRQ(ierr)
      geoms(i1+voff:voff+Ndim) = sum(vertex_coord,dim=2)/Nvertices

      ! and use 3 coordinates to compute normal
      ! print *,'normal of face', iface,'::', compute_normal_3d(vertex_coord(1,:),vertex_coord(2,:),vertex_coord(3,:))
      geoms(voff+Ndim+i1: voff+Ndim*2) = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,2),vertex_coord(:,3))

      deallocate(vertex_coord)
    enddo

    call VecRestoreArrayReadF90(coordinates, coords, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

    call PetscObjectViewFromOptions(plex%geomVec, PETSC_NULL_VEC, "-show_dm_geom_vec", ierr); CHKERRQ(ierr)
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

    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices

    ! Create Default Section
    call PetscSectionCreate(comm, section, ierr); CHKERRQ(ierr)
    call PetscSectionSetNumFields(section, 1, ierr); CHKERRQ(ierr)
    call PetscSectionSetChart(section, pStart, pEnd, ierr); CHKERRQ(ierr)

    do i = cStart, cEnd-1
      call PetscSectionSetDof(section, i, cdof, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(section, i, 0, cdof, ierr); CHKERRQ(ierr)
    enddo
    do i = fStart, fEnd-1
      call PetscSectionSetDof(section, i, fdof, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(section, i, 0, fdof, ierr); CHKERRQ(ierr)
    enddo
    do i = eStart, eEnd-1
      call PetscSectionSetDof(section, i, edof, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(section, i, 0, edof, ierr); CHKERRQ(ierr)
    enddo
    do i = vStart, vEnd-1
      call PetscSectionSetDof(section, i, vdof, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(section, i, 0, vdof, ierr); CHKERRQ(ierr)
    enddo

    call PetscSectionSetUp(section, ierr); CHKERRQ(ierr)

    call PetscSectionGetStorageSize(section, section_size, ierr);

    call PetscObjectSetName(section, trim(sectionname), ierr);CHKERRQ(ierr)
    call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, '-show_'//trim(sectionname), ierr); CHKERRQ(ierr)
  end subroutine

  subroutine setup_edir_dmplex(plex, dm)
    type(t_plexgrid), intent(inout) :: plex
    type(tDM), allocatable :: dm
    type(tPetscSection) :: edirSection
    integer(mpiint) :: ierr

    if(allocated(dm)) stop 'called setup_edir_dmplex on an already allocated DM'
    allocate(dm)

    call DMClone(plex%dm, dm, ierr); call CHKERR(ierr)

    call PetscObjectSetName(dm, 'plex_edir', ierr);CHKERRQ(ierr)
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

    call PetscObjectSetName(dm, 'plex_abso', ierr);CHKERRQ(ierr)
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

    call mpi_comm_rank(comm, myid, ierr); CHKERRQ(ierr)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices

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

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

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
      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); CHKERRQ(ierr)

      call DMPlexGetCone(plex%abso_dm, icell, faces_of_cell, ierr); CHKERRQ(ierr) ! Get Faces of cell
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(edir_section, faces_of_cell(iface), edir_offset, ierr); CHKERRQ(ierr)

        call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); CHKERRQ(ierr)
        cell_center = geoms(1+geom_offset:3+geom_offset)

        call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); CHKERRQ(ierr)
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
      call DMPlexRestoreCone(plex%abso_dm, icell, faces_of_cell,ierr); CHKERRQ(ierr)
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

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

    call DMGetDefaultSection(plex%edir_dm, sec, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(plex%edir_dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(plex%edir_dm, i1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMGetLabel(plex%edir_dm, "Face Position", faceposlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(plex%edir_dm, "Vertical Index", zindexlabel, ierr); CHKERRQ(ierr)

    call DMCreateMatrix(plex%edir_dm, A, ierr); CHKERRQ(ierr)

    do icell = cStart, cEnd-1

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); CHKERRQ(ierr) ! Get Faces of cell

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
            call PetscSectionGetOffset(sec, faces_of_cell(idst), irow, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
            call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces

            call MatSetValue(A, irow, icol, -dir2dir(idst, iface), INSERT_VALUES, ierr); CHKERRQ(ierr)
          enddo

        endif

      enddo

      ! Set diagonal entries
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(sec, faces_of_cell(iface), irow, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
        call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
        call MatSetValue(A, irow, icol, one, INSERT_VALUES, ierr); CHKERRQ(ierr)
      enddo

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); CHKERRQ(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_VEC, '-show_A', ierr); CHKERRQ(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

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

    call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); CHKERRQ(ierr)
    cell_center = geoms(1+geom_offset:3+geom_offset)

    do iface = 1, size(faces_of_cell)
      call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); CHKERRQ(ierr)
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
      call DMLabelGetValue(faceposlabel, faces_of_cell(iface), ifacepos, ierr); CHKERRQ(ierr)
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

    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(1)), izindex(1), ierr); CHKERRQ(ierr)
    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(2)), izindex(2), ierr); CHKERRQ(ierr)
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
    function icell_icon_2_plex(icon, plex, icell, k)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nz) stop 'icell_icon_2_plex :: vertical index k out of range'
        if(icell.lt.i1 .or. icell.gt.icon%Nfaces) stop 'icell_icon_2_plex :: icon cell index out of range'
      endif
      icell_icon_2_plex = (k-i1)*icon%Nfaces + icell - i1
    end function

    !> @brief return the vertical and icon base grid cell index for a given dmplex cell index.
    function icell_plex_2_icon(icon, plex, icell)
      type(t_icongrid), intent(in) :: icon      !< @param[in] icon mesh object, holding info about number of grid cells
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 0 up to cEnd (size of DMPlex cells)
      integer(iintegers) :: icell_plex_2_icon(2)!< @param[out] return icon cell index for base grid(starting with 1, to Nfaces) and vertical index(1 at top of domain)
      if(ldebug) then
        if(icell.lt.plex%cStart .or. icell.ge.plex%cEnd) stop 'icell_plex_2_icon, dmplex cell index out of range'
      endif
      icell_plex_2_icon(1) = modulo(icell, icon%Nfaces) + i1
      icell_plex_2_icon(2) = icell / icon%Nfaces + i1
      if(ldebug) then
        if(icell_plex_2_icon(2).lt.i1 .or. icell_plex_2_icon(2).gt.plex%Nz) stop 'icell_icon_2_plex :: vertical index k out of range'
        if(icell_plex_2_icon(1).lt.i1 .or. icell_plex_2_icon(1).gt.icon%Nfaces) stop 'icell_icon_2_plex :: icon cell index out of range'
      endif

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
