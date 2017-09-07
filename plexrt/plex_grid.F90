module m_icon_plexgrid
#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload
  use m_helper_functions, only: CHKERR, compute_normal_3d, spherical_2_cartesian, norm, &
    determine_normal_direction, angle_between_two_vec, rad2deg, hit_plane
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, zero, one, &
    default_str_len

  implicit none

  private
  public :: t_plexgrid, read_icon_grid_file, load_plex_from_file, &
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
    Vec :: geomVec

    real(ireals) :: sunpos(3) ! cartesian direction of sun rays in a global reference system

    character(len=8) :: boundary_label = 'boundary'

    ! Index counters on plex:
    integer(iintegers) :: pStart, pEnd ! points
    integer(iintegers) :: cStart, cEnd ! cells
    integer(iintegers) :: fStart, fEnd ! faces
    integer(iintegers) :: eStart, eEnd ! edges
    integer(iintegers) :: vStart, vEnd ! vertices

    integer(iintegers) :: Nfaces2d, Nedges2d, Nvertices2d ! number of entries in base icon grid
    integer(iintegers) :: Nz = 1 ! Number of layers, in 2D set to one

    integer(iintegers), allocatable, dimension(:,:) :: icon_vertex_of_cell, icon_edge_of_cell, icon_edge_vertices
    integer(iintegers), allocatable, dimension(:) :: icon_cell_index, icon_edge_index, icon_vertex_index
    real(ireals), allocatable, dimension(:) :: icon_cartesian_x_vertices, icon_cartesian_y_vertices, icon_cartesian_z_vertices
    integer(iintegers), allocatable, dimension(:) :: icon_cell_sea_land_mask ! "sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cell"


  end type

  contains
    subroutine read_icon_grid_file(fname, plexgrid)
      character(len=*),intent(in) :: fname
      type(t_plexgrid),intent(inout) :: plexgrid
      character(default_str_len) :: varname(2)
      integer(mpiint) :: ierr

      if( allocated(plexgrid%icon_vertex_of_cell) .or. allocated(plexgrid%icon_edge_of_cell) ) then
        print *,'Icon plexgrid already loaded...'
      else
        if (ldebug) print *,'Reading Icon plexgrid File:', trim(fname)
        varname(1) = fname

        varname(2) = 'vertex_of_cell'; call ncload(varname, plexgrid%icon_vertex_of_cell, ierr); call CHKERR(ierr)
        varname(2) = 'edge_of_cell'  ; call ncload(varname, plexgrid%icon_edge_of_cell  , ierr); call CHKERR(ierr)
        varname(2) = 'edge_vertices' ; call ncload(varname, plexgrid%icon_edge_vertices , ierr); call CHKERR(ierr)
        varname(2) = 'cell_index'    ; call ncload(varname, plexgrid%icon_cell_index    , ierr); call CHKERR(ierr)
        varname(2) = 'edge_index'    ; call ncload(varname, plexgrid%icon_edge_index    , ierr); call CHKERR(ierr)
        varname(2) = 'vertex_index'  ; call ncload(varname, plexgrid%icon_vertex_index  , ierr); call CHKERR(ierr)
        varname(2) = 'cell_sea_land_mask'  ; call ncload(varname, plexgrid%icon_cell_sea_land_mask  , ierr); call CHKERR(ierr)
        varname(2) = 'cartesian_x_vertices'; call ncload(varname, plexgrid%icon_cartesian_x_vertices, ierr); call CHKERR(ierr)
        varname(2) = 'cartesian_y_vertices'; call ncload(varname, plexgrid%icon_cartesian_y_vertices, ierr); call CHKERR(ierr)
        varname(2) = 'cartesian_z_vertices'; call ncload(varname, plexgrid%icon_cartesian_z_vertices, ierr); call CHKERR(ierr)
      endif

      if (ldebug) then
        print *,'shape vertex of cell', shape(plexgrid%icon_vertex_of_cell), size(plexgrid%icon_vertex_of_cell),'::', &
          minval(plexgrid%icon_vertex_of_cell), maxval(plexgrid%icon_vertex_of_cell)

        print *,'shape edge of cell', shape(plexgrid%icon_edge_of_cell), size(plexgrid%icon_edge_of_cell)
        print *,'shape cell_index',   shape(plexgrid%icon_cell_index  )
        print *,'shape edge_index',   shape(plexgrid%icon_edge_index  )
        print *,'shape vertex_index', shape(plexgrid%icon_vertex_index)
      endif
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
        call DMPlexDistribute(plex%dm, 0_mpiint, sf, dmdist, ierr); call CHKERR(ierr)
        call DMDestroy(plex%dm, ierr); call CHKERR(ierr)
        plex%dm   = dmdist
      endif
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
        !print *,'cell',icell,'norm of face', faces_of_cell(iface),'::',face_normal

        ! Then determine if the face is src(in line with the sun vec) or if it is destination(contra sun direction)
        angle_to_sun = rad2deg(angle_between_two_vec(face_normal, plex%sunpos))

        if(angle_to_sun.lt.90 .or. angle_to_sun.ge.270) then
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

    type(tPetscSection) :: s

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(mpiint) :: ierr

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: iface, irow, icol, icell, isrc, idst
    type(PetscInt) :: ifacepos_target, ifacepos_src

    real(ireals) :: coeff(1), angle_to_sun

    type(tDMLabel) :: faceposlabel, zindexlabel

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    integer(iintegers) :: geom_offset

    real(ireals) :: cell_center(3), face_normal(3), face_center(3)
    real(ireals) :: ray_loc(3), ray_dir(3), distance, min_distance
    integer(iintegers) :: min_distance_idst

    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)

    call DMGetDefaultSection(plex%edir_dm, s, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(plex%edir_dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(plex%edir_dm, i1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMGetLabel(plex%edir_dm, "Face Position", faceposlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(plex%edir_dm, "Vertical Index", zindexlabel, ierr); CHKERRQ(ierr)

    call DMCreateMatrix(plex%edir_dm, A, ierr); CHKERRQ(ierr)

    do icell = cStart, cEnd-1
      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); CHKERRQ(ierr) ! Get Faces of cell
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); CHKERRQ(ierr)
        cell_center = geoms(1+geom_offset:3+geom_offset)

        call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); CHKERRQ(ierr)
        face_center = geoms(1+geom_offset: 3+geom_offset)
        face_normal = geoms(4+geom_offset: 6+geom_offset)

        ! Determine the inward normal vec for the face
        face_normal = face_normal * determine_normal_direction(face_normal, face_center, cell_center)
        !print *,'cell',icell,'norm of face', faces_of_cell(iface),'::',face_normal

        ! Then determine if the face is src(in line with the sun vec) or if it is destination(contra sun direction)
        angle_to_sun = rad2deg(angle_between_two_vec(face_normal, plex%sunpos))
        ! print *,'angle to sun', int(angle_to_sun)

        if(angle_to_sun.lt.90 .or. angle_to_sun.ge.270) then
          lsrc(iface) = .True.
        else
          lsrc(iface) = .False.
        endif
      enddo

      !print *,'cell', icell, 'face', faces_of_cell, ' is src? ',lsrc
      !print *,'cell', icell, 'face', faces_of_cell-fStart, ' is src? ',lsrc

      do isrc = 1, size(lsrc)
        if (.not.lsrc(isrc)) cycle ! this is not a destination

        ! Now we have to find the contributions of all sources to each of the destination:
        min_distance = huge(min_distance)

        do idst = 1, size(lsrc)
          if (lsrc(idst)) cycle ! this is not a destination

          ! at this point we have a src face and a destination face
          ! lets find out if we can cast a ray from the center of the src and hit the destination
          call PetscSectionGetOffset(geomSection, faces_of_cell(isrc), geom_offset, ierr); CHKERRQ(ierr)
          ray_loc = geoms(1+geom_offset: 3+geom_offset) ! center of src face
          ray_dir = plex%sunpos

          call PetscSectionGetOffset(geomSection, faces_of_cell(idst), geom_offset, ierr); CHKERRQ(ierr)
          face_center = geoms(1+geom_offset: 3+geom_offset) ! dst face which we try to hit
          face_normal = geoms(4+geom_offset: 6+geom_offset)

          distance = hit_plane(ray_loc, ray_dir, face_center, face_normal)
          if(distance.le.min_distance) then
            min_distance = distance
            min_distance_idst = idst
          endif


          ! call DMLabelGetValue(faceposlabel, faces_of_cell(idst), ifacepos_target, ierr); CHKERRQ(ierr)
          ! call DMLabelGetValue(faceposlabel, faces_of_cell(isrc), ifacepos_src, ierr); CHKERRQ(ierr)
          ! if(ifacepos_target.eq.TOP_BOT_FACE .and. ifacepos_src.eq.TOP_BOT_FACE) then
          !   coeff = .5
          ! else if (ifacepos_target.eq.SIDE_FACE .and. ifacepos_src.eq.TOP_BOT_FACE) then
          !   coeff = .1
          ! else if (ifacepos_target.eq.TOP_BOT_FACE .and. ifacepos_src.eq.SIDE_FACE) then
          !   coeff = .1
          ! else
          !   coeff = zero
          ! endif

          ! call PetscSectionGetOffset(s, faces_of_cell(idst), irow, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
          ! call PetscSectionGetOffset(s, faces_of_cell(isrc), icol, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
          ! !print *,idst, isrc, 'setting val @',faces_of_cell(idst), faces_of_cell(isrc),'offsets',irow,icol

          ! call MatSetValue(A, irow, icol, -coeff, INSERT_VALUES, ierr); CHKERRQ(ierr)
        enddo

         call PetscSectionGetOffset(s, faces_of_cell(isrc), icol, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
         call PetscSectionGetOffset(s, faces_of_cell(min_distance_idst), irow, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
         coeff = .9 ! exp(- min_distance * 1e-4)
         print *,'distance', min_distance_idst, min_distance, coeff

         call MatSetValue(A, irow, icol, -coeff, INSERT_VALUES, ierr); CHKERRQ(ierr)
      enddo

      ! Set diagonal entries
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(s, faces_of_cell(iface), irow, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
        call PetscSectionGetOffset(s, faces_of_cell(iface), icol, ierr); CHKERRQ(ierr) ! this is the offset of the neighboring faces
        call MatSetValue(A, irow, icol, one, INSERT_VALUES, ierr); CHKERRQ(ierr)
      enddo

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); CHKERRQ(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_VEC, '-show_A', ierr); CHKERRQ(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
  end subroutine

end module
