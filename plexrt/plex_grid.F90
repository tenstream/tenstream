module m_plex_grid
#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload

  use m_helper_functions, only: CHKERR, compute_normal_3d, approx, strF2C, distance, &
    triangle_area_by_vertices, swap, determine_normal_direction, &
    vec_proj_on_plane, cross_3d, rad2deg, is_between, &
    angle_between_two_vec, angle_between_two_normed_vec, &
    resize_arr, get_arg, normalize_vec, &
    imp_bcast, itoa, ftoa, imp_allreduce_max, &
    rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world

  use m_data_parameters, only : ireals, irealLUT, ireal_params, &
    iintegers, mpiint, zero, one, pi, &
    i0, i1, i2, i3, i4, i5, i6, i7, i8, default_str_len

  use m_icon_grid, only : t_icongrid, ICONULL

  use m_tenstream_options, only: read_commandline_options, twostr_ratio
  use m_LUT_param_phi, only: param_phi_param_theta_from_phi_and_theta_withnormals, &
    LUT_wedge_aspect_zx
  use m_optprop_parameters, only: param_eps

  implicit none

  private
  public :: t_plexgrid, &
    icell_icon_2_plex, iface_top_icon_2_plex, update_plex_indices, &
    compute_face_geometry, &
    setup_cell1_dmplex, setup_edir_dmplex, setup_ediff_dmplex, setup_abso_dmplex, &
    print_dmplex, ncvar2d_to_globalvec, facevec2cellvec, &
    orient_face_normals_along_sundir, compute_wedge_orientation, is_solar_src, &
    compute_local_wedge_ordering, compute_local_vertex_coordinates, &
    get_inward_face_normal, create_plex_section, setup_plexgrid, &
    get_consecutive_vertical_cell_idx, get_top_bot_face_of_cell, gen_test_mat, &
    TOAFACE, BOTFACE, SIDEFACE, INNERSIDEFACE, destroy_plexgrid, &
    determine_diff_incoming_outgoing_offsets, get_normal_of_first_TOA_face, &
    interpolate_horizontal_face_var_onto_vertices, get_horizontal_faces_around_vertex, &
    atm_dz_to_vertex_heights, dmplex_set_new_section

  logical, parameter :: ldebug=.False.

  type :: t_plexgrid
    integer(mpiint) :: comm

    type(tDM), allocatable :: dm2d
    type(tDM), allocatable :: dm
    type(tDM), allocatable :: cell1_dm
    type(tDM), allocatable :: horizface1_dm
    type(tDM), allocatable :: edir_dm
    type(tDM), allocatable :: ediff_dm
    type(tDM), allocatable :: abso_dm
    type(tDM), allocatable :: geom_dm
    type(tDM), allocatable :: wedge_orientation_dm
    type(tDM), allocatable :: srfc_boundary_dm
    type(tDM), allocatable :: rayli_dm
    type(tDM), allocatable :: nca_dm
    type(tVec), allocatable :: geomVec ! see compute_face_geometry for details
    type(tVec), allocatable :: wedge_orientation ! see compute_wedge_orientation

    ! Index counters on plex:
    integer(iintegers) :: pStart, pEnd ! points
    integer(iintegers) :: cStart, cEnd ! cells
    integer(iintegers) :: fStart, fEnd ! faces
    integer(iintegers) :: eStart, eEnd ! edges
    integer(iintegers) :: vStart, vEnd ! vertices

    integer(iintegers) :: Nlay = 1 ! Number of layers
    integer(iintegers) :: Ncells   ! Number of cells in DMPlex (wedge form)
    integer(iintegers) :: Nfaces   ! Number of faces (for each cell 2 for top/bottom and 3 sideward-faces)
    integer(iintegers) :: Nedges   ! Number of edges (each cell has 3 at top/bottom faces plus 3 vertically)
    integer(iintegers) :: Nverts   ! Number of vertices (each cell has 3 at top/bottom faces)

    ! In consecutive numbering in Plex chart we first have the cells until offset_faces-1
    integer(iintegers) :: offset_faces           ! then the top/bot faces
    integer(iintegers) :: offset_faces_sides     ! and next the side faces
    integer(iintegers) :: offset_edges
    integer(iintegers) :: offset_edges_vertical
    integer(iintegers) :: offset_vertices

    logical,allocatable :: ltopfacepos(:)                ! TOP_BOT_FACE or SIDE_FACE of faces and edges, fStart..eEnd-1
    integer(iintegers),allocatable :: zindex(:)          ! vertical layer / level of cells/faces/edges/vertices , pStart..pEnd-1
    real(ireals), allocatable :: hhl1d(:)                ! 1D vertical height levels which were used in 2D_to_3D extrusion
    integer(iintegers),allocatable :: localiconindex(:)  ! local index of face, edge, vertex on icongrid, pStart, pEnd-1, i.e. on each rank from 1..Ncell, 1..Nedges etc..
    integer(iintegers),allocatable :: globaliconindex(:) ! global index of face, edge, vertex on icongrid, pStart, pEnd-1, i.e. for each rank has the indices of the global icon grid as it is read from nc

    type(tDMLabel), allocatable :: boundarylabel        ! 1 if boundary of local mesh
    type(tDMLabel), allocatable :: domainboundarylabel  ! TOAFACE if top, SIDEFACE if side face, BOTFACE if bot face, -1 otherwise
    type(tDMLabel), allocatable :: ownerlabel           ! rank that posesses this element
  end type

  integer(iintegers), parameter :: TOAFACE=1, SIDEFACE=2, BOTFACE=3, INNERSIDEFACE=4

  contains

    subroutine destroy_plexgrid(plex)
      type(t_plexgrid), intent(inout) :: plex
      integer(mpiint) :: ierr

      call dealloc_dmlabel(plex%boundarylabel)
      call dealloc_dmlabel(plex%domainboundarylabel)
      call dealloc_dmlabel(plex%ownerlabel)

      call dealloc_vec(plex%geomVec)
      call dealloc_vec(plex%wedge_orientation)

      call dealloc_dm(plex%dm)
      call dealloc_dm(plex%cell1_dm)
      call dealloc_dm(plex%horizface1_dm)
      call dealloc_dm(plex%edir_dm)
      call dealloc_dm(plex%ediff_dm)
      call dealloc_dm(plex%abso_dm)
      call dealloc_dm(plex%geom_dm)
      call dealloc_dm(plex%wedge_orientation_dm)
      call dealloc_dm(plex%srfc_boundary_dm)
      call dealloc_dm(plex%rayli_dm)
      call dealloc_dm(plex%nca_dm)

      plex%pStart = -1; plex%pEnd = -1
      plex%cStart = -1; plex%cEnd = -1
      plex%fStart = -1; plex%fEnd = -1
      plex%eStart = -1; plex%eEnd = -1
      plex%vStart = -1; plex%vEnd = -1
      plex%Nlay = -1
      plex%Ncells = -1
      plex%Nverts = -1
      plex%Nedges = -1
      plex%Nverts = -1

      plex%offset_faces = -1
      plex%offset_faces_sides = -1
      plex%offset_edges = -1
      plex%offset_edges_vertical = -1
      plex%offset_vertices = -1

      if(allocated(plex%ltopfacepos    )) deallocate(plex%ltopfacepos    )
      if(allocated(plex%zindex         )) deallocate(plex%zindex         )
      if(allocated(plex%localiconindex )) deallocate(plex%localiconindex )
      if(allocated(plex%globaliconindex)) deallocate(plex%globaliconindex)

      plex%comm = -1
      contains
        subroutine dealloc_vec(vec)
          type(tVec), allocatable, intent(inout) :: vec
          if(allocated(vec)) then
            call VecDestroy(vec, ierr); call CHKERR(ierr)
            deallocate(vec)
          endif
        end subroutine
        subroutine dealloc_dmlabel(label)
          type(tDMLabel), allocatable, intent(inout) :: label
          if(allocated(label)) then
            ! call DMLabelDestroy(label, ierr); call CHKERR(ierr) ! dont destroy the label here, DMDestroy takes care of it
            deallocate(label)
          endif
        end subroutine
        subroutine dealloc_dm(dm)
          type(tDM), allocatable, intent(inout) :: dm
          if(allocated(dm)) then
            call DMDestroy(dm, ierr); call CHKERR(ierr)
            deallocate(dm)
          endif
        end subroutine
    end subroutine

    subroutine plex_set_ltopfacepos(dm, ltopfacepos)
      type(tDM), intent(in) :: dm
      logical, allocatable, intent(inout) :: ltopfacepos(:) ! TOP_BOT_FACE or SIDE_FACE of faces and edges, fStart..eEnd-1
      integer(iintegers) :: fStart, fEnd, eStart, eEnd
      integer(iintegers) :: iface
      integer(iintegers), pointer :: edges_of_face(:)
      integer(mpiint) :: ierr

      if(allocated(ltopfacepos)) call CHKERR(1_mpiint, 'Dont call plex_set_ltopfacepos on already allocated array!')
      if(ldebug) print *,'plex_set_ltopfacepos... start'

      call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr)

      allocate(ltopfacepos(fStart:eEnd-1), source=.False.)

      do iface = fStart, fEnd-1
        call DMPlexGetCone(dm, iface, edges_of_face, ierr); call CHKERR(ierr)
        select case (size(edges_of_face, kind=iintegers))
        case (i3)
          ltopfacepos(iface) = .True.
          ltopfacepos(edges_of_face) = .True.

        case (i4)
          ltopfacepos(iface) = .False.
          ! side edges remain .false.(from allocation source)

        case default
          call CHKERR(int(size(edges_of_face), mpiint), &
            'Dont know this type of face element with '//itoa(size(edges_of_face))//' edges')
        end select
        call DMPlexRestoreCone(dm, iface, edges_of_face, ierr); call CHKERR(ierr)
      enddo

      if(ldebug) print *,'plex_set_ltopfacepos... end'
    end subroutine

    subroutine setup_plexgrid(dm2d, dm3d, Nlay, zindex, plex, hhl)
      type(tDM), intent(in) :: dm2d, dm3d
      integer(iintegers), intent(in) :: Nlay, zindex(:)
      type(t_plexgrid), allocatable, intent(inout) :: plex
      real(ireals), optional :: hhl(:)

      integer(iintegers) :: pStart, pEnd
      integer(mpiint) :: ierr


      if(allocated(plex)) call CHKERR(1_mpiint, 'Dont call setup_plexgrid on already allocated object')
      allocate(plex)

      call PetscObjectGetComm(dm3d, plex%comm, ierr); call CHKERR(ierr)
      call read_commandline_options(plex%comm)

      allocate(plex%dm2d)
      call DMClone(dm2d, plex%dm2d, ierr); call CHKERR(ierr)

      allocate(plex%dm)
      call DMClone(dm3d, plex%dm, ierr); call CHKERR(ierr)
      call DMPlexGetChart(plex%dm, pStart, pEnd, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(plex%dm, i3, plex%cStart, plex%cEnd, ierr); call CHKERR(ierr) ! cells
      call DMPlexGetDepthStratum(plex%dm, i2, plex%fStart, plex%fEnd, ierr); call CHKERR(ierr) ! faces / edges
      call DMPlexGetDepthStratum(plex%dm, i1, plex%eStart, plex%eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetDepthStratum(plex%dm, i0, plex%vStart, plex%vEnd, ierr); call CHKERR(ierr) ! vertices

      plex%Nlay = Nlay
      allocate(plex%zindex(pStart:pEnd-1), source=zindex)

      call plex_set_ltopfacepos(plex%dm, plex%ltopfacepos)

      call label_domain_boundary(plex%dm, plex%ltopfacepos, plex%zindex, &
        plex%boundarylabel, plex%domainboundarylabel, plex%ownerlabel)

      call compute_face_geometry(plex, plex%geom_dm)

      call setup_srfc_boundary_dm(plex, plex%srfc_boundary_dm)
      call setup_cell1_dmplex(plex%dm, plex%cell1_dm)

      if(present(hhl)) then
        allocate(plex%hhl1d(size(hhl)), source=hhl)
      endif
    end subroutine

    subroutine gen_test_mat(dm)
      type(tDM), intent(inout) :: dm
      type(tMat) :: A
      integer(mpiint) :: ierr

      call dmplex_set_new_section(dm, 'face_test_section', i1, [i0], [i1], [i0], [i0])

      call DMCreateMatrix(dm, A, ierr); call CHKERR(ierr)
      call MatDestroy(A, ierr); call CHKERR(ierr)
      call CHKERR(1_mpiint, 'DEBUG')
    end subroutine

    subroutine label_domain_boundary(dm, ltopfacepos, zindex, boundarylabel, domainboundarylabel, ownerlabel)
      type(tDM), intent(in) :: dm
      logical, allocatable, intent(in) :: ltopfacepos(:)
      integer(iintegers), allocatable, intent(in) :: zindex(:)
      type(tDMLabel), allocatable, intent(inout) :: boundarylabel, domainboundarylabel, ownerlabel

      type(tDM) :: facedm
      type(tPetscSection) :: facesection
      type(tVec) :: lVec, gVec
      real(ireals), pointer :: xv(:)
      type(tIS) :: boundary_ids

      type(tPetscSF) :: sf
      integer(iintegers) :: nleaves, nroots
      integer(iintegers), pointer :: myidx(:) ! list of my indices that we do not own
      type(PetscSFNode), pointer  :: remote(:) ! rank and remote idx of those points

      integer(iintegers), pointer :: xbndry_iface(:)
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: i, iface, voff, lv
      integer(mpiint) :: comm, myid, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if(.not.allocated(ltopfacepos)) call CHKERR(1_mpiint, 'ltopfacepos has to be allocated')
      if(.not.allocated(zindex)) call CHKERR(1_mpiint, 'zindex has to be allocated')

      if(.not.allocated(ownerlabel)) then
        allocate(ownerlabel)
        call DMCreateLabel(dm, "owner", ierr); call CHKERR(ierr)
      endif
      call DMGetLabel(dm, "owner", ownerlabel, ierr); call CHKERR(ierr)

      if(.not.allocated(boundarylabel)) then
        allocate(boundarylabel)
        call DMCreateLabel(dm, "Boundary", ierr); call CHKERR(ierr)
      endif
      call DMGetLabel(dm, "Boundary", boundarylabel, ierr); call CHKERR(ierr)
      call DMPlexMarkBoundaryFaces(dm, i1, boundarylabel, ierr); call CHKERR(ierr)
      ! The following code should do the same as DMPlexMarkBoundaryFaces...
      ! and is placed here for historic reasons to debug a PETSc issue
      !   call DMLabelSetDefaultValue(boundarylabel, i0, ierr); call CHKERR(ierr)
      !   call DMPlexGetDepthStratum (dm, i2, fStart, fEnd, ierr); call CHKERR(ierr)
      !   do iface = fStart, fEnd-1
      !     call DMPlexGetSupportSize(dm, iface, i, ierr); call CHKERR(ierr)
      !     if(i.lt.i2) then
      !       call DMLabelSetValue(boundarylabel, iface, i1, ierr); call CHKERR(ierr)
      !     endif
      !   enddo

      call DMLabelSetDefaultValue(ownerlabel, int(myid, kind=iintegers), ierr); call CHKERR(ierr)

      call DMGetPointSF(dm, sf, ierr); call CHKERR(ierr)
      call PetscSFGetGraph(sf, nroots, nleaves, myidx, remote, ierr); call CHKERR(ierr)

      do i = 1, nleaves
        call DMLabelSetValue(ownerlabel, myidx(i), remote(i)%rank, ierr); call CHKERR(ierr)
      enddo
      myidx=>NULL(); remote=>NULL()

      call DMClone(dm, facedm, ierr); call CHKERR(ierr)

      call dmplex_set_new_section(facedm, 'Face_Section', i1, [i0], [i1], [i0], [i0])  ! Contains 1 dof on each side
      call DMGetSection(facedm, facesection, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(facedm, gVec, ierr); call CHKERR(ierr)
      call VecSet(gVec, zero, ierr); call CHKERR(ierr)

      call DMGetLocalVector(facedm, lVec, ierr); call CHKERR(ierr)
      call VecSet(lVec, zero, ierr); call CHKERR(ierr)

      call DMGetStratumIS(dm, "Boundary", i1, boundary_ids, ierr); call CHKERR(ierr)
      call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)

      if (boundary_ids.ne.PETSC_NULL_IS) then
        call ISGetIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)
        do i = 1, size(xbndry_iface)
          iface = xbndry_iface(i)

          call DMLabelGetValue(boundarylabel, iface, lv, ierr); call CHKERR(ierr)
          call PetscSectionGetOffset(facesection, iface, voff, ierr); call CHKERR(ierr)
          xv(i1+voff) = real(lv, ireals)
        enddo
        call ISRestoreIndicesF90(boundary_ids, xbndry_iface, ierr); call CHKERR(ierr)
      endif

      call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)

      call DMLocalToGlobalBegin(facedm, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (facedm, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(facedm, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (facedm, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum (facedm, i2, fStart, fEnd, ierr); call CHKERR(ierr)

      if(.not.allocated(domainboundarylabel)) then
        allocate(domainboundarylabel)
        call DMCreateLabel(dm, "DomainBoundary", ierr); call CHKERR(ierr)
      endif
      call DMGetLabel(dm, "DomainBoundary", domainboundarylabel, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(lVec, xv, ierr); call CHKERR(ierr)
      do iface = fStart, fEnd-1
        call PetscSectionGetOffset(facesection, iface, voff, ierr); call CHKERR(ierr)
        select case(nint(xv(i1+voff), iintegers))
        case(i1) ! if the additive val is 1 it must be at the domain edge
          if(ltopfacepos(iface)) then
            if(zindex(iface).eq.1) then
              call DMLabelSetValue(domainboundarylabel, iface, TOAFACE, ierr); call CHKERR(ierr)
            else
              call DMLabelSetValue(domainboundarylabel, iface, BOTFACE, ierr); call CHKERR(ierr)
            endif
          else
            call DMLabelSetValue(domainboundarylabel, iface, SIDEFACE, ierr); call CHKERR(ierr)
          endif
        case(i2) ! if the additive val is 2 it must be at the inner domain edge
          call DMLabelSetValue(domainboundarylabel, iface, INNERSIDEFACE, ierr); call CHKERR(ierr)
        end select
      enddo
      call VecRestoreArrayReadF90(lVec, xv, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(facedm, lVec, ierr); call CHKERR(ierr)
      call DMRestoreGlobalVector(facedm, gVec, ierr); call CHKERR(ierr)
      call DMDestroy(facedm, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_Boundary_DM', ierr); call CHKERR(ierr)
    end subroutine

    subroutine facevec2cellvec(faceVec_dm, global_faceVec, vecshow_string)
      type(tDM), intent(in) :: faceVec_dm
      type(tVec),intent(in) :: global_faceVec
      character(len=*), intent(in), optional :: vecshow_string

      type(tVec) :: cellVec, global_cellVec, faceVec
      real(ireals), pointer :: xcellVec(:)
      real(ireals), pointer :: xfaceVec(:)

      type(tDM) :: celldm
      type(tPetscSection) :: cellSection, faceVecSection
      integer(iintegers) :: icell, cStart, cEnd, cell_offset, iface, faceVec_offset
      integer(iintegers),pointer :: faces_of_cell(:)
      integer(iintegers) :: idof, num_dof, max_num_dof

      integer(mpiint) :: myid, comm, ierr
      character(len=default_str_len) :: faceVecname, cellVecname
      logical :: option_is_set

      call PetscObjectGetComm(faceVec_dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      call PetscObjectGetName(global_faceVec, faceVecname, ierr); call CHKERR(ierr)
      cellVecname = 'fV2cV_'//trim(faceVecname)

      call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-show_'//trim(cellVecname), &
        option_is_set, ierr); call CHKERR(ierr)
      if(.not.option_is_set) return

      if(ldebug.and.myid.eq.0) print *,'facevec2cellvec :: starting..'//trim(faceVecname)

      call DMClone(faceVec_dm, celldm, ierr); ; call CHKERR(ierr)

      call DMGetSection(faceVec_dm, faceVecSection, ierr); call CHKERR(ierr)

      call DMGetLocalVector(faceVec_dm, faceVec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(faceVec_dm, global_faceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (faceVec_dm, global_faceVec, INSERT_VALUES, faceVec, ierr); call CHKERR(ierr)

      call DMPlexGetHeightStratum(celldm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

      ! count number of dof on faces
      max_num_dof = 0
      call DMPlexGetHeightStratum(celldm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      do icell = cStart, cEnd-1
        num_dof = 0
        call DMPlexGetCone(celldm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
        do iface = 1, size(faces_of_cell)
          call PetscSectionGetDof(faceVecSection, faces_of_cell(iface), idof, ierr); call CHKERR(ierr)
          num_dof = num_dof + idof
        enddo
        call DMPlexRestoreCone(celldm, icell, faces_of_cell,ierr); call CHKERR(ierr)
        max_num_dof = max(max_num_dof, num_dof)
      enddo

      ! global max of dof is required to easily plot cell values in visit via xdmf if we do not want to write a custom xmf description
      num_dof = max_num_dof
      call imp_allreduce_max(comm, num_dof, max_num_dof)

      if(ldebug.and.myid.eq.0) then
        print *,'max number of dof on faces per cell is:', max_num_dof
      endif

      call dmplex_set_new_section(celldm, 'Faces_to_Cells_Section', i1, [max_num_dof], [i0], [i0], [i0])

      call DMGetSection(celldm, cellSection, ierr); call CHKERR(ierr)
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
          call PetscSectionGetDof(faceVecSection, faces_of_cell(iface), num_dof, ierr); call CHKERR(ierr)
          do idof = 1, num_dof
            cell_offset = cell_offset+1
            xcellVec(cell_offset) = xfaceVec(faceVec_offset+idof)
          enddo
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

      call PetscObjectSetName(global_cellVec, cellVecname, ierr); call CHKERR(ierr)

      call PetscObjectGetName(global_cellVec, cellVecname, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(global_cellVec, PETSC_NULL_VEC, '-show_'//trim(cellVecname), ierr); call CHKERR(ierr)
      if(present(vecshow_string)) then
        call PetscObjectViewFromOptions(global_cellVec, PETSC_NULL_VEC, trim(vecshow_string), ierr); call CHKERR(ierr)
      endif

      call DMRestoreGlobalVector(celldm, global_cellVec, ierr); call CHKERR(ierr)
      call DMDestroy(celldm, ierr); call CHKERR(ierr)
      if(ldebug.and.myid.eq.0) print *,'facevec2cellvec :: end'//trim(faceVecname)
    end subroutine

    subroutine compute_face_geometry(plex, dm)
      type(t_plexgrid), intent(inout) :: plex
      type(tDM), intent(inout), allocatable :: dm

      type(tVec) :: coordinates
      integer(iintegers) :: cStart, cEnd, icell
      integer(iintegers) :: fStart, fEnd, iface
      integer(iintegers) :: eStart, eEnd, iedge
      integer(iintegers) :: vStart, vEnd, ivert

      integer(iintegers), pointer :: transclosure(:)
      real(ireals), pointer :: coords(:) ! pointer to coordinates vec

      integer(iintegers), target :: vertices6(6)
      integer(iintegers), pointer :: vertices(:)

      type(tPetscSection) :: coordSection, geomSection
      integer(iintegers) :: Ndim, voff0, voff1, voff2, voff3

      real(ireals) :: vertex_coord(3,6) ! shape (Ndim, Nvertices)
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec

      real(ireals) :: area_top, area_bot, volume
      real(ireals) :: centroid_top_face(3), centroid_bot_face(3), height
      integer(iintegers) :: iface_up, iface_dn

      integer(mpiint) :: ierr

      if(allocated(dm)) call CHKERR(1_mpiint, 'called compute_face_geometry on an already allocated geom_dm')
      allocate(dm)
      call DMClone(plex%dm, dm, ierr); ; call CHKERR(ierr)
      call DMGetCoordinateDim(dm, Ndim, ierr); call CHKERR(ierr)
      call DMGetCoordinateSection(dm, coordSection, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(coordSection, PETSC_NULL_SECTION, "-show_dm_coord_section", ierr); call CHKERR(ierr)
      ! Geometry Vec Contains 4 Fields:
      ! field 0: 3 dof for centroid on cells and faces
      ! field 1: 3 for normal vecs on faces
      ! field 2: 1 dof on cells, faces and edges for volume, area, length
      ! field 3: dz on cells
      call dmplex_set_new_section(dm, 'Geometry Section', i4, &
        [i3, i0, i1, i1], &
        [i3, i3, i1, i0], &
        [i0, i0, i1, i0], &
        [i0, i0, i0, i0] )
      call DMGetSection(dm, geomSection, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(geomSection, PETSC_NULL_SECTION, "-show_dm_geom_section", ierr); call CHKERR(ierr)

      call DMGetCoordinatesLocal(dm, coordinates, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_dm_coord", ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(dm, i3, cStart, cEnd, ierr); call CHKERR(ierr)  ! cells
      call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
      call DMPlexGetDepthStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! vertices
      call DMPlexGetDepthStratum(dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      if(.not.allocated(plex%geomVec)) then
        allocate(plex%geomVec)
        call DMCreateLocalVector(dm, plex%geomVec,ierr); call CHKERR(ierr)
        call PetscObjectSetName(plex%geomVec, 'geomVec', ierr);call CHKERR(ierr)
      endif
      call VecGetArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      ! First Set lengths of edges
      do iedge = eStart, eEnd-1
        call DMPlexGetCone(dm, iedge, vertices, ierr); call CHKERR(ierr)

        ! Get the coordinates of vertices
        do ivert=1,size(vertices)
          call PetscSectionGetOffset(coordSection, vertices(ivert), voff0, ierr); call CHKERR(ierr)
          vertex_coord(:, ivert) = coords(i1+voff0:voff0+Ndim)
        enddo

        !print *,'edge length:',iedge,'::', distance(vertex_coord(:,1), vertex_coord(:,2))
        call PetscSectionGetFieldOffset(geomSection, iedge, i2, voff2, ierr); call CHKERR(ierr)
        geoms(i1+voff2) = distance(vertex_coord(:,1), vertex_coord(:,2))

        call DMPlexRestoreCone(dm, iedge, vertices, ierr); call CHKERR(ierr)
      enddo

      ! Then define geom info for faces
      do iface = fStart, fEnd-1
        call PetscSectionGetFieldOffset(geomSection, iface, i0, voff0, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(geomSection, iface, i1, voff1, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(geomSection, iface, i2, voff2, ierr); call CHKERR(ierr)
        call compute_face_geometry_info(dm, iface, &
          centroid=geoms(i1+voff0:voff0+Ndim), &
          normal=geoms(i1+voff1:voff1+Ndim), &
          area=geoms(i1+voff2))
        !print *,'iface', iface, 'center', geoms(i1+voff0:voff0+Ndim)
        !print *,'iface', iface, 'normal', geoms(i1+voff1:voff1+Ndim)
        !print *,'iface', iface, 'area  ', geoms(i1+voff2)
      enddo

      ! Last but not least define geom info for cells
      do icell = cStart, cEnd-1
        call DMPlexGetTransitiveClosure(dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
        !print *,'cell transclosure:',icell, transclosure
        if(size(transclosure).ne.21*2) then
          print *,'len of transclosure with size', size(transclosure), 'not supported -- is this a wedge?'
          call CHKERR(1_mpiint, 'geometry not supported -- is this a wedge?')
        endif

        vertices6 = transclosure(size(transclosure)-11:size(transclosure):2)

        do ivert=1,size(vertices6)
          call PetscSectionGetOffset(coordSection, vertices6(ivert), voff0, ierr); call CHKERR(ierr)
          vertex_coord(1:Ndim, ivert) = coords(i1+voff0:Ndim+voff0)
        enddo

        !print *,'centroid of cell:',icell,'::', sum(vertex_coord(1:Ndim,:),dim=2)/size(vertices6)
        call PetscSectionGetOffset(geomSection, icell, voff0, ierr); call CHKERR(ierr)
        geoms(i1+voff0:voff0+Ndim) = sum(vertex_coord(1:Ndim,:),dim=2)/size(vertices6)

        ! Compute volume of wedges
        iface_up = transclosure(3)
        iface_dn = transclosure(5)

        call PetscSectionGetFieldOffset(geomSection, iface_up, i2, voff2, ierr); call CHKERR(ierr)
        area_top = geoms(i1+voff2)
        call PetscSectionGetFieldOffset(geomSection, iface_dn, i2, voff2, ierr); call CHKERR(ierr)
        area_bot = geoms(i1+voff2)

        call PetscSectionGetFieldOffset(geomSection, iface_up, i0, voff0, ierr); call CHKERR(ierr)
        centroid_top_face = geoms(i1+voff0:voff0+Ndim)
        call PetscSectionGetFieldOffset(geomSection, iface_dn, i0, voff0, ierr); call CHKERR(ierr)
        centroid_bot_face = geoms(i1+voff0:voff0+Ndim)
        height = distance(centroid_top_face, centroid_bot_face)

        !print *,'cell top area', area_top, 'down_area', area_bot, 'height', height

        ! Volume of a "Pyramidenstumpf" -- de.wikipedia.org/wiki/Pyramidenstumpf
        volume = height/3._ireals * (area_top + sqrt(area_top*area_bot) + area_bot)

        call PetscSectionGetFieldOffset(geomSection, icell, i2, voff2, ierr); call CHKERR(ierr)
        geoms(i1+voff2) = volume
        !print *,'cell volume', geoms(i1+voff2)

        call DMPlexRestoreTransitiveClosure(dm, icell, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)

        ! dz as the distance between top and bot face centroids
        call PetscSectionGetFieldOffset(geomSection, icell, i3, voff3, ierr); call CHKERR(ierr)
        geoms(i1+voff3) = height
      enddo

      call VecRestoreArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(plex%geomVec, PETSC_NULL_VEC, "-show_dm_geom_vec", ierr); call CHKERR(ierr)
    end subroutine

    subroutine dmplex_set_new_section(dm, sectionname, numfields, cdof, fdof, edof, vdof, fieldnames)
      type(tDM), intent(inout) :: dm
      character(len=*), intent(in) :: sectionname
      integer(iintegers), intent(in) :: numfields
      integer(iintegers), intent(in) :: cdof(:), fdof(:), edof(:), vdof(:) ! dim=numfields
      character(len=*), intent(in), optional :: fieldnames(:)
      type(tPetscSection) :: section
      integer(mpiint) :: ierr
      integer(iintegers) :: ifield
      logical :: luseCone, luseClosure

      call DMGetBasicAdjacency(dm, luseCone, luseClosure, ierr); call CHKERR(ierr)
      call create_plex_section(dm, sectionname, numfields, cdof, fdof, edof, vdof, section, fieldnames)

      call DMSetLocalSection(dm, section, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(section, ierr); call CHKERR(ierr)

      do ifield = 0, numfields-1
        call DMSetAdjacency(dm, ifield, luseCone, luseClosure, ierr); call CHKERR(ierr)
      enddo
    end subroutine

    subroutine create_plex_section(dm, sectionname, numfields, cdof, fdof, edof, vdof, section, fieldnames)
      type(tDM), intent(in) :: dm
      character(len=*) :: sectionname
      integer(iintegers), intent(in) :: numfields
      integer(iintegers), intent(in) :: cdof(:), fdof(:), edof(:), vdof(:) ! dim=numfields
      type(tPetscSection), intent(inout) :: section
      character(len=*), intent(in), optional :: fieldnames(:)

      integer(iintegers) :: i, depth, ifield, section_size, sum_cdof, sum_fdof, sum_edof, sum_vdof

      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      integer(mpiint) :: comm, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)

      ! This is code to manually create a section, e.g. as in DMPlexCreateSection
      sum_cdof = sum(cdof)
      sum_fdof = sum(fdof)
      sum_edof = sum(edof)
      sum_vdof = sum(vdof)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepth(dm, depth, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
      call DMPlexGetDepthStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetDepthStratum(dm, i3, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

      pEnd = pStart
      if(sum_cdof.gt.i0) pEnd = max(pEnd, cEnd)
      if(sum_fdof.gt.i0) pEnd = max(pEnd, fEnd)
      if(sum_edof.gt.i0) pEnd = max(pEnd, eEnd)
      if(sum_vdof.gt.i0) pEnd = max(pEnd, vEnd)

      pStart = pEnd
      if(sum_vdof.gt.i0) pStart = min(pStart, vStart)
      if(sum_edof.gt.i0) pStart = min(pStart, eStart)
      if(sum_fdof.gt.i0) pStart = min(pStart, fStart)
      if(sum_cdof.gt.i0) pStart = min(pStart, cStart)

      if(depth.lt.i1) &
        call CHKERR(int(sum_edof, mpiint), 'DMPlex does not have edges in the stratum.. cant create section')

      if(depth.lt.i2) &
        call CHKERR(int(sum_fdof, mpiint), 'DMPlex does not have faces in the stratum.. cant create section')

      if(depth.lt.i3) &
        call CHKERR(int(sum_cdof, mpiint), 'DMPlex does not have cells in the stratum.. cant create section')

      ! Create Default Section
      call PetscSectionCreate(comm, section, ierr); call CHKERR(ierr)
      call PetscSectionSetNumFields(section, numfields, ierr); call CHKERR(ierr)
      call PetscSectionSetChart(section, pStart, pEnd, ierr); call CHKERR(ierr)

      if(sum_cdof.gt.i0) then
        do i = cStart, cEnd-i1
          call PetscSectionSetDof(section, i, sum_cdof, ierr); call CHKERR(ierr)
          do ifield = 1, numfields
            call PetscSectionSetFieldDof(section, i, ifield-1, cdof(ifield), ierr); call CHKERR(ierr)
          enddo
        enddo
      endif
      if(sum_fdof.gt.i0) then
        do i = fStart, fEnd-i1
          call PetscSectionSetDof(section, i, sum_fdof, ierr); call CHKERR(ierr)
          do ifield = 1, numfields
            call PetscSectionSetFieldDof(section, i, ifield-1, fdof(ifield), ierr); call CHKERR(ierr)
          enddo
        enddo
      endif
      if(sum_edof.gt.i0) then
        do i = eStart, eEnd-i1
          call PetscSectionSetDof(section, i, sum_edof, ierr); call CHKERR(ierr)
          do ifield = 1, numfields
            call PetscSectionSetFieldDof(section, i, ifield-1, edof(ifield), ierr); call CHKERR(ierr)
          enddo
        enddo
      endif
      if(sum_vdof.gt.i0) then
        do i = vStart, vEnd-i1
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

    subroutine setup_cell1_dmplex(orig_dm, dm)
      type(tDM), intent(in) :: orig_dm
      type(tDM), allocatable, intent(inout) :: dm
      integer(mpiint) :: ierr

      if(allocated(dm)) call CHKERR(1_mpiint, 'called setup_cell1_dmplex on an already allocated DM')
      allocate(dm)

      call DMClone(orig_dm, dm, ierr); call CHKERR(ierr)
      call PetscObjectSetName(dm, 'plex_cell1_dm', ierr);call CHKERR(ierr)
      call dmplex_set_new_section(dm, 'cell_section', i1, [i1], [i0], [i0], [i0])
    end subroutine

    subroutine setup_edir_dmplex(plex, orig_dm, top_streams, side_streams, dof_per_stream, dm)
      type(t_plexgrid), intent(in) :: plex
      type(tDM), intent(in) :: orig_dm
      integer(iintegers), intent(in) :: top_streams, side_streams, dof_per_stream
      type(tDM), allocatable, intent(inout) :: dm
      type(tPetscSection) :: section
      integer(mpiint) :: ierr
      !logical :: luseCone, luseClosure


      if(allocated(dm)) call CHKERR(1_mpiint, 'called setup_edir_dmplex on an already allocated DM')
      allocate(dm)

      call DMClone(orig_dm, dm, ierr); call CHKERR(ierr)

      call PetscObjectSetName(dm, 'plex_direct_radiation', ierr);call CHKERR(ierr)
      call DMSetOptionsPrefix(dm, 'dir', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      call gen_face_section(dm, top_streams=top_streams, side_streams=side_streams, dof_per_stream=dof_per_stream, &
        section=section, geomdm=plex%geom_dm, geomVec=plex%geomVec, aspect_constraint=twostr_ratio)

      call DMSetLocalSection(dm, section, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, '-show_section', ierr); call CHKERR(ierr)
      call PetscSectionDestroy(section, ierr); call CHKERR(ierr)

      call DMSetFromOptions(dm, ierr); call CHKERR(ierr)
    end subroutine

    !> @brief setup the section on which diffuse radiation lives, e.g. 2 dof on top/bot faces and 4 dof on side faces
    !> @details the section has 2 fields, one for incoming radiation and one for outgoing.
    !    The direction for the fields is given as:
    !      * field 0: radiation travelling vertical
    !      * field 1: radiation travelling horizontally
    !
    !      * component 0: radiation travelling from cell_a to cell_b
    !      * component 1: radiation going from cell_b to cell_a
    !      where id of cell_b is larger id(cell_a)
    subroutine setup_ediff_dmplex(plex, orig_dm, top_streams, side_streams, dof_per_stream, dm)
      type(t_plexgrid), intent(in) :: plex
      type(tDM), intent(in) :: orig_dm
      integer(iintegers), intent(in) :: top_streams, side_streams, dof_per_stream
      type(tDM), allocatable, intent(inout) :: dm
      type(tPetscSection) :: section
      integer(mpiint) :: ierr

      if(allocated(dm)) call CHKERR(1_mpiint, 'called setup_ediff_dmplex on an already allocated DM')
      allocate(dm)

      call DMClone(orig_dm, dm, ierr); call CHKERR(ierr)

      call PetscObjectSetName(dm, 'plex_ediff', ierr);call CHKERR(ierr)
      call DMSetOptionsPrefix(dm, 'diff', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      call gen_face_section(dm, top_streams=top_streams, side_streams=side_streams, dof_per_stream=dof_per_stream, &
        section=section, geomdm=plex%geom_dm, geomVec=plex%geomVec, aspect_constraint=twostr_ratio)

      call DMSetLocalSection(dm, section, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, '-show_diff_section', ierr); call CHKERR(ierr)
      call PetscSectionDestroy(section, ierr); call CHKERR(ierr)
      call DMSetFromOptions(dm, ierr); call CHKERR(ierr)
    end subroutine

    subroutine gen_face_section(dm, top_streams, side_streams, dof_per_stream, section, &
        geomdm, geomVec, aspect_constraint)
      type(tDM), intent(in) :: dm
      integer(iintegers), intent(in) :: top_streams, side_streams, dof_per_stream ! number of streams
      type(tPetscSection), intent(inout) :: section
      type(tDM), intent(in), optional :: geomdm
      type(tVec), intent(in), optional :: geomVec
      real(ireals), intent(in), optional :: aspect_constraint

      logical :: l1d
      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:)
      real(ireals) :: face_area, dz, aspect
      integer(iintegers) :: ic, geom_offset, num_constrained, num_unconstrained
      integer(iintegers), pointer :: cells_of_face(:), faces_of_cell(:)

      integer(iintegers) :: iface, fStart, fEnd, istream, num_edges_of_face, tot_top_dof, tot_side_dof
      integer(mpiint) :: comm, ierr

      if(    any([present(geomdm), present(geomVec), present(aspect_constraint)]) .and. &
        .not.all([present(geomdm), present(geomVec), present(aspect_constraint)])) then
        call CHKERR(1_mpiint, 'if you want to constraint dofs with the aspect ratio '// &
          'you have to provide both, the geom dm and a value for the constraint')
      endif

      if(present(aspect_constraint)) then
        call DMGetSection(geomdm, geomSection, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(geomVec, geoms, ierr); call CHKERR(ierr)
      endif

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call PetscSectionCreate(comm, section, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call PetscSectionSetNumFields(section, top_streams+side_streams, ierr); call CHKERR(ierr)

      do istream = 1, top_streams
        call PetscSectionSetFieldName(section, istream-i1, 'Etop'//itoa(istream), ierr); call CHKERR(ierr)
      enddo
      do istream = 1, side_streams
        call PetscSectionSetFieldName(section, top_streams+istream-i1, 'Eside'//itoa(istream), ierr); call CHKERR(ierr)
      enddo

      tot_top_dof = top_streams*dof_per_stream
      tot_side_dof = side_streams*dof_per_stream

      do istream = 1, top_streams+side_streams
        call PetscSectionSetFieldComponents(section, istream-i1, dof_per_stream, ierr); call CHKERR(ierr)
      enddo

      num_constrained = 0
      num_unconstrained = 0
      call PetscSectionSetChart(section, fStart, fEnd, ierr); call CHKERR(ierr)
      do iface = fStart,  fEnd-1
        call DMPlexGetConeSize(dm, iface, num_edges_of_face, ierr); call CHKERR(ierr)

        if(num_edges_of_face.eq.i3) then
          call PetscSectionSetDof(section, iface, tot_top_dof, ierr); call CHKERR(ierr)
          do istream = 1, top_streams
            call PetscSectionSetFieldDof(section, iface, istream-i1, dof_per_stream, ierr); call CHKERR(ierr)
          enddo

        else
          l1d = .False.
          if(present(aspect_constraint)) then
            l1d = .True.

            call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
            face_area = geoms(i1+geom_offset)

            call DMPlexGetSupport(geomdm, iface, cells_of_face, ierr); call CHKERR(ierr)
            do ic = 1, size(cells_of_face)
              call PetscSectionGetFieldOffset(geomSection, cells_of_face(ic), i3, geom_offset, ierr); call CHKERR(ierr)
              dz = geoms(i1+geom_offset)

              call DMPlexGetCone(geomdm, cells_of_face(ic), faces_of_cell, ierr); call CHKERR(ierr)
              call PetscSectionGetFieldOffset(geomSection, faces_of_cell(1), i2, geom_offset, ierr); call CHKERR(ierr)
              face_area = geoms(i1+geom_offset) ! top face area
              call DMPlexRestoreCone(geomdm, cells_of_face(ic), faces_of_cell, ierr); call CHKERR(ierr)

              aspect = LUT_wedge_aspect_zx(face_area, dz)
              if(aspect.lt.aspect_constraint) l1d=.False.
            enddo
            call DMPlexRestoreSupport(geomdm, iface, cells_of_face, ierr); call CHKERR(ierr)
          endif

          if(.not.l1d) then
            call PetscSectionSetDof(section, iface, tot_side_dof, ierr); call CHKERR(ierr)
            do istream = 1, side_streams
              call PetscSectionSetFieldDof(section, iface, top_streams+istream-i1, dof_per_stream, ierr); call CHKERR(ierr)
            enddo
            num_unconstrained = num_unconstrained + 1
          else
            !print *,'Have 1D aspect constraint on face', iface, ':', face_area, dz, dx, '=>', aspect, aspect_constraint, l1d
            num_constrained = num_constrained + 1
          endif
        endif
      enddo
      call PetscSectionSetUp(section, ierr); call CHKERR(ierr)
      if(ldebug) then
        if(num_unconstrained.eq.0) then
          print *,'Have '//itoa(num_constrained)//' constrained dofs : 100%'
        else
          print *,'Have '//itoa(num_constrained)//' constrained dofs :'//&
            ftoa(real(num_constrained, ireals)*100._ireals/real(num_unconstrained, ireals))//' %'
        endif
      endif

      if(present(aspect_constraint)) then
        call VecRestoreArrayReadF90(geomVec, geoms, ierr); call CHKERR(ierr)
      endif
    end subroutine

    subroutine setup_srfc_boundary_dm(plex, dm)
      type(t_plexgrid), intent(inout) :: plex
      type(tDM), allocatable, intent(inout) :: dm
      type(tPetscSection) :: section
      integer(mpiint) :: ierr

      if(allocated(dm)) call CHKERR(1_mpiint, 'called setup_srfc_boundary_dm on an already allocated DM')
      if(.not.allocated(plex%geom_dm)) call CHKERR(1_mpiint, 'plex%geom_dm has to be allocated first')
      allocate(dm)

      call DMClone(plex%dm, dm, ierr); call CHKERR(ierr)

      call PetscObjectSetName(dm, 'srfc_boundary_dm', ierr);call CHKERR(ierr)

      call gen_section(section, Ndof=i1)

      call DMSetSection(dm, section, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(section, ierr); call CHKERR(ierr)
    contains
      subroutine gen_section(section, Ndof)
        type(tPetscSection), intent(inout) :: section
        integer(iintegers), intent(in) :: Ndof
        type(tIS) :: srfc_ids
        integer(iintegers), pointer :: xi(:)
        integer(iintegers) :: i, iface, fStart, fEnd
        integer(mpiint) :: comm, ierr

        call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
        call PetscSectionCreate(comm, section, ierr); call CHKERR(ierr)
        call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
        call PetscSectionSetNumFields(section, i1, ierr); call CHKERR(ierr)
        call PetscSectionSetFieldComponents(section, i0, i1, ierr); call CHKERR(ierr)
        call PetscSectionSetChart(section, fStart, fEnd, ierr); call CHKERR(ierr)

        call DMGetStratumIS(plex%dm, 'DomainBoundary', BOTFACE, srfc_ids, ierr); call CHKERR(ierr)
        if (srfc_ids.eq.PETSC_NULL_IS) then ! dont have surface points
        else
          call ISGetIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
          do i = 1, size(xi)
            iface = xi(i)
            call PetscSectionSetDof(section, iface, Ndof, ierr); call CHKERR(ierr)
            call PetscSectionSetFieldDof(section, iface, i0, Ndof, ierr); call CHKERR(ierr)
          enddo
          call ISRestoreIndicesF90(srfc_ids, xi, ierr); call CHKERR(ierr)
        endif
        call PetscSectionSetUp(section, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(section, PETSC_NULL_SECTION, "-show_plex_srfc_boundary_section", ierr); call CHKERR(ierr)
      end subroutine
    end subroutine

    subroutine determine_diff_incoming_outgoing_offsets(plex, ediffdm, icell, incoming_offsets, outgoing_offsets)
      type(t_plexgrid), intent(in) :: plex
      type(tDM), intent(in) :: ediffdm
      integer(iintegers), intent(in) :: icell
      integer(iintegers), allocatable, intent(inout) :: incoming_offsets(:), outgoing_offsets(:)

      integer(iintegers), pointer :: faces_of_cell(:)
      integer(iintegers), pointer :: cells_of_face(:)
      integer(iintegers) :: i, iface, istream, neigh_cell, offset_a, offset_b
      integer(iintegers) :: j_incoming, j_outgoing, num_dof, num_fields, idof
      integer(iintegers) :: boundarylabelval, owner
      type(tPetscSection) :: section
      integer(mpiint) :: comm, myid, ierr

      call DMGetSection(ediffdm, section, ierr); call CHKERR(ierr)
      call PetscSectionGetNumFields(section, num_fields, ierr); call CHKERR(ierr)

      call DMPlexGetCone(ediffdm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      num_dof = 0
      do i = 1, size(faces_of_cell)
        iface = faces_of_cell(i)
        call PetscSectionGetDof(section, iface, idof, ierr); call CHKERR(ierr)
        num_dof = num_dof + idof
      enddo

      if(allocated(incoming_offsets)) then
        if(size(incoming_offsets).ne.num_dof/2) deallocate(incoming_offsets)
      endif
      if(allocated(outgoing_offsets)) then
        if(size(outgoing_offsets).ne.num_dof/2) deallocate(outgoing_offsets)
      endif

      if(.not.allocated(incoming_offsets)) allocate(incoming_offsets(num_dof/2))
      if(.not.allocated(outgoing_offsets)) allocate(outgoing_offsets(num_dof/2))

      j_incoming = 1; j_outgoing = 1

      do i = 1, size(faces_of_cell)
        iface = faces_of_cell(i)
        call DMPlexGetSupport(ediffdm, iface, cells_of_face, ierr); call CHKERR(ierr) ! Get Faces of cell
        if(size(cells_of_face).eq.1) then
          ! This is either because we are at the outer domain or this is a local mesh boundary with a neighboring process
          call DMLabelGetValue(plex%domainboundarylabel, iface, boundarylabelval, ierr); call CHKERR(ierr)
          if(boundarylabelval.ne.-i1) then ! This is a global boundary face, i.e. at the side top or bottom of the domain
            neigh_cell = -1
          else
            call DMLabelGetValue(plex%ownerlabel, iface, owner, ierr); call CHKERR(ierr)
            call PetscObjectGetComm(ediffdm, comm, ierr); call CHKERR(ierr)
            call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
            if(owner.eq.myid) then
              neigh_cell = -1
            else
              neigh_cell = huge(icell)
            endif
          endif
        else
          if(cells_of_face(1).eq.icell) then
            neigh_cell = cells_of_face(2)
          else
            neigh_cell = cells_of_face(1)
          endif
        endif
        call DMPlexRestoreSupport(ediffdm, iface, cells_of_face, ierr); call CHKERR(ierr) ! Get Faces of cell

        do istream = 1, num_fields
          call PetscSectionGetFieldDof(section, iface, istream-i1, num_dof, ierr); call CHKERR(ierr)
          if(num_dof.eq.2) then
            call PetscSectionGetFieldOffset(section, iface, istream-i1, offset_a, ierr); call CHKERR(ierr)
            offset_b = offset_a+1
            !print *,icell,'iface', iface, 'stream', istream, 'offset', offset_a, offset_b

            if(neigh_cell.lt.icell) then ! see definition of directions in setup_ediff_dmplex
              incoming_offsets(j_incoming) = offset_a
              j_incoming = j_incoming + 1
              outgoing_offsets(j_outgoing) = offset_b
              j_outgoing = j_outgoing + 1
            else
              incoming_offsets(j_incoming) = offset_b
              j_incoming = j_incoming + 1
              outgoing_offsets(j_outgoing) = offset_a
              j_outgoing = j_outgoing + 1
            endif
          endif
        enddo
      enddo

      call DMPlexRestoreCone(ediffdm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
    end subroutine

    subroutine setup_abso_dmplex(orig_dm, dm)
      type(tDM), intent(in) :: orig_dm
      type(tDM), allocatable, intent(inout) :: dm
      integer(mpiint) :: ierr

      if(allocated(dm)) stop 'called setup_abso_dmplex on an already allocated DM'
      allocate(dm)

      call DMClone(orig_dm, dm, ierr); call CHKERR(ierr)

      call PetscObjectSetName(dm, 'plex_abso', ierr);call CHKERR(ierr)
      call DMSetOptionsPrefix(dm, 'abso', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      call dmplex_set_new_section(dm, 'absorption_section', i1, [i1], [i0], [i0], [i0])
      call DMSetFromOptions(dm, ierr); call CHKERR(ierr)
    end subroutine

    !> @brief reorient face normal so that they point towards the same direction as the sun does
    subroutine orient_face_normals_along_sundir(plex, sundir)
      type(t_plexgrid), intent(in) :: plex
      real(ireals), intent(in) :: sundir(3)

      integer(iintegers) :: iface

      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      integer(iintegers) :: fStart, fEnd, geom_offset
      real(ireals) :: face_normal(3), mu

      integer(mpiint) :: ierr

      call DMPlexGetDepthStratum(plex%dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! cells

      call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      do iface = fStart, fEnd-1
        call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)
        face_normal = geoms(geom_offset+i4: geom_offset+i6)

        mu = dot_product(sundir/norm2(sundir), face_normal)
        if(mu.lt.zero) then ! normal is in the opposite direction of the sun -> turn it around
          geoms(geom_offset+i4: geom_offset+i6) = -face_normal
        endif
      enddo

      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
    end subroutine

    !> @brief get face normal which points inwards into a cell
    subroutine get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)
      integer(iintegers), intent(in) :: iface, icell
      type(tPetscSection), intent(in) :: geomSection
      real(ireals), intent(in), pointer :: geoms(:)
      real(ireals), intent(out) :: face_normal(:)

      real(ireals) :: cell_center(3), face_center(3)
      integer(iintegers) :: voff0, voff1, n
      integer(mpiint) :: ierr

      call PetscSectionGetFieldOffset(geomSection, icell, i0, voff0, ierr); call CHKERR(ierr)
      cell_center = geoms(voff0+i1: voff0+i3)

      call PetscSectionGetFieldOffset(geomSection, iface, i0, voff0, ierr); call CHKERR(ierr)
      call PetscSectionGetFieldOffset(geomSection, iface, i1, voff1, ierr); call CHKERR(ierr)
      face_center = geoms(voff0+i1: voff0+i3)
      face_normal = geoms(voff1+i1: voff1+i3) ! DEBUG this is the 3D Normal
      ! face_normal = geoms(voff1+i1+i3: voff1+i3+i3) ! DEBUG use 2d normal instead

      if(ldebug) then
        if(.not.approx(one, norm2(face_normal))) then
          print *,'cell', icell, 'face', iface, 'Face Normal:', face_normal
          call CHKERR(1_mpiint, 'face_normal not normed :( '//ftoa(face_normal))
        endif
      endif

      n = determine_normal_direction(face_normal, face_center, cell_center)

      face_normal = face_normal * real(n, ireals)
    end subroutine

    !> @brief create a vector that holds all the wedge orientation ordering information
    !> as well as zenith and azimuth angles
    !> , i.e. in short, everything that we need to know to translate between DMPlex face ordering and the LUT's
    subroutine compute_wedge_orientation(plex, sundir, wedge_orientation_dm, wedge_orientation)
      type(t_plexgrid), intent(in) :: plex
      real(ireals), intent(in) :: sundir(3)
      type(tDM), allocatable, intent(inout) :: wedge_orientation_dm
      type(tVec), allocatable, intent(inout) :: wedge_orientation

      type(tPetscSection) :: geomSection, wedgeSection

      integer(iintegers) :: icell, iface, wedge_offset1, wedge_offset2, wedge_offset3, wedge_offset4
      integer(iintegers) :: upper_face, base_face, left_face, right_face, bottom_face
      logical :: lsrc(5)

      real(ireals) :: zenith, azimuth, param_phi, param_theta, Cx, Cy, aspect_zx
      real(ireals),pointer :: xv(:), geoms(:)

      integer(mpiint) :: ierr

      if(.not.allocated(wedge_orientation_dm)) then
        allocate(wedge_orientation_dm)
        call DMClone(plex%edir_dm, wedge_orientation_dm, ierr); ; call CHKERR(ierr)
        ! wedge_orientation Vec Contains 4 Fields:
        ! 5 dof on cells for zenith, azimuth, param_phi, param_theta, aspect_zx
        ! 5 dof for permutation of faces from faces_of_cell to BoxMonteCarlo face ordering
        ! 5 dof that determine if a face is src or destination with respect to solar radiation (src=1, dst=0)
        ! 2 dof on cells for C point coordinates in local boxmc geometry(2D coords)
        call dmplex_set_new_section(wedge_orientation_dm, 'Wedge_Orientation_Section', i4, &
          [i5, i5, i5, i2], &
          [i0, i0, i0, i0], &
          [i0, i0, i0, i0], &
          [i0, i0, i0, i0])
      endif
      call DMGetSection(wedge_orientation_dm, wedgeSection, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(wedgeSection, PETSC_NULL_SECTION, '-show_WedgeSection', ierr); call CHKERR(ierr)

      if(.not.allocated(wedge_orientation)) then
        allocate(wedge_orientation)
        call DMCreateGlobalVector(wedge_orientation_dm, wedge_orientation, ierr); call CHKERR(ierr)
        call PetscObjectSetName(wedge_orientation, 'WedgeOrient', ierr);call CHKERR(ierr)
      endif
      call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call VecGetArrayF90(wedge_orientation, xv, ierr); call CHKERR(ierr)

      do icell = plex%cStart, plex%cEnd-1
        call compute_local_wedge_ordering(plex, icell, &
          geomSection, geoms, sundir, &
          lsrc, zenith, azimuth, param_phi, param_theta, &
          Cx, Cy, aspect_zx, &
          upper_face, bottom_face, &
          base_face, left_face, right_face)

        call PetscSectionGetFieldOffset(wedgeSection, icell, i0, wedge_offset1, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(wedgeSection, icell, i1, wedge_offset2, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(wedgeSection, icell, i2, wedge_offset3, ierr); call CHKERR(ierr)
        call PetscSectionGetFieldOffset(wedgeSection, icell, i3, wedge_offset4, ierr); call CHKERR(ierr)

        xv(wedge_offset1+i1: wedge_offset1+i5) = [zenith, azimuth, param_phi, param_theta, aspect_zx]

        ! set iwedge_plex2bmc, index mapping from plex faces_of_cell to bmc
        xv(wedge_offset2+i1: wedge_offset2+i5) = real([upper_face, base_face, left_face, right_face, bottom_face], ireals)
        do iface = 1, size(lsrc)
          if(lsrc(iface)) then
            xv(wedge_offset3+iface) = one
          else
            xv(wedge_offset3+iface) = -one
          endif
        enddo

        xv(wedge_offset4+i1: wedge_offset4+i2) = [Cx, Cy]
      enddo

      call VecRestoreArrayF90(wedge_orientation, xv, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(wedge_orientation, PETSC_NULL_VEC, '-show_WedgeOrient', ierr); call CHKERR(ierr)
    end subroutine

    !> @brief translate the ordering of faces in the DMPlex to the ordering we assume in the box-montecarlo routines
    subroutine compute_local_wedge_ordering(plex, icell, geomSection, geoms, in_sundir, &
        lsrc, zenith, azimuth, param_phi, param_theta, Cx, Cy, aspect_zx, &
        upper_face, bottom_face, base_face, left_face, right_face)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: icell
      type(tPetscSection) :: geomSection
      real(ireals), intent(in), pointer :: geoms(:) ! pointer to coordinates vec
      real(ireals), intent(in) :: in_sundir(3)
      real(ireals), intent(out) :: zenith, azimuth, param_phi, param_theta
      real(ireals), intent(out) :: Cx, Cy, aspect_zx
      integer(iintegers), intent(out) :: upper_face, bottom_face, base_face, left_face, right_face
      logical, intent(out) :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

      integer(iintegers), pointer :: faces_of_cell(:)!, edges_of_face(:)
      integer(iintegers) :: iface, izindex(2) !,i

      real(ireals) :: face_normals(3,5)
      real(ireals) :: sundir(3), proj_angles_to_sun(3), proj_sundir(3)

      real(ireals) :: e_x(3) ! unit vectors of local coord system in which we compute the transfer coefficients
      real(ireals) :: MrotWorld2Local(3,3) ! Rotation Matrix into the local wedge space, (ex, ey, ez)
      ! Normal vectors in local wedgemc geometry.
      ! For even sided triangles, this is smth. like left: [.5, -.866] or right [-.5, -.866]
      real(ireals) :: local_normal_base(3), local_normal_left(3), local_normal_right(3)

      integer(iintegers) :: side_faces(3), top_faces(2) ! indices in faces_of_cell which give the top/bot and side faces via labeling
      integer(iintegers) :: iside_faces, itop_faces, ibase_face, iright_face ! indices to fill above arrays

      real(ireals) :: side_face_normal_projected_on_upperface(3,3)

      real(ireal_params) :: rparam_phi, rparam_theta, n2(3), n3(3), n4(3)

      integer(iintegers) :: geom_offset, iedge
      integer(iintegers), target :: points(2)
      integer(iintegers), pointer :: ppoints(:), coveredPoints(:)
      real(ireals) :: dz, Atop

      integer(mpiint) :: ierr

      sundir = in_sundir / norm2(in_sundir)

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      ! get numbering for 2 top/bot faces and 3 side faces which indicate the position in the faces_of_cell vec
      iside_faces = 1; itop_faces = 1
      do iface = 1, size(faces_of_cell)
        if(plex%ltopfacepos(faces_of_cell(iface))) then
          top_faces(itop_faces) = iface
          itop_faces = itop_faces+1
        else
          side_faces(iside_faces) = iface
          iside_faces = iside_faces+1
        endif
      enddo

      izindex(1) = plex%zindex(faces_of_cell(top_faces(1)))
      izindex(2) = plex%zindex(faces_of_cell(top_faces(2)))

      if(izindex(1).gt.izindex(2)) then
        upper_face = top_faces(2)
        bottom_face = top_faces(1)
      else
        upper_face = top_faces(1)
        bottom_face = top_faces(2)
      endif

      ! Now we have to find the solar azimuth, zenith angles
      do iface = 1, size(faces_of_cell)
        call get_inward_face_normal(faces_of_cell(iface), icell, geomSection, geoms, face_normals(:, iface))
        lsrc(iface) = is_solar_src(face_normals(:,iface), sundir)
      enddo

      zenith = angle_between_two_vec(sundir, face_normals(:, upper_face))
      proj_sundir = vec_proj_on_plane(sundir, face_normals(:,upper_face))
      call normalize_vec(proj_sundir, ierr) !; call CHKERR(ierr, 'bad proj_sundir'//ftoa(proj_sundir))

      do iface=1,size(side_faces)
        side_face_normal_projected_on_upperface(:, iface) = &
          vec_proj_on_plane(face_normals(:,side_faces(iface)), face_normals(:,upper_face))
        call normalize_vec(side_face_normal_projected_on_upperface(:, iface), &
                           side_face_normal_projected_on_upperface(:, iface), ierr)
        call CHKERR(ierr, 'bad side face normal '//ftoa(side_face_normal_projected_on_upperface(:, iface)))
        !if(norm2(proj_sundir).lt.epsilon(zero)) then
        !  proj_angles_to_sun(iface) = zero
        !  lsrc(side_faces(iface)) = .False.
        !else

        proj_angles_to_sun(iface) = &
          angle_between_two_normed_vec(proj_sundir, side_face_normal_projected_on_upperface(:, iface))

        !endif
        !print *,'iface',iface, ':', lsrc(side_faces(iface)), '->', proj_angles_to_sun(iface), rad2deg(proj_angles_to_sun(iface))
      enddo

      if(norm2(proj_sundir).gt.10*epsilon(zenith)) then ! only do the azimuth computation if zenith is larger than 0 deg
        ibase_face = minloc(proj_angles_to_sun,dim=1)
        base_face  = side_faces(ibase_face)

        ! Local unit vec on upperface, pointing towards '+x'
        e_x = cross_3d(face_normals(:, upper_face), side_face_normal_projected_on_upperface(:, ibase_face))

        azimuth = proj_angles_to_sun(ibase_face)
        azimuth = azimuth * sign(one, dot_product(proj_sundir, e_x))

        left_face  = side_faces(modulo(ibase_face,size(side_faces, kind=iintegers))+i1)
        iright_face = modulo(ibase_face+i1,size(side_faces, kind=iintegers))+i1
        right_face = side_faces(iright_face)

      else ! if sun is directly on top, i.e. zenith 0, just pick the first one
        ibase_face = 1
        do iface = 2, size(side_faces)
          if(lsrc(side_faces(iface))) ibase_face = iface
        enddo
        e_x = cross_3d(face_normals(:, upper_face), side_face_normal_projected_on_upperface(:, ibase_face))
        base_face  = side_faces(ibase_face)
        left_face  = side_faces(modulo(ibase_face,size(side_faces, kind=iintegers))+i1)
        iright_face = modulo(ibase_face+i1,size(side_faces, kind=iintegers))+i1
        right_face = side_faces(iright_face)

        !print *,'side_face_normal_projected_on_upperface ibase_face', side_face_normal_projected_on_upperface(:, ibase_face)
        !print *,'face normal upper face', -face_normals(:, upper_face)
        !print *,'e_x', e_x
        !print *,'side_face_normal_projected_on_upperface', side_face_normal_projected_on_upperface(:, iright_face)

        azimuth = zero
        zenith = zero
      endif

      if(dot_product(e_x, side_face_normal_projected_on_upperface(:, iright_face)).gt.zero) then
        call swap(right_face, left_face)
      endif

      !print *,'norm proj_sundir', norm2(proj_sundir),':', ibase_face, rad2deg(azimuth), ':', rad2deg(zenith),&
      !  '::',norm2(e_x), norm2(side_face_normal_projected_on_upperface(:, iright_face))

      ! Determine edge length of edge between base face and upper face triangle
      ppoints => points
      points = [faces_of_cell(upper_face), faces_of_cell(base_face)]
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      iedge = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      ! Determine mean edge length of upper face triangle
      !call DMPlexGetCone(plex%edir_dm, faces_of_cell(upper_face), edges_of_face, ierr); call CHKERR(ierr)
      !dx = 0
      !do i=1,size(edges_of_face)
      !  iedge = edges_of_face(i)
      !  call PetscSectionGetOffset(geomSection, iedge, geom_offset, ierr); call CHKERR(ierr)
      !  dx = dx + geoms(geom_offset+i1)
      !  !print *,'Upper Edgelength', faces_of_cell(upper_face), ':', i, geoms(geom_offset+i1)
      !enddo
      !dx = dx / size(edges_of_face)
      !call DMPlexRestoreCone(plex%edir_dm, faces_of_cell(upper_face), edges_of_face, ierr); call CHKERR(ierr)

      !call PetscSectionGetOffset(geomSection, iedge, geom_offset, ierr); call CHKERR(ierr)
      !dx = geoms(geom_offset+i1)
      !print *,'base face edgelength', iedge, dx

      call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr)
      dz = geoms(i1+geom_offset)
      call PetscSectionGetFieldOffset(geomSection, faces_of_cell(upper_face), i2, geom_offset, ierr); call CHKERR(ierr)
      Atop = geoms(i1+geom_offset)
      aspect_zx = LUT_wedge_aspect_zx(Atop, dz)

      ! Compute local vertex coordinates for upper face
      call compute_local_vertex_coordinates(plex, &
        faces_of_cell(upper_face), &
        faces_of_cell(base_face), &
        faces_of_cell(left_face), &
        faces_of_cell(right_face), &
        Cx, Cy)

      MrotWorld2Local = rotation_matrix_world_to_local_basis(&
        e_x, &
        side_face_normal_projected_on_upperface(:, ibase_face), &
        -face_normals(:, upper_face))

      local_normal_base = matmul(MrotWorld2Local, face_normals(:, base_face))
      local_normal_left = matmul(MrotWorld2Local, face_normals(:, left_face))
      local_normal_right= matmul(MrotWorld2Local, face_normals(:, right_face))

      ! renormalize because of precision gain/loss from ireals to ireal_params
      call normalize_vec(real(local_normal_base , ireal_params), n2, ierr)
      call normalize_vec(real(local_normal_left , ireal_params), n3, ierr)
      call normalize_vec(real(local_normal_right, ireal_params), n4, ierr)

      call param_phi_param_theta_from_phi_and_theta_withnormals(&
        n2, n3, n4, &
        real(Cx, ireal_params), real(Cy, ireal_params), &
        real(azimuth, ireal_params), real(zenith, ireal_params), &
        rparam_phi, rparam_theta, ierr); call CHKERR(ierr)

      !print *,'normals base', local_normal_base
      !print *,'normals left', local_normal_left
      !print *,'normals righte', local_normal_right
      !print *,'Cx/Cy', Cx, Cy
      !print *,'rparam_phi, rparam_theta',rparam_phi, rparam_theta

      if(is_between(rparam_theta, real(param_eps, ireal_params), 1._ireal_params)) then ! only if baseface should be src
        ierr = 0
        if(rparam_phi.lt.-1_irealLUT-param_eps .and. .not.lsrc(left_face)) then
          ierr = 1
          print *,'param_phi < -1 but left face is not src', rparam_phi, rparam_theta, base_face, left_face, right_face, lsrc
        endif

        if(rparam_phi.gt.-1_irealLUT+param_eps .and.      lsrc(left_face)) then
          ierr = 2
          print *,'param_phi > -1 but left face is not dst', rparam_phi, rparam_theta, base_face, left_face, right_face, lsrc
        endif

        if(rparam_phi.gt.1_irealLUT+param_eps .and. .not.lsrc(right_face)) then
          ierr = 3
          print *,'param_phi > 1 but right face is not src', rparam_phi, rparam_theta, base_face, left_face, right_face, lsrc
        endif
        if(rparam_phi.lt.1_irealLUT-param_eps .and. lsrc(right_face)) then
          ierr = 4
          print *,'param_phi < 1 but right face is not dst', rparam_phi, rparam_theta, base_face, left_face, right_face, lsrc
        endif
        if(ierr.ne.0) then
          print *,'angles between local coords', &
            dot_product(e_x, side_face_normal_projected_on_upperface(:, ibase_face)), &
            dot_product(e_x, -face_normals(:, upper_face))
          print *,'proj_angles_to_sun', rad2deg(proj_angles_to_sun)
          print *,'azimuth, zenith ', azimuth, zenith, rad2deg(azimuth), rad2deg(zenith), &
            'param phi/theta', rparam_phi, rparam_theta
          print *,'Cx, Cy', Cx, Cy
          print *,'local sundir      ', matmul(MrotWorld2Local, sundir)
          print *,'local normal top  ', matmul(MrotWorld2Local, face_normals(:, upper_face))
          print *,'local normal base ', local_normal_base, norm2(local_normal_base),  ':', &
            dot_product(local_normal_base, matmul(MrotWorld2Local, sundir))
          print *,'local normal left ', local_normal_left, norm2(local_normal_left),  ':', &
            dot_product(local_normal_left, matmul(MrotWorld2Local, sundir))
          print *,'local normal right', local_normal_right,norm2(local_normal_right), ':', &
            dot_product(local_normal_right, matmul(MrotWorld2Local, sundir))

          do iface = 1, size(faces_of_cell)
            print *,'iface', iface, 'solar_src', is_solar_src(face_normals(:,iface), sundir), &
              ':', sundir/norm2(sundir), ',', &
              face_normals(:, iface)/norm2(face_normals(:, iface)), &
              '=', dot_product(sundir/norm2(sundir), face_normals(:, iface)/norm2(face_normals(:, iface)))
          enddo
          print *, 'local azimuth ', rad2deg(angle_between_two_vec(&
            vec_proj_on_plane(local_normal_base,               matmul(MrotWorld2Local, face_normals(:, upper_face))), &
            vec_proj_on_plane(matmul(MrotWorld2Local, sundir), matmul(MrotWorld2Local, face_normals(:, upper_face))) &
            ))

        endif
        call CHKERR(ierr, 'found bad param_phi')
      endif

      ! Snap param_phi to the correct side
      if(approx(rparam_phi, -1._ireal_params, real(param_eps, kind(rparam_phi)))) then
        if(lsrc(left_face)) then
          rparam_phi = -1._ireal_params-param_eps
        else
          rparam_phi = -1._ireal_params+param_eps
        endif
      elseif(approx(rparam_phi, +1._ireal_params, real(param_eps, kind(rparam_phi)))) then
        if(lsrc(right_face)) then
          rparam_phi = 1._ireal_params+param_eps
        else
          rparam_phi = 1._ireal_params-param_eps
        endif
      endif

      param_phi   = real(rparam_phi, ireals)
      param_theta = real(rparam_theta, ireals)

      if(ldebug) then
        if(param_theta.gt.zero .and. zenith.le.pi/2 .and. .not.lsrc(base_face)) then
          print *,'azimuth', rad2deg(azimuth), 'zenith', rad2deg(zenith)
          print *,'param_phi', param_phi, 'param_theta', param_theta
          print *,'ibase_face', ibase_face, 'baseface', base_face, 'lsrc', lsrc
          print *,'proj_normal', proj_sundir, '::norm', norm2(proj_sundir)
          print *,'face_normals(:,base_face)', face_normals(:,base_face)
          print *,'face_normals(:,left_face)', face_normals(:,left_face)
          print *,'face_normals(:,right_face)', face_normals(:,right_face)
          call CHKERR(1_mpiint, 'base face is not a src! this should not be the case!')
        endif
        if(rad2deg(azimuth).lt.-90 .or. rad2deg(azimuth).gt.90) then
          print *,'ibase_face', ibase_face
          print *,'proj_normal', proj_sundir, '::norm', norm2(proj_sundir)
          print *,'face_normals(:,base_face)', face_normals(:,base_face)
          print *,'face_normals(:,left_face)', face_normals(:,left_face)
          print *,'face_normals(:,right_face)', face_normals(:,right_face)

          ierr = int(rad2deg(azimuth), mpiint)
          call CHKERR(ierr, 'local azimuth greater than 90 deg. something must have gone wrong with the base face selection!')
        endif
      endif


      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)
    end subroutine

    subroutine compute_local_vertex_coordinates(plex, upperface, baseface, leftface, rightface, Cx, Cy)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: upperface, baseface, leftface, rightface
      real(ireals), intent(out) :: Cx, Cy  ! vertex coordinates of C point in local boxmc geometry

      type(tPetscSection) :: coordSection
      type(tVec) :: coordinates
      real(ireals), pointer :: coords(:) ! pointer to coordinates vec

      integer(iintegers), target :: points(2)
      integer(iintegers), pointer :: ppoints(:), coveredPoints(:)
      integer(iintegers) :: iedges(3), ivertices(3)
      integer(iintegers) :: i, voff
      real(ireals) :: vertex_coord(3,3), AB(3), BC(3), AC(3), a, b, c, mu

      integer(mpiint) :: ierr

      ppoints => points
      points = [upperface, baseface]
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      iedges(1) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      points = [upperface, leftface]
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      iedges(2) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      points = [upperface, rightface]
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      iedges(3) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      points = iedges(1:2)
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      ivertices(1) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      points = iedges([1,3])
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      ivertices(2) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      points = iedges(2:3)
      call DMPlexGetMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)
      ivertices(3) = coveredPoints(1)
      call DMPlexRestoreMeet(plex%edir_dm, i2, ppoints, coveredPoints, ierr); call CHKERR(ierr)

      call DMGetCoordinateSection(plex%edir_dm, coordSection, ierr); call CHKERR(ierr)
      call DMGetCoordinatesLocal(plex%edir_dm, coordinates, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)
      do i=1,3
        call PetscSectionGetOffset(coordSection, ivertices(i), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, i) = coords(voff+i1:voff+i3)
      enddo
      call VecRestoreArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)

      AB = vertex_coord(:, 2) - vertex_coord(:, 1)
      BC = vertex_coord(:, 3) - vertex_coord(:, 2)
      AC = vertex_coord(:, 3) - vertex_coord(:, 1)

      c = norm2(AB)
      a = norm2(BC)
      b = norm2(AC)

      !print *,'A', ivertices(1), vertex_coord(:, 1),':AB', AB,' ::', a
      !print *,'B', ivertices(2), vertex_coord(:, 2),':BC', BC,' ::', b
      !print *,'C', ivertices(3), vertex_coord(:, 3),':AC', AC,' ::', c

      !nAB = cross_3d(vertex_coord(:, 1), AB)
      !nAB = nAB / norm2(nAB)
      !print *,'nAB', nAB

      ! law of cosines to get angle between AB and CA
      mu = (b**2 + c**2 - a**2) / (2*b*c)

      Cx = mu * b
      Cy = sqrt(one - mu**2) * b

      Cx = Cx / c
      Cy = Cy / c

      !print *,'local_coords C', Cx, Cy
    end subroutine

    function is_solar_src(face_normal, sundir)
      real(ireals),intent(in) :: face_normal(3), sundir(3)
      logical :: is_solar_src
      real(ireals) :: mu

      mu = dot_product(sundir/norm2(sundir), face_normal/norm2(face_normal))
      is_solar_src = mu.gt.zero
      !print *,'is_solar_src', face_normal, sundir, ':', mu, '->', is_solar_src,'::', rad2deg(acos(mu))
    end function

    subroutine print_dmplex(comm, dm)
      integer(mpiint), intent(in) :: comm
      type(tDM),intent(in) :: dm

      integer(mpiint) :: i, myid, numnodes, ierr
      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(dm, i3, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call DMPlexGetDepthStratum(dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
      call DMPlexGetDepthStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetDepthStratum(dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      do i = 0, numnodes-1
      if(myid.eq.i) then
        print *,myid,'pStart,End :: ',pStart, pEnd
        print *,myid,'cStart,End :: ',cStart, cEnd
        print *,myid,'fStart,End :: ',fStart, fEnd
        print *,myid,'eStart,End :: ',eStart, eEnd
        print *,myid,'vStart,End :: ',vStart, vEnd
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
      enddo
    end subroutine

    subroutine ncvar2d_to_globalvec(plexgrid, filename, varname, gvec, timeidx, cell_ao_2d, cell_ao_3d)
      type(t_plexgrid), intent(in) :: plexgrid
      character(len=*), intent(in) :: filename, varname
      type(tVec), allocatable, intent(inout) :: gvec
      integer(iintegers), intent(in), optional :: timeidx
      AO, optional, intent(in), target :: cell_ao_2d, cell_ao_3d ! mapping into 2D or 3D plex cells on rank0

      type(tDM) :: celldm
      type(tVec) :: rank0Vec
      type(tVecScatter) :: scatter_context
      real(ireals), pointer :: xloc(:)=>null()
      integer(mpiint) :: myid, ierr
      integer(iintegers) :: ic, icell, icell_global, k
      real(ireals), allocatable :: arr(:,:)
      real(ireals), allocatable :: arr3d(:,:,:)
      character(len=default_str_len) :: ncgroups(2)

      if(.not.(present(cell_ao_2d) .or. present(cell_ao_3d))) call CHKERR(1_mpiint, 'have to provide one of cell AO`s')
      if(present(cell_ao_2d) .and. present(cell_ao_3d)) call CHKERR(1_mpiint, 'please provide only one cell AO')

      call mpi_comm_rank(plexgrid%comm, myid, ierr); call CHKERR(ierr)
      if(ldebug.and.myid.eq.0) print *,'Loading from nc file: ',trim(filename), ' varname: ', trim(varname)

      call print_dmplex(plexgrid%comm, plexgrid%dm)

      call DMClone(plexgrid%dm, celldm, ierr); call CHKERR(ierr)
      call dmplex_set_new_section(celldm, 'Cell_Section', i1, [i1], [i0], [i0], [i0])

      ! Now lets get vectors!
      if(.not.allocated(gvec)) then
        allocate(gvec)
        call DMCreateGlobalVector(celldm, gvec, ierr); call CHKERR(ierr)
        call PetscObjectSetName(gvec, 'VecfromNC_'//trim(varname), ierr);call CHKERR(ierr)
      endif

      call VecScatterCreateToZero(gvec, scatter_context, rank0Vec, ierr); call CHKERR(ierr)
      if(ldebug.and.present(cell_ao_2d)) then
        call AOView(cell_ao_2d, PETSC_VIEWER_STDOUT_WORLD, ierr); call CHKERR(ierr)
      endif

      if(myid.eq.0) then
        ncgroups(1) = trim(filename)
        ncgroups(2) = trim(varname)
        if(ldebug) print *,'plex_grid::ncvar2d_to_globalvec : Loading file ', trim(ncgroups(1)), ':', trim(ncgroups(2))
        call ncload(ncgroups, arr, ierr)
        if(ierr.ne.0) then !try to load 3D data
          if(ldebug) print *,'Could not load 2D array Data, trying 3D var'
          call ncload(ncgroups, arr3d, ierr); call CHKERR(ierr, 'Could not load Data from NetCDF')

          k = get_arg(i1, timeidx)
          if(k.lt.i1 .or. k.gt.ubound(arr3d,3)) &
            call CHKERR(1_mpiint, 'invalid time index, shape of arr:' &
            //itoa(size(arr3d,1))//itoa(size(arr3d,2))//itoa(size(arr3d,3)))

          allocate(arr(size(arr3d,1), size(arr3d,2)), source=arr3d(:,:,k))
          deallocate(arr3d)
        endif

        if(ldebug) print *,'plex_grid::ncvar2d_to_globalvec : shape ',trim(varname), shape(arr)

        call VecGetArrayF90(rank0Vec,xloc,ierr); call CHKERR(ierr)
        do icell=1, size(arr, dim=1)
          if(present(cell_ao_2d)) then
            icell_global = icell-1
            ! icell_global will be mapped to petsc index assuming the AO maps cells in C ordering from 0 to Ncells-1
            call AOApplicationToPetsc(cell_ao_2d, i1, icell_global, ierr); call CHKERR(ierr)
            do k=1,size(arr, dim=2)
              ic = icell_global*size(arr, dim=2) + k
              xloc(ic) = arr(icell,k)
            enddo
          else
            do k=1,size(arr, dim=2)
              icell_global = icell*size(arr, dim=2) + k -i1
              call AOApplicationToPetsc(cell_ao_3d, i1, icell_global, ierr); call CHKERR(ierr)
              xloc(icell_global) = arr(icell,k)
            enddo
          endif
        enddo
        call VecRestoreArrayF90(rank0Vec,xloc,ierr) ;call CHKERR(ierr)

      endif ! rank 0

      if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: scatter reverse....'
      call VecScatterBegin(scatter_context, rank0Vec, gvec, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
      call VecScatterEnd  (scatter_context, rank0Vec, gvec, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

      call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)
      call VecDestroy(rank0Vec, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(gvec, PETSC_NULL_VEC, '-show_ncvar_to_globalvec', ierr); call CHKERR(ierr)
    end subroutine

    !> @brief return the dmplex cell index for an icon base grid cell index
    function icell_icon_2_plex(plex, icell, k)
      type(t_plexgrid), intent(in) :: plex      !< @param[in] dmplex mesh object, holding info about number of grid cells
      integer(iintegers),intent(in) :: icell    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
      integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
      integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay) then
          print *,'icell_icon_2_plex :: inp',icell, k
          stop 'icell_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1) then
          print *,'icell_icon_2_plex :: inp',icell, k
          stop 'icell_icon_2_plex :: icon cell index out of range'
        endif
      endif
      icell_icon_2_plex = (k-i1) + (icell-i1)*plex%Nlay
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
        offset = icon%parNfaces(owner) * plex%Nlay
      else
        offset = plex%offset_faces
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay+1) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner)
          stop 'iface_top_icon_2_plex :: vertical index k out of range'
        endif
        if(icell.lt.i1) then
          print *,'iface_top_icon_2_plex :: inp', icell, k, present(owner)
          stop 'iface_top_icon_2_plex :: icon cell index out of range'
        endif
      endif
      iface_top_icon_2_plex = offset + (k-i1) + (icell-i1)*(plex%Nlay+1)
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
        offset = Nfaces * plex%Nlay + (plex%Nlay + 1) * Nfaces
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        offset = plex%offset_faces_sides
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay) then
          print *,'iface_side_icon_2_plex :: inp', iedge, k, present(owner), '::', Nfaces, Nedges, offset
          stop 'iface_side_icon_2_plex :: vertical index k out of range'
        endif
        if(iedge.lt.i1 .or. iedge.gt.Nedges) then
          print *,'iface_side_icon_2_plex :: inp', iedge, k, present(owner), '::', Nfaces, Nedges, offset
          stop 'iface_side_icon_2_plex :: icon edge index out of range'
        endif
      endif
      iface_side_icon_2_plex = offset + (k-i1) + (iedge-i1)*plex%Nlay
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
        offset = Nfaces * plex%Nlay + Nfaces * (plex%Nlay+1) + Nedges * plex%Nlay
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        offset = plex%offset_edges
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay+1) stop 'iedge_top_icon_2_plex :: vertical index k out of range'
        if(iedge.lt.i1 .or. iedge.gt.Nedges) stop 'iedge_top_icon_2_plex :: icon cell index out of range'
      endif

      iedge_top_icon_2_plex = offset + (k-i1) + (iedge-i1)*(plex%Nlay+i1)
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
        offset = Nfaces * plex%Nlay + Nfaces * (plex%Nlay+1) + Nedges * plex%Nlay + Nedges * (plex%Nlay+i1)
      else
        Nfaces = icon%Nfaces
        Nedges = icon%Nedges
        Nvertices = icon%Nvertices
        offset = plex%offset_edges_vertical
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay) stop 'iedge_side_icon_2_plex :: vertical index k out of range'
        if(ivertex.lt.i1 .or. ivertex.gt.Nvertices) stop 'iedge_side_icon_2_plex :: icon vertex index out of range'
      endif

      iedge_side_icon_2_plex = offset + (k-i1) + (ivertex-i1)*plex%Nlay
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
        offset = Nfaces * plex%Nlay + Nfaces * (plex%Nlay+1) + Nedges * plex%Nlay + Nedges * (plex%Nlay+1) + Nvertices * plex%Nlay
      else
        Nfaces    = icon%Nfaces
        Nedges    = icon%Nedges
        Nvertices = icon%Nvertices
        offset    = plex%offset_vertices
      endif
      if(ldebug) then
        if(k.lt.i1 .or. k.gt.plex%Nlay+1) then
          print *,'ivertex_icon_2_plex error! input was',ivertex, k
          stop 'ivertex_side_icon_2_plex :: vertical index k out of range'
        endif
        if(ivertex.lt.i1 .or. ivertex.gt.Nvertices) stop 'ivertex_side_icon_2_plex :: icon vertex index out of range'
      endif

      ivertex_icon_2_plex = offset + (k-i1) + (ivertex-i1)*(plex%Nlay+i1)
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

    subroutine get_consecutive_vertical_cell_idx(plex, startcell, idx)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in) :: startcell
      integer(iintegers), allocatable, intent(inout) :: idx(:)
      integer(iintegers), allocatable :: dbg_idx(:)
      integer(iintegers) :: k

      if(.not.allocated(idx)) then
        allocate(idx(plex%Nlay))
      else
        call CHKERR(int(size(idx)-plex%Nlay, mpiint), 'wrong size of cell idx')
      endif

      idx(1) = startcell
      do k=2,plex%Nlay
        idx(k) = idx(k-1) + 1
      enddo

      if(ldebug) then
        ! make sure that startcell is at TOA
        call CHKERR(int(plex%zindex(startcell)-1, mpiint), 'startcell has to be at TOA to use this routine')

        call get_vertical_cell_idx(plex%dm, startcell, plex%Nlay, dbg_idx)
        if(.not.all(idx.eq.dbg_idx)) then
          print *,'consecutive idx', idx
          print *,'plex_scan   idx', dbg_idx
          call CHKERR(1_mpiint, 'consecutive and dmplex scan vertical cell idx set does not give the same results')
        endif
      endif
    end subroutine

    subroutine get_vertical_cell_idx(dm, startcell, idx_maxsize, idx)
      type(tDM), intent(in) :: dm
      integer(iintegers), intent(in) :: startcell, idx_maxsize
      integer(iintegers), allocatable, intent(out) :: idx(:)

      integer(iintegers) :: icell, iface, num_edges_of_face, next_cell, neigh_cells(2)
      integer(iintegers), pointer :: faces_of_cell(:), cell_support(:)

      integer(iintegers) :: i, j, k

      integer(mpiint) :: ierr

      if(ldebug) then
        if(idx_maxsize.lt.1) call CHKERR(1_mpiint, 'idx_maxsize should be a positive int, dont you want to get a result?')
        icell = startcell
        call DMPlexGetCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
        call DMPlexGetConeSize(dm, faces_of_cell(1), num_edges_of_face, ierr); call CHKERR(ierr)
        call CHKERR(int(num_edges_of_face-i3, mpiint), 'we assume that first entry of faces is top face, should have 3 edges')
        call DMPlexGetConeSize(dm, faces_of_cell(2), num_edges_of_face, ierr); call CHKERR(ierr)
        call CHKERR(int(num_edges_of_face-i3, mpiint), 'we assume that secondentry of faces is bot face, should have 3 edges')

        ! we assume that we start at one end of the columns
        ! i.e. one of the cell supports has to be the startcell
        iface = faces_of_cell(1)
        call DMPlexGetSupport(dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        call CHKERR(int(minval(abs(cell_support-startcell)), mpiint), &
          'startcell has to be a cell at either end of the column')
        call DMPlexRestoreSupport(dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        call DMPlexRestoreCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      endif

      allocate(idx(idx_maxsize), source=-i1)
      idx(1) = startcell
      k = 2

      next_cell = startcell
      do
        neigh_cells = neighboring_cells_of_cell(dm, next_cell)
        j = k
        do i = 1, size(neigh_cells)
          if(neigh_cells(i).lt.i0) cycle
          if(any(neigh_cells(i).eq.idx(max(i1,[k-i1, k-i2])))) cycle ! already have this neigh cell in index set
          if(ldebug) then
            if(k.gt.idx_maxsize) &
              call CHKERR(int(k, mpiint), 'found more cells than we have space to store in idx('//itoa(idx_maxsize)//')')
            !print *,'Adding new cell ', neigh_cells(i), '=>', idx
          endif
          idx(k) = neigh_cells(i)
          k = k+1
          next_cell = neigh_cells(i)
        enddo
        if(j.eq.k) then ! didnt add any new cells, exiting loop
          exit
        endif
      enddo
      call resize_arr(k-1, idx)
      !call PetscSortInt(k-1, idx, ierr); call CHKERR(ierr)
      !print *,'cell idx in column of cell',startcell,':',idx
    end subroutine

    function neighboring_cells_of_cell(dm, icell)
      type(tDM), intent(in) :: dm
      integer(iintegers), intent(in) :: icell
      integer(iintegers) :: neighboring_cells_of_cell(2)
      integer(iintegers), pointer :: faces_of_cell(:), cell_support(:)
      integer(iintegers) :: i, j, k, num_edges_of_face
      integer(mpiint) :: ierr

      neighboring_cells_of_cell(:) = -1
      k = 1
      call DMPlexGetCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      !print *,icell, 'Faces of cell:', faces_of_cell
      do i = 1, size(faces_of_cell)
        call DMPlexGetConeSize(dm, faces_of_cell(i), num_edges_of_face, ierr); call CHKERR(ierr)
        !print *,faces_of_cell(i), 'Number of edges', num_edges_of_face
        if(num_edges_of_face.eq.3) then ! horizontal face
          call DMPlexGetSupport(dm, faces_of_cell(i), cell_support, ierr); call CHKERR(ierr)
          !print *,faces_of_cell(i), 'cell_support', cell_support
          do j = 1, size(cell_support)
            if(cell_support(j).ne.icell) then
              neighboring_cells_of_cell(k) = cell_support(j)
              k = k+1
            endif
          enddo
          call DMPlexRestoreSupport(dm, faces_of_cell(i), cell_support, ierr); call CHKERR(ierr)
        endif
      enddo
      call DMPlexRestoreCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      !print *,'Neighboring cells of ', icell, ':', neighboring_cells_of_cell
    end function

    subroutine get_top_bot_face_of_cell(dm, icell, topface, botface)
      type(tDM), intent(in) :: dm
      integer(iintegers), intent(in) :: icell
      integer(iintegers), intent(out) :: topface, botface
      integer(iintegers), pointer :: faces_of_cell(:)
      real(ireals) :: centroid_top(3), centroid_bot(3), normal(3), area
      integer(mpiint) :: ierr

      call DMPlexGetCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      topface = min(faces_of_cell(1),faces_of_cell(2))
      botface = max(faces_of_cell(1),faces_of_cell(2))
      call DMPlexRestoreCone(dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

      if(ldebug) then
        call compute_face_geometry_info(dm, topface, centroid_top, normal, area)
        call compute_face_geometry_info(dm, botface, centroid_bot, normal, area)
        if(norm2(centroid_top).lt.norm2(centroid_bot)) then
          print *,'Centroid top/bot', centroid_top, ':', centroid_bot, '=>', norm2(centroid_top), norm2(centroid_bot)
          call CHKERR(1_mpiint, 'Hmpf, I guessed wrong, I confused top and bot face of cell')
        endif
      endif
    end subroutine

    subroutine compute_face_geometry_info(dm, iface, centroid, normal, area)
      type(tDM), intent(in) :: dm
      integer(iintegers), intent(in) :: iface
      real(ireals), intent(out) :: centroid(3), normal(3), area

      integer(iintegers), pointer :: transclosure(:)
      real(ireals), pointer :: coords(:) ! pointer to coordinates vec

      integer(iintegers), target :: vertices3(3), vertices4(4)
      integer(iintegers), pointer :: vertices(:)
      real(ireals), allocatable :: vertex_coord(:,:) ! shape (Nvertices, dims)
      real(ireals) :: AB(3), CD(3), dotp!, normals(3,2)

      type(tVec) :: coordinates
      type(tPetscSection) :: coordSection
      integer(iintegers) :: ivert, Ndim, Nvertices, voff, Nedges
      integer(mpiint) :: ierr

      call DMGetCoordinateDim(dm, Ndim, ierr); call CHKERR(ierr)
      call DMGetCoordinateSection(dm, coordSection, ierr); call CHKERR(ierr)
      call DMPlexGetConeSize(dm, iface, Nedges, ierr); call CHKERR(ierr)
      call DMGetCoordinatesLocal(dm, coordinates, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)

      vertices => NULL()
      Nvertices = 0
      select case(Nedges)
      case (i3)
        vertices => vertices3
        Nvertices = i3

      case (i4)
        vertices => vertices4
        Nvertices = i4

      case default
        vertices => vertices3 ! to get rid of -Werror=maybe-uninitialized
        print *,'Computation of face normal not implemented for faces with ',Nedges,'edges/vertices'
        call CHKERR(1_mpiint, 'Invalid number of edges for face normal computation')
      end select

      call DMPlexGetTransitiveClosure(dm, iface, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)
      !print *,'transclosure', iface,'::',Nedges,'::',transclosure

      vertices = transclosure(3+2*Nedges:size(transclosure):2) ! indices come from: 1 for faceid, Nedge, and then vertices, each with two entries, one for index, and one for orientation
      !print *,'iface',iface,'vertices', vertices

      call DMPlexRestoreTransitiveClosure(dm, iface, PETSC_TRUE, transclosure, ierr); call CHKERR(ierr)

      ! Get the coordinates of vertices
      allocate(vertex_coord(Ndim, Nvertices))
      do ivert=1,size(vertices)
        call PetscSectionGetOffset(coordSection, vertices(ivert), voff, ierr); call CHKERR(ierr)
        vertex_coord(:, ivert) = coords(i1+voff:voff+Ndim)
        !print *,'iface',iface,'vertex',vertices(ivert),'::',coords(1+voff:voff+Ndim)
      enddo

      !print *,'centroid of face:',iface,'::', sum(vertex_coord,dim=2)/real(Nvertices,ireals)
      centroid = sum(vertex_coord,dim=2) / real(Nvertices, kind=ireals)

      ! and use 3 coordinates to compute normal
      !normals(:,1) = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,2),vertex_coord(:,3))
      !if(Nvertices.gt.3) then
      !  normals(:,2) = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,3),vertex_coord(:,Nvertices))
      !  if(dot_product(normals(:,1),normals(:,2)).le.zero) normals(:,2) = -normals(:,2)
      !  normals(:,1) = (normals(:,1) + normals(:,2))/2
      !  normals(:,1) = normals(:,1) / norm2(normals(:,1))
      !  !print *,'normal of face', iface,'::', normals(:,2)
      !endif
      !
      !normals(:,2) = normals(:,2) * sign(one, normals(:,1))
      !print *,'normal of face', iface,'::', normals(:,1)
      !print *,'-------------------------------'

      normal = compute_normal_3d(vertex_coord(:,1),vertex_coord(:,2),vertex_coord(:,3))
      if(.not.approx(norm2(normal), one, 10*epsilon(one))) &
        call CHKERR(1_mpiint, 'face normal not normed :( '//ftoa(normal)//' ( '//ftoa(norm2(normal))//' )')

      if(Nvertices.eq.3) then
        area = triangle_area_by_vertices(vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3))
        !print *,'face triangle area:', geoms(i1+voff+Ndim*2)
      else
        AB = vertex_coord(:,2) - vertex_coord(:,1)
        CD = vertex_coord(:,4) - vertex_coord(:,3)
        dotp = dot_product(AB/norm2(AB),CD/norm2(CD))
        if(dotp.gt.zero) then
          !print *,'swapping vertex coordinates because dot_product>0', dotp
          call swap(vertex_coord(:,1),vertex_coord(:,2))
        endif
        area = triangle_area_by_vertices(vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3)) + &
          triangle_area_by_vertices(vertex_coord(:,3), vertex_coord(:,4), vertex_coord(:,1))
        !print *,'face rectangle area:', geoms(i1+voff+Ndim*2)
      endif

      deallocate(vertex_coord)
      call VecRestoreArrayReadF90(coordinates, coords, ierr); call CHKERR(ierr)
    end subroutine

    function get_normal_of_first_TOA_face(plex, facenr)
      type(t_plexgrid), intent(in) :: plex
      integer(iintegers), intent(in), optional :: facenr
      real(ireals) :: get_normal_of_first_TOA_face(3)
      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      real(ireals) :: cell_center(3), face_center(3)
      real(ireals),allocatable :: face_normal(:)

      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xitoa(:), cell_support(:)
      integer(iintegers) :: geom_offset, iface, icell

      integer(mpiint) :: myid, comm, ierr

      call PetscObjectGetComm(plex%dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        if(.not.allocated(plex%geom_dm)) stop 'get_normal_of_first_TOA_face::needs allocated geom_dm first'
        call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
        call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)

        if (toa_ids.eq.PETSC_NULL_IS) then ! dont have TOA points
          call CHKERR(1_mpiint, 'This didnt work, we tried to set the sundir according to first face on rank 0'// &
                                'but it seems he does not have TOA faces')
        else
          call ISGetIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
          if(present(facenr)) then
            if(facenr.lt.1.or.facenr.gt.size(xitoa)) call CHKERR(1_mpiint, 'bad facenr: has to be in range '//itoa(shape(xitoa)))
            iface = xitoa(facenr)
          else
            iface = xitoa(1) ! first face of TOA faces
          endif

          call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
          call PetscSectionGetFieldOffset(geomSection, iface, i0, geom_offset, ierr); call CHKERR(ierr)

          face_center = geoms(i1+geom_offset:geom_offset+i3)
          allocate(face_normal(3))
          call PetscSectionGetFieldOffset(geomSection, iface, i1, geom_offset, ierr); call CHKERR(ierr)
          face_normal = geoms(i1+geom_offset:geom_offset+3)

          call DMPlexGetSupport(plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell

          call PetscSectionGetFieldOffset(geomSection, icell, i0, geom_offset, ierr); call CHKERR(ierr)
          cell_center = geoms(i1+geom_offset:geom_offset+i3)

          ! Determine the inward normal vec for the face
          face_normal = face_normal * real(determine_normal_direction(face_normal, face_center, cell_center), kind=ireals)

          call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

          call ISRestoreIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
        endif
      endif

      call imp_bcast(plex%comm, face_normal, 0_mpiint)
      get_normal_of_first_TOA_face = face_normal

    end function

    subroutine atm_dz_to_vertex_heights(atm_dz, dm3d, TOA_height)
      real(ireals), intent(in) :: atm_dz(:,:)              ! shape(Nlay, Ncol) values start at the surface
      type(tDM), intent(inout) :: dm3d
      real(ireals), intent(in), optional :: TOA_height

      type(tDM) :: facedm, vertdm
      type(tPetscSection) :: face_section, vert_section, coord_section
      type(tVec) :: level_heights_vec, coordinates
      type(tVec), allocatable :: vertvec

      real(ireals), pointer :: xlvl_hgt(:), xvert(:), coords(:)

      integer(iintegers) :: Nlay, Ncol, Nvert, face_section_size
      integer(iintegers) :: vStart, vEnd, fStart, fEnd
      integer(iintegers) :: k, icol, iface, srfc_vert, toa_face
      integer(iintegers) :: voff, coff

      real(ireals) :: srfc_distance, max_TOA_height

      integer(mpiint) :: ierr

      max_TOA_height = get_arg(120e3_ireals, TOA_height) ! height at which the mesh is homogenized, from there counting down

      call DMClone(dm3d, facedm, ierr); call CHKERR(ierr)
      call gen_face_section(facedm, i1, i0, i1, face_section)

      call DMSetSection(facedm, face_section, ierr); call CHKERR(ierr)

      call DMClone(dm3d, vertdm, ierr); call CHKERR(ierr)
      call dmplex_set_new_section(vertdm, 'vert_section', i1, [i0], [i0], [i0], [i1])
      call DMGetLocalSection(vertdm, vert_section, ierr); call CHKERR(ierr)

      Nlay = size(atm_dz, 1)
      Ncol = size(atm_dz, 2)
      call DMPlexGetDepthStratum(facedm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! vertices
      call PetscSectionGetStorageSize(face_section, face_section_size, ierr); call CHKERR(ierr)

      call CHKERR(int(face_section_size-Ncol*(Nlay+1), mpiint), &
        'bad face section size or bad atm_dz size '//itoa(face_section_size)//' vs '//itoa(Ncol*(Nlay+1)))

      call DMPlexGetDepthStratum(vertdm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
      Nvert = vEnd - vStart

      call DMGetLocalVector(facedm, level_heights_vec, ierr); call CHKERR(ierr)
      call VecGetArrayF90    (level_heights_vec, xlvl_hgt, ierr); call CHKERR(ierr)

      do icol = 1, Ncol
        toa_face = fStart + (icol-i1)*(Nlay+1)

        call PetscSectionGetOffset(face_section, toa_face, voff, ierr); call CHKERR(ierr)
        xlvl_hgt(i1+voff) = max_TOA_height

        do k = 1, Nlay
          iface = toa_face + k
          xlvl_hgt(i1+voff+k) = xlvl_hgt(i1+voff+k-1) - atm_dz(Nlay-(k-1),icol)
        enddo
        !print *,'srfc_face', srfc_face, icol,':', xlvl_hgt(i1+voff:i1+voff-Nlay:-1)
      enddo
      ! lvl_hgt has now the height starting from srfc height and add up the atm_dz's

      call VecRestoreArrayF90(level_heights_vec, xlvl_hgt, ierr); call CHKERR(ierr)

      if(ldebug) then
        call PetscObjectSetName(level_heights_vec, 'level_heights_vec', ierr); call CHKERR(ierr)
        call facevec2cellvec(facedm, level_heights_vec)
      endif

      call interpolate_horizontal_face_var_onto_vertices(facedm, level_heights_vec, vertdm, vertvec)

      call DMRestoreLocalVector(facedm, level_heights_vec, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(face_section, ierr); call CHKERR(ierr)
      call DMDestroy(facedm, ierr); call CHKERR(ierr)

      ! Now generate new coordinates for 3D DMPlex
      call DMGetCoordinateSection(dm3d, coord_section, ierr); call CHKERR(ierr)
      call DMGetCoordinatesLocal(dm3d, coordinates, ierr); call CHKERR(ierr)
      call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

      call VecGetArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)
      do srfc_vert = vStart+Nlay, vEnd-1, Nlay+1
        call PetscSectionGetOffset(coord_section, srfc_vert, coff, ierr); call CHKERR(ierr)

        call PetscSectionGetOffset(vert_section, srfc_vert, voff, ierr); call CHKERR(ierr)

        srfc_distance = norm2(coords(i1+coff:coff+i3)) ! distance from origin till vertex at surface

        do k = Nlay, 0, -1
          call PetscSectionGetOffset(coord_section, srfc_vert-k, coff, ierr); call CHKERR(ierr)
          coords(i1+coff:coff+i3) = coords(i1+coff:coff+i3) / norm2(coords(i1+coff:coff+i3)) &
            * (srfc_distance + xvert(i1+voff-k))
          !print *,'xvert', srfc_vert, xvert(i1+voff-k), ':coords', coords(i1+coff:coff+i3),'fac',(srfc_distance + xvert(i1+voff-k))
        enddo
      enddo

      call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(vert_section, ierr); call CHKERR(ierr)
      call DMDestroy(vertdm, ierr); call CHKERR(ierr)
    end subroutine

    ! takes the average of horizontal face values around a vertex and builds the mean
    subroutine interpolate_horizontal_face_var_onto_vertices(facedm, facevec, vertdm, vertvec)
      type(tDM), intent(in) :: facedm, vertdm
      type(tVec), intent(in) :: facevec
      type(tVec), allocatable, intent(inout) :: vertvec ! Local vector on vertices

      type(tVec) :: Numvec, gVec

      type(tPetscSection) :: face_section, vert_section

      integer(iintegers) :: ivert, vStart, vEnd, j, iface, voff, foff
      integer(iintegers), allocatable :: faces_around_vert(:)

      real(ireals), pointer :: xface(:), xvert(:), xNum(:)
      integer(mpiint) :: ierr

      !integer(mpiint) :: comm, myid
      !call PetscObjectGetComm(facedm, comm, ierr); call CHKERR(ierr)
      !call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      call DMGetSection(facedm, face_section, ierr); call CHKERR(ierr)
      call DMGetSection(vertdm, vert_section, ierr); call CHKERR(ierr)

      if(.not.allocated(vertvec)) then
        allocate(vertvec)
        call DMCreateLocalVector(vertdm, vertvec, ierr); call CHKERR(ierr)
      endif
      call VecSet(vertvec, zero, ierr); call CHKERR(ierr)

      call DMGetLocalVector(vertdm, Numvec, ierr); call CHKERR(ierr)
      call VecSet(Numvec, zero, ierr); call CHKERR(ierr)

      call DMPlexGetDepthStratum(vertdm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      call VecGetArrayReadF90(facevec, xface, ierr); call CHKERR(ierr)
      call VecGetArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)
      call VecGetArrayF90(Numvec, xNum, ierr); call CHKERR(ierr)

      do ivert = vStart, vEnd-1
        call get_horizontal_faces_around_vertex(facedm, ivert, faces_around_vert)

        ! add to Nsumvec, this will be used to take the average in the end
        call PetscSectionGetOffset(vert_section, ivert, voff, ierr); call CHKERR(ierr)
        xNum(i1+voff) = size(faces_around_vert)

        do j=1,size(faces_around_vert)
          iface = faces_around_vert(j)
          call PetscSectionGetOffset(face_section, iface, foff, ierr); call CHKERR(ierr)
          xvert(i1+voff) = xvert(i1+voff) + xface(i1+foff)
          !print *,myid,'ivert', ivert, 'adding', iface, xface(i1+foff), '=', xvert(i1+voff), 'xnum', xNum(i1+voff)
        enddo
      enddo

      call VecRestoreArrayF90(Numvec, xNum, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)

      call VecRestoreArrayReadF90(facevec, xface, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(vertdm, gVec, ierr); call CHKERR(ierr)

      call VecSet(gVec, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(vertdm, Numvec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (vertdm, Numvec, ADD_VALUES, gVec, ierr)
      call DMGlobalToLocalBegin(vertdm, gVec, INSERT_VALUES, Numvec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (vertdm, gVec, INSERT_VALUES, Numvec, ierr)

      call VecSet(gVec, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(vertdm, vertvec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd  (vertdm, vertvec, ADD_VALUES, gVec, ierr)
      call DMGlobalToLocalBegin(vertdm, gVec, INSERT_VALUES, vertvec, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (vertdm, gVec, INSERT_VALUES, vertvec, ierr)

      call DMRestoreGlobalVector(vertdm, gVec, ierr); call CHKERR(ierr)

      ! call VecPointwiseDivide(vertvec, vertvec, Numvec, ierr); call CHKERR(ierr)
      ! Take the average
      call VecGetArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(Numvec, xNum, ierr); call CHKERR(ierr)
      xvert = xvert / xNum
      call VecRestoreArrayReadF90(Numvec, xNum, ierr); call CHKERR(ierr)
      call VecRestoreArrayF90(vertvec, xvert, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(vertdm, Numvec, ierr); call CHKERR(ierr)
    end subroutine

    subroutine get_horizontal_faces_around_vertex(dm, ivert, idx)
      type(tDM), intent(in) :: dm
      integer(iintegers),intent(in) :: ivert
      integer(iintegers), allocatable, intent(out) :: idx(:)

      type(tDMLabel) :: depthlabel
      integer(iintegers) :: i, j, numfaces, numSupport
      logical :: lcontains
      integer(iintegers), pointer :: transclosure(:)
      integer(mpiint) :: ierr

      call DMPlexGetDepthLabel(dm, depthLabel, ierr); call CHKERR(ierr)

      call DMPlexGetTransitiveClosure(dm, ivert, PETSC_FALSE, transclosure, ierr); call CHKERR(ierr)
      numfaces = 0
      do i=1,size(transclosure),2
        call DMLabelStratumHasPoint(depthlabel, i2, transclosure(i), lcontains, ierr); call CHKERR(ierr)
        if(lcontains) then
          call DMPlexGetConeSize(dm, transclosure(i), numSupport, ierr); call CHKERR(ierr)
          if(numSupport.eq.i3) then ! if face has 3 edges
            numfaces = numfaces + 1
          else
            transclosure(i) = -100 -transclosure(i)
          endif
        else
          transclosure(i) = -100 -transclosure(i)
        endif
      enddo
      allocate(idx(numfaces))
      j=1
      do i=1,size(transclosure),2
        if(transclosure(i).ge.0) then
          idx(j) = transclosure(i)
          j = j+1
        endif
      enddo
      call DMPlexRestoreTransitiveClosure(dm, ivert, PETSC_FALSE, transclosure, ierr); call CHKERR(ierr)
    end subroutine
  end module
