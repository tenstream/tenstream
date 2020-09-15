module m_pprts2plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    & one, i0

  use m_helper_functions, only: CHKWARN, CHKERR, approx, toStr, &
    & is_inrange, &
    & ndarray_offsets, ind_1d_to_nd, ind_nd_to_1d

  use m_buildings, only: t_plex_buildings, t_pprts_buildings, &
    & PPRTS_TOP_FACE,                     &
    & PPRTS_BOT_FACE,                     &
    & PPRTS_LEFT_FACE,                    &
    & PPRTS_RIGHT_FACE,                   &
    & PPRTS_REAR_FACE,                    &
    & PPRTS_FRONT_FACE

  use m_pprts_base, only : t_coord

  use m_plex_grid, only: t_plexgrid, &
    & get_inward_face_normal

  implicit none

contains

  subroutine pprts_buildings_to_plex(plex, pprts_buildings, plex_buildings, ierr)
    type(t_plexgrid), intent(in) :: plex
    type(t_pprts_buildings), intent(in) :: pprts_buildings
    type(t_plex_buildings), intent(inout), allocatable :: plex_buildings
    integer(mpiint), intent(out) :: ierr
    integer(iintegers) :: i, k, icell, iface
    integer(iintegers) :: nr_plex_faces, ndidx(4), cell_offsets(3)
    integer(iintegers) :: cStart, cEnd, fStart, fEnd

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:)
    integer(iintegers) :: geom_offset

    ierr = 0

    if(.not.allocated(pprts_buildings%albedo)) ierr = 1
    if(.not.allocated(pprts_buildings%iface)) ierr = 2
    call CHKERR(ierr, 'bad input data in pprts_buildings')

    ! count number of needed new faces
    associate(P => pprts_buildings)
      nr_plex_faces = 0
      do i = lbound(P%iface,1), ubound(P%iface,1)

        call ind_1d_to_nd(P%da_offsets, P%iface(i), ndidx)
        select case(ndidx(1))

        case (PPRTS_TOP_FACE, PPRTS_BOT_FACE)
          nr_plex_faces = nr_plex_faces + 2

        case (PPRTS_LEFT_FACE, PPRTS_RIGHT_FACE, PPRTS_REAR_FACE, PPRTS_FRONT_FACE)
          nr_plex_faces = nr_plex_faces + 1

        case default
          ierr = -1
          call CHKERR(ierr, 'dont know the faceidx')
        end select

      enddo
    end associate

    if(.not.allocated(plex_buildings)) then
      allocate(plex_buildings)
      allocate(plex_buildings%albedo(nr_plex_faces))
      allocate(plex_buildings%iface (nr_plex_faces))
    else
      call CHKERR(int(size(plex_buildings%iface )-2*size(pprts_buildings%iface ),mpiint), "wrong size plex_buildings%iface ")
      call CHKERR(int(size(plex_buildings%albedo)-2*size(pprts_buildings%albedo),mpiint), "wrong size plex_buildings%albedo")
    endif

    call DMPlexGetDepthStratum(plex%dm, 3_iintegers, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetDepthStratum(plex%dm, 2_iintegers, fStart, fEnd, ierr); call CHKERR(ierr) ! faces

    call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    associate(P => pprts_buildings, T=>plex_buildings)
      cell_offsets = P%da_offsets(2:4)/P%da_offsets(2)
      cell_offsets(2:3) = cell_offsets(2:3) * 2 ! two cells in x direction
      !print *,'cell offsets', cell_offsets

      k = 1
      do i = lbound(P%iface,1), ubound(P%iface,1)

        call ind_1d_to_nd(P%da_offsets, P%iface(i), ndidx)
        icell = cStart + ind_nd_to_1d(cell_offsets, ndidx(2:4)) - 1
        call PetscSectionGetFieldOffset(geomSection, icell, i0, geom_offset, ierr); call CHKERR(ierr)
        !print *,'cell', icell, 'centroid', geoms(i1+geom_offset:i3+geom_offset)

        ! First plex cell
        iface = find_face_idx_by_orientation(plex, icell, ndidx(1))
        if (iface.ge.i0) then
          T%iface(k) = iface
          T%albedo(k) = P%albedo(i)
          k = k+1
        endif

        icell = icell + plex%Nlay
        call PetscSectionGetFieldOffset(geomSection, icell, i0, geom_offset, ierr); call CHKERR(ierr)
        !print *,'cell', icell, 'centroid', geoms(i1+geom_offset:i3+geom_offset)

        iface = find_face_idx_by_orientation(plex, icell, ndidx(1))
        if (iface.ge.i0) then
          T%iface(k) = iface
          T%albedo(k) = P%albedo(i)
          k = k+1
        endif
      enddo
    end associate
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief: find a plex face index with a given plex cell and a pprts_face_id (only for rectangular_plex_grids)
  function find_face_idx_by_orientation(plex, icell, fidx) result(iface)
    type(t_plexgrid), intent(in) :: plex
    integer(iintegers), intent(in) :: icell, fidx
    integer(iintegers) :: iface
    integer(iintegers) :: i
    real(ireals) :: face_normal(3)
    integer(iintegers), pointer :: faces_of_cell(:)
    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:)
    integer(mpiint) :: ierr

    iface = -1

    call DMGetSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMPlexGetCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)

    main: block
      select case(fidx)
      case (PPRTS_TOP_FACE)
        iface = faces_of_cell(1)
        exit main

      case (PPRTS_BOT_FACE)
        iface = faces_of_cell(2)
        exit main
      end select

      do i=3,size(faces_of_cell)
        iface = faces_of_cell(i)
        call get_inward_face_normal(iface, icell, geomSection, geoms, face_normal)
        !print *,'iface', iface, 'normal', face_normal

        select case(fidx)
        case (PPRTS_LEFT_FACE)
          if(approx(face_normal(1),  one, .001_ireals)) exit main
        case (PPRTS_RIGHT_FACE)
          if(approx(face_normal(1), -one, .001_ireals)) exit main
        case (PPRTS_REAR_FACE)
          if(approx(face_normal(2),  one, .001_ireals)) exit main
        case (PPRTS_FRONT_FACE)
          if(approx(face_normal(2), -one, .001_ireals)) exit main

        case default
          ierr = -1
          call CHKERR(ierr, 'dont know the faceidx '//toStr(fidx))
        end select
      enddo
      iface = -1
    end block main

    !print *,'looking for', fidx, 'in cell', icell, 'faces', faces_of_cell, '=>', iface

    call DMPlexRestoreCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
  end function

  !> @brief: find plex cell_ids that correspond to a pprts cell (only for rectangular_plex_grids)
  subroutine pprts_cell_to_plex_cell_idx(pprts_C, pprts_cell, plex, plex_cell, ierr)
    type(t_coord), intent(in) :: pprts_C
    integer(iintegers), intent(in) :: pprts_cell(3) ! k,i,j
    type(t_plexgrid), intent(in) :: plex
    integer(iintegers), intent(out) :: plex_cell(2)

    integer(iintegers) :: plex_cell_offsets(3)
    integer(iintegers) :: i, cStart, cEnd
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    call CHKERR(int(pprts_C%glob_zm - plex%Nlay, mpiint), &
      & 'vertical axis of pprts_dmda (Nz='//toStr(pprts_C%glob_zm)//')'// &
      & ' and plex (Nz='//toStr(plex%Nlay)//') dont match')
    call ndarray_offsets([plex%Nlay, 2*pprts_C%glob_xm, pprts_C%glob_ym], plex_cell_offsets)

    associate( k=>pprts_cell(1), i=>pprts_cell(2), j=>pprts_cell(3) )
      plex_cell(1) = ind_nd_to_1d(plex_cell_offsets, [k, 2*i, j], cstyle=.True.)
      plex_cell(2) = ind_nd_to_1d(plex_cell_offsets, [k, 2*i+1, j], cstyle=.True.)
    end associate

    call DMPlexGetDepthStratum(plex%dm, 3_iintegers, cStart, cEnd, ierr); call CHKERR(ierr)
    do i=1,2
      if(.not.is_inrange(plex_cell(i), cStart, cEnd-1)) then
        ierr = ierr+i
        call CHKWARN(ierr, 'when looking for pprts_cell '//toStr(pprts_cell)// &
          & ' found icell ('//toStr(plex_cell(i))//')'// &
          & ' but is not in range cStart/cEnd-1 ('//toStr(cStart)//' - '//toStr(cEnd-1)//')')
      endif
    enddo
    call CHKERR(ierr, 'bad cell idxs found')
  end subroutine
end module
