module m_buildings
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: iintegers, ireals, mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & toStr, &
    & ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets, &
    & get_arg

  implicit none

  private
  public :: &
    & t_pprts_buildings, t_plex_buildings,&
    & init_buildings,                     &
    & faceidx_by_cell_plus_offset,        &
    & check_buildings_consistency,        &
    & PPRTS_TOP_FACE,                     &
    & PPRTS_BOT_FACE,                     &
    & PPRTS_LEFT_FACE,                    &
    & PPRTS_RIGHT_FACE,                   &
    & PPRTS_REAR_FACE,                    &
    & PPRTS_FRONT_FACE

  type, abstract :: t_buildings
    real(ireals), allocatable :: albedo(:)      ! albedo on a face, dim(Nbuilding_faces)
  end type

  type, extends(t_buildings) :: t_pprts_buildings
    ! local face indices in pprts mesh,
    ! use faceidx_by_cell_plus_offset to translate from cell idx + face_ids (from PPRTS_XXX_FACE)
    integer(iintegers), allocatable :: iface(:)

    ! offsets from iface idx sets to da indices, da_size e.g. [6, zm, xm, ym]
    integer(iintegers), allocatable :: da_offsets(:)

    ! irradiances on surfaces, if allocated -> dim(size(iface))
    real(ireals), allocatable :: edir(:), incoming(:), outgoing(:)
  end type

  type, extends(t_buildings) :: t_plex_buildings
    integer(iintegers), allocatable :: iface(:)! face indices dim(Nbuilding_faces) valid in range [fStart, fEnd-1]
  end type

  integer, parameter :: &
    & PPRTS_TOP_FACE  =1, & ! z+0
    & PPRTS_BOT_FACE  =2, & ! z+1
    & PPRTS_LEFT_FACE =3, & ! x+0
    & PPRTS_RIGHT_FACE=4, & ! x+1
    & PPRTS_REAR_FACE =5, & ! y+0
    & PPRTS_FRONT_FACE=6    ! y+1

contains

  subroutine init_buildings(buildings, da_sizes, ierr)
    type(t_pprts_buildings), intent(inout), allocatable :: buildings
    integer(iintegers), intent(in) :: da_sizes(:)
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if(allocated(buildings)) ierr = 1
    call CHKERR(ierr, 'buildings already allocated')

    allocate(buildings)
    allocate(buildings%da_offsets(size(da_sizes)))
    call ndarray_offsets(da_sizes, buildings%da_offsets)
  end subroutine

  !> @brief: find face index in a pprts dmda
  !> @details: this a helper routine to determine a flattened index in a 3D DMDA domain as in pprts
  !> \n        k,i,j are cell indices in 3D starting to count at 1
  !> \n        face_id has to be one of [PPRTS_TOP_FACE, ... ,PPRTS_FRONT_FACE]
  function faceidx_by_cell_plus_offset(da_offsets, k, i, j, face_id) result(fidx)
    integer(iintegers), intent(in) :: da_offsets(:)
    integer(iintegers), intent(in) :: k, i, j, face_id
    integer(iintegers) :: fidx
    fidx = ind_nd_to_1d(da_offsets, [face_id, k, i, j])
  end function

  subroutine check_buildings_consistency(buildings, Nz, Nx, Ny, ierr)
    type(t_pprts_buildings), intent(inout), allocatable :: buildings
    integer(iintegers) :: Nx, Ny, Nz ! domain size, i.e. cell count
    integer(mpiint), intent(out) :: ierr

    logical, allocatable :: sides(:,:,:,:)
    integer(iintegers) :: m, idx(4), i,j,k

    ierr = 0

    call CHKERR(int(buildings%da_offsets(2)-6_iintegers, mpiint), &
      & 'da_offsets(2) should be 6 as in 6 sides but found '//toStr(buildings%da_offsets(2)))
    call CHKERR(int(buildings%da_offsets(3)-6*Nz, mpiint), &
      & 'da_offsets(3) should be 6*Nz but found '//toStr(buildings%da_offsets(3))// &
      & ' instead of '//toStr(Nz*6))

    call CHKERR(int(buildings%da_offsets(4)-6*Nz*Nx, mpiint), &
      & 'da_offsets(3) should be 6*Nz*Nx but found '//toStr(buildings%da_offsets(4))// &
      & ' instead of '//toStr(6*Nz*Nx))

    allocate(sides(6, Nz, Nx, Ny))
    sides = .False.

    do m = 1, size(buildings%iface)
      call ind_1d_to_nd(buildings%da_offsets, buildings%iface(m), idx)
      associate(d => idx(1), k => idx(2), i => idx(3), j => idx(4))
        sides(d,k,i,j) = .True.
      end associate
    enddo

    do j = 1, Ny
      do i = 1, Nx
        do k = 1, Nz
          if(any(sides(:,k,i,j))) then
            if(.not.all(sides(:,k,i,j))) then
              ierr = ierr + 1
              print *,'in cell k,i,j:', k,i,j, &
                & ' we found an semi open box'// new_line('')// &
                & ' - this may work for direct radiation or rayli results'// &
                & ' but diffuse radiation in the TenStream is probably wrong!', &
                & ' sides:', sides(:,k,i,j)
            endif
          endif
        enddo
      enddo
    enddo
  end subroutine
end module
