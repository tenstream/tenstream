module m_buildings
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: iintegers, ireals, mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & ind_nd_to_1d, ndarray_offsets, &
    & get_arg

  implicit none

  private
  public :: &
    & t_pprts_buildings, t_plex_buildings,&
    & init_buildings,                     &
    & faceidx_by_cell_plus_offset,        &
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

    ! offsets from iface idx sets to da indices, i.e. [6, zm, xm, ym]
    integer(iintegers), allocatable :: da_offsets(:)
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
end module
