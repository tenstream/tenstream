module m_icon_plex_utils

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint

  use m_helper_functions, only: chkerr

  use m_plex_grid, only: print_dmplex

  implicit none

  private
  public :: dmplex_2D_to_3D

  logical, parameter :: ldebug=.True.

  contains

    subroutine dmplex_2D_to_3D(dm2d, hhl, dm3d)
      type(tDM), intent(in) :: dm2d
      real(ireals), intent(in) :: hhl(:) ! height levels of interfaces, those will be added to base height of 2D elements
      type(tDM), intent(out) :: dm3d

      integer(mpiint) :: comm, ierr

      call PetscObjectGetComm(dm2d, comm, ierr); call CHKERR(ierr)

      call print_dmplex(comm, dm2d)

    end subroutine
end module
