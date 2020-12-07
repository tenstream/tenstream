module m_gen_regular_plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, init_mpi_data_parameters, &
    zero, one, i0, i1, i2, i3, i4, i5

  use m_helper_functions, only: itoa, CHKERR

  use m_icon_plex_utils, only: create_2d_regular_plex, &
    dmplex_2D_to_3D, dump_ownership

  implicit none

  integer(iintegers), parameter :: TOP_BOT_FACE=1, SIDE_FACE=2

  PetscErrorCode :: ierr

  contains
    subroutine init()
      type(tDM) :: dm2d, dm2d_dist, dm3d
      character(len=*), parameter :: default_options = '&
        & -default_option_show_plex hdf5:regular.h5 &
        & -default_option_show_plex3d hdf5:regular3d.h5'

      real(ireals), parameter :: hhl(2) = [zero, one]
      integer(iintegers), allocatable :: zindex(:)
      integer(mpiint) :: myid, numnodes
      integer(iintegers) :: Nx, Ny
      PetscBool :: lflg

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHKERRQ(ierr)
      call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr); call CHKERR(ierr)

      call init_mpi_data_parameters(PETSC_COMM_WORLD)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
      if(lflg.eqv.PETSC_FALSE) Nx = 2
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
      if(lflg.eqv.PETSC_FALSE) Ny = 3

      call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(PETSC_COMM_WORLD, numnodes, ierr); call CHKERR(ierr)

      call create_2d_regular_plex(PETSC_COMM_WORLD, Nx, Ny, dm2d, dm2d_dist)
      call dmplex_2D_to_3D(dm2d_dist, size(hhl, kind=iintegers), hhl, [zero,zero,-huge(zero)], dm3d, zindex)

      call PetscObjectViewFromOptions(dm2d_dist, PETSC_NULL_DM, "-default_option_show_plex", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm2d_dist, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, "-default_option_show_plex3d", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, "-show_plex3d", ierr); call CHKERR(ierr)

      call dump_ownership(dm3d, '-show_plex3d_ownership')

      call DMDestroy(dm2d, ierr);call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr);call CHKERR(ierr)
      call DMDestroy(dm3d, ierr);call CHKERR(ierr)
      call PetscFinalize(ierr)
    end subroutine


  end module

  program main
    use m_gen_regular_plex

    call init()

  end program
