module m_gen_fish_plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, init_mpi_data_parameters, &
    zero, i0, i1, i2, i3, i4, i5

  use m_helper_functions, only: itoa, CHKERR

  use m_icon_plex_utils, only: create_2d_fish_plex

  implicit none

  integer(iintegers), parameter :: TOP_BOT_FACE=1, SIDE_FACE=2

  PetscErrorCode :: ierr

  contains
    subroutine init()
      type(tDM) :: dm
      character(len=*), parameter :: default_options = '&
        & -default_option_show_plex hdf5:fish.h5'

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

      call create_2d_fish_plex(dm, Nx, Ny, .False.)

      call PetscObjectViewFromOptions(dm, PETSC_NULL_VEC, "-default_option_show_plex", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm, PETSC_NULL_VEC, "-show_plex", ierr); call CHKERR(ierr)

      call DMDestroy(dm, ierr);CHKERRQ(ierr)
      call PetscFinalize(ierr)
    end subroutine


  end module

  program main
    use m_gen_fish_plex

    call init()

  end program
