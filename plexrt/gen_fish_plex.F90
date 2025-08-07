module m_gen_fish_plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters, &
                               zero, one, i0, i1, i2, i3, i4, i5

  use m_helper_functions, only: CHKERR, get_petsc_opt

  use m_icon_plex_utils, only: create_2d_fish_plex, &
                               dmplex_2D_to_3D, dump_ownership

  implicit none

  integer(iintegers), parameter :: TOP_BOT_FACE = 1, SIDE_FACE = 2

  PetscErrorCode :: ierr

contains
  subroutine init()
    type(tDM) :: dm2d, dm2d_dist, dm3d
    character(len=*), parameter :: default_options = '&
      & -default_option_show_plex hdf5:fish.h5 &
      & -default_option_show_plex3d hdf5:fish3d.h5'

    real(ireals), parameter :: hhl(2) = [zero, one]
    integer(iintegers), allocatable :: zindex(:)
    integer(mpiint) :: myid, numnodes
    integer(iintegers) :: Nx, Ny
    PetscBool :: lflg

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr); call CHKERR(ierr)

    call init_mpi_data_parameters(PETSC_COMM_WORLD)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
    if (lflg .eqv. PETSC_FALSE) Nx = 2
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
    if (lflg .eqv. PETSC_FALSE) Ny = 3

    call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(PETSC_COMM_WORLD, numnodes, ierr); call CHKERR(ierr)

    call create_2d_fish_plex(PETSC_COMM_WORLD, Nx, Ny, dm2d, dm2d_dist)
    call dmplex_2D_to_3D(dm2d_dist, size(hhl, kind=iintegers), hhl, [zero, zero, -huge(zero) * 1e-1_ireals], dm3d, zindex)

    call PetscObjectViewFromOptions(PetscObjectCast(dm2d_dist), PETSC_NULL_OBJECT, "-default_option_show_plex", ierr)
    call CHKERR(ierr)
    call PetscObjectViewFromOptions(PetscObjectCast(dm2d_dist), PETSC_NULL_OBJECT, "-show_plex", ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(PetscObjectCast(dm3d), PETSC_NULL_OBJECT, "-default_option_show_plex3d", ierr)
    call CHKERR(ierr)
    call PetscObjectViewFromOptions(PetscObjectCast(dm3d), PETSC_NULL_OBJECT, "-show_plex3d", ierr); call CHKERR(ierr)

    call dump_ownership(dm3d, '-show_plex3d_ownership')

    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)
    call DMDestroy(dm3d, ierr); call CHKERR(ierr)
    call PetscFinalize(ierr)
  end subroutine

end module

program main
  use m_gen_fish_plex

  call init()

end program
