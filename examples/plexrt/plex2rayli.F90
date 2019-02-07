module m_ex_plex2rayli

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    default_str_len, init_mpi_data_parameters, &
    zero, one, i0

  use m_helper_functions, only: CHKERR, reverse, itoa, cstr

  use m_plex2rayli, only: dm3d_to_rayli_dmplex

  use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D

  implicit none

  contains
    subroutine fish_plex_to_rayli(comm, Nx, Ny, Nz, dz)
      MPI_Comm, intent(in)           :: comm
      integer(iintegers), intent(in) :: Nx, Ny, Nz
      real(ireals), intent(in)       :: dz

      type(tDM)          :: dm2d, dm3d, dmrayli
      real(ireals)       :: hhl(Nz)
      integer(iintegers) :: k
      integer(iintegers), allocatable :: zindex(:)
      integer(mpiint)    :: myid, numnodes, ierr

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
      call create_2d_fish_plex(Nx, Ny, dm2d)

      hhl(1) = zero
      do k=2,Nz
        hhl(k) = hhl(k-1) + dz
      enddo
      hhl = reverse(hhl)

      call dmplex_2D_to_3D(dm2d, Nz, hhl, dm3d, zindex)
      call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, "-show_dm3d", ierr); call CHKERR(ierr)
      call dm3d_to_rayli_dmplex(dm3d, dmrayli)
      call PetscObjectViewFromOptions(dmrayli, PETSC_NULL_DM, "-show_rayli", ierr); call CHKERR(ierr)

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm3d, ierr); call CHKERR(ierr)
      call DMDestroy(dmrayli, ierr); call CHKERR(ierr)
    end subroutine
end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_ex_plex2rayli, only: fish_plex_to_rayli
  use m_data_parameters, only : ireals, iintegers, mpiint, default_str_len, init_mpi_data_parameters
  use m_helper_functions, only: CHKERR
  use m_tenstream_options, only : read_commandline_options
  implicit none

  character(len=default_str_len) :: outfile
  logical :: lflg
  integer(mpiint) :: ierr
  character(len=10*default_str_len) :: default_options
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dz

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
  call init_mpi_data_parameters(PETSC_COMM_WORLD)
  call read_commandline_options(PETSC_COMM_WORLD)

  Nx = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
  Ny = 3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
  Nz = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg,ierr) ; call CHKERR(ierr)
  dz = 1._ireals/Nz
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

  default_options='-polar_coords no'
  default_options=trim(default_options)//' -show_rayli hdf5:'//trim(outfile)
  !default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'

  print *,'Adding default Petsc Options:', trim(default_options)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

  call fish_plex_to_rayli(PETSC_COMM_WORLD, Nx, Ny, Nz, dz)

  call mpi_barrier(PETSC_COMM_WORLD, ierr)
  call PetscFinalize(ierr)
end program
