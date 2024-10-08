program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only: iintegers, mpiint, ireals
  use m_example_pprts_rrtm_iterations, only: example_rrtm_lw_sw
  use m_helper_functions, only: get_petsc_opt

  implicit none

  integer(mpiint) :: ierr, myid
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dx, dy
  logical :: lflg

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  Nx = 3; Ny = 3; Nz = 5
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr)

  dx = 500
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr)
  dy = dx
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr)

  call example_rrtm_lw_sw(Nx, Ny, Nz, dx, dy)

  call mpi_finalize(ierr)
end program
