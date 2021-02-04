program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only : iintegers, mpiint, ireals
  use m_example_pprts_rrtm_lw_sw, only: ex_pprts_rrtm_lw_sw

  implicit none

  integer(mpiint) :: ierr, myid
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals)       :: dx, dy, phi0, theta0, albedo_th, albedo_sol
  logical :: lflg

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

  call read_commandline_options(PETSC_COMM_WORLD)

  Nx=3; Ny=3; Nz=5
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr)

  dx = 500
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr)
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr)

  phi0 = 180._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr)
  theta0 = 0._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr)
  albedo_th = 0._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", albedo_th, lflg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag_thermal", albedo_th, lflg, ierr)
  albedo_sol = .1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", albedo_sol, lflg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag_solar", albedo_sol, lflg, ierr)

  if (myid.eq.0) print *,'Running rrtm_lw_sw example with grid size:', Nx, Ny, Nz

  call ex_pprts_rrtm_lw_sw(Nx, Ny, Nz, dx, dy, phi0, theta0, albedo_th, albedo_sol)

  call mpi_finalize(ierr)
end program
