program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only: iintegers, mpiint, ireals, default_str_len
  use m_example_pprts_rrtm_lw_sw, only: ex_pprts_rrtm_lw_sw
  use m_helper_functions, only: CHKERR, get_petsc_opt, deallocate_allocatable

  implicit none

  integer(mpiint) :: ierr, comm, myid
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dx, dy, phi0, theta0, albedo_th, albedo_sol, vlwc, viwc
  character(len=default_str_len) :: atm_filename
  logical :: lthermal, lsolar, lflg
  real(ireals), allocatable, dimension(:, :, :) :: fdir, fdn, fup, fdiv

  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  call mpi_comm_rank(comm, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call read_commandline_options(comm)

  Nx = 3; Ny = 3; Nz = 5
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr); call CHKERR(ierr)

  dx = 500
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dy = dx
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)

  phi0 = 180._ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 0._ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)
  albedo_th = 0._ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag", albedo_th, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_thermal", albedo_th, lflg, ierr); call CHKERR(ierr)
  albedo_sol = .1_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag", albedo_sol, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_solar", albedo_sol, lflg, ierr); call CHKERR(ierr)

  atm_filename = 'atm.dat'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-solar", lsolar, lflg, ierr); call CHKERR(ierr)

  vlwc = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-lwc", vlwc, lflg, ierr); call CHKERR(ierr)

  viwc = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-iwc", viwc, lflg, ierr); call CHKERR(ierr)

  if (myid .eq. 0) print *, 'Running rrtm_lw_sw example with grid size:', Nx, Ny, Nz

  call ex_pprts_rrtm_lw_sw(comm, &
    & Nx, Ny, Nz, dx, dy, &
    & phi0, theta0, albedo_th, albedo_sol, &
    & lthermal, lsolar, atm_filename, &
    & fdir, fdn, fup, fdiv, &
    & vlwc, viwc)

  call deallocate_allocatable(fdir)
  call deallocate_allocatable(fdn)
  call deallocate_allocatable(fup)
  call deallocate_allocatable(fdiv)

  call mpi_finalize(ierr)
end program
