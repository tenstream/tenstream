program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only: iintegers, mpiint, ireals, default_str_len, share_dir
  use m_example_pprts_specint_lw_sw, only: ex_pprts_specint_lw_sw
  use m_helper_functions, only: CHKERR, get_petsc_opt, deallocate_allocatable
  use m_netcdfIO, only: ncwrite

  implicit none

  integer(mpiint) :: ierr, comm, myid
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dx, dy, phi0, theta0, albedo_th, albedo_sol, vlwc, viwc
  character(len=default_str_len) :: atm_filename, outfile, specint
  logical :: lthermal, lsolar, lflg, lhave_outfile
  real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso

  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  call mpi_comm_rank(comm, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call read_commandline_options(comm)

  specint = 'no_default_set'
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-specint", specint, lflg, ierr); call CHKERR(ierr)

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

  atm_filename = share_dir//'tenstream_default.atm'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-solar", lsolar, lflg, ierr); call CHKERR(ierr)

  vlwc = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-lwc", vlwc, lflg, ierr); call CHKERR(ierr)

  viwc = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-iwc", viwc, lflg, ierr); call CHKERR(ierr)

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', outfile, lhave_outfile, ierr); call CHKERR(ierr)
  ! if(.not.lhave_outfile) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  if (myid .eq. 0) print *, 'Running repwvl_lw_sw example with grid size:', Nx, Ny, Nz

  if (lhave_outfile) then
    call ex_pprts_specint_lw_sw(specint, &
      & comm, &
      & Nx, Ny, Nz, dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & gedir, gedn, geup, gabso, &
      & vlwc, viwc, &
      & outfile=outfile )
  else
    call ex_pprts_specint_lw_sw(specint, &
      & comm, &
      & Nx, Ny, Nz, dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & gedir, gedn, geup, gabso, &
      & vlwc, viwc )
  endif

  call deallocate_allocatable(gedir)
  call deallocate_allocatable(gedn)
  call deallocate_allocatable(geup)
  call deallocate_allocatable(gabso)

  call mpi_finalize(ierr)
end program
