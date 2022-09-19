program main
#include "petsc/finclude/petsc.h"
  use petsc

  use m_examples_plex_specint_fish, only: ex_plex_specint_fish

  use m_data_parameters, only: &
    & ireals, iintegers, mpiint, &
    & default_str_len, &
    & init_mpi_data_parameters

  use m_helper_functions, only: &
    & CHKERR, &
    & angle_between_two_vec, &
    & deg2rad, rad2deg, &
    & get_petsc_opt, &
    & spherical_2_cartesian, &
    & rotation_matrix_around_axis_vec

  use m_tenstream_options, only: read_commandline_options

  implicit none

  character(len=default_str_len) :: atm_filename, outfile, specint
  logical :: lflg, lthermal, lsolar, lregular_mesh, lverbose, ladd_rayli_opts
  integer(mpiint) :: comm, ierr
  character(len=10*default_str_len) :: default_options, rayli_options
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dx, dz, Ag_th, Ag_sol, sundir(3), lwc, phi, theta
  real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)

  default_options = ''
  !default_options = trim(default_options)//' -plexrt_view_geometry'
  !default_options=trim(default_options)//' -plexrt_view_optprop'
  !default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_iwc hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

  print *, 'Adding default Petsc Options:', trim(default_options)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr); call CHKERR(ierr)

  rayli_options = '-use_regular_mesh -plexrt_use_rayli -rayli_cyclic_bc'
  ladd_rayli_opts = .false.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-rayli_opts', &
                     ladd_rayli_opts, lflg, ierr); call CHKERR(ierr)
  if (ladd_rayli_opts) then
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, rayli_options, ierr); call CHKERR(ierr)
  end if

  comm = PETSC_COMM_WORLD
  call init_mpi_data_parameters(comm)
  call read_commandline_options(comm)

  specint = 'no_default_set'
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-specint", specint, lflg, ierr); call CHKERR(ierr)

  Nx = 2
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  Ny = 3
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  Nz = 2
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr); call CHKERR(ierr)
  dx = 1e2_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dz = 1e2_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr); call CHKERR(ierr)
  Ag_th = .1_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_th", Ag_th, lflg, ierr); call CHKERR(ierr)
  Ag_sol = .1_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_sol", Ag_sol, lflg, ierr); call CHKERR(ierr)

  phi = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-phi", phi, lflg, ierr); call CHKERR(ierr)
  theta = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-theta", theta, lflg, ierr); call CHKERR(ierr)

  sundir = spherical_2_cartesian(phi, theta)

  lwc = 0
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-lwc", lwc, lflg, ierr); call CHKERR(ierr)

  lverbose = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-verbose', lverbose, lflg, ierr); call CHKERR(ierr)

  lregular_mesh = .false.
  call get_petsc_opt(PETSC_NULL_CHARACTER, &
                     '-use_regular_mesh', lregular_mesh, lflg, ierr); call CHKERR(ierr)

  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-thermal', lthermal, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-solar', lsolar, lflg, ierr); call CHKERR(ierr)

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need to supply a atm filename... please call with -atm <fname_of_atm_file.dat>')

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

  call ex_plex_specint_fish(&
      & specint, &
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_th, &
      & Ag_sol, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

  call PetscFinalize(ierr)
end program
