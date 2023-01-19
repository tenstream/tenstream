program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi, only: MPI_COMM_WORLD

  use m_data_parameters, only: &
    & default_str_len, &
    & finalize_mpi, &
    & iintegers, &
    & init_mpi_data_parameters, &
    & ireals, &
    & mpiint, &
    & share_dir

  use m_examples_pprts_specint_tree, only: ex_pprts_specint_tree
  use m_helper_functions, only: CHKERR, toStr, get_petsc_opt
  use m_netcdfio, only: ncwrite
  use m_tenstream_options, only: read_commandline_options

  implicit none

  character(len=default_str_len) :: outfile, atm_filename, specint
  integer(iintegers) :: Nx, Ny, Nlay, icollapse, Ntree_height
  real(ireals) :: dx, dy
  real(ireals) :: phi0, theta0
  real(ireals) :: Ag_solar, Ag_thermal
  real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays on rank 0

  character(len=10*default_str_len) :: rayli_options
  character(len=default_str_len) :: groups(2), dimnames(3)
  logical :: lflg, lverbose, lrayli_opts, lsolar, lthermal, lfile_exists, lhave_outfile, luse_usgs_db
  integer(mpiint) :: comm, myid, ierr

  call mpi_init(ierr)
  comm = mpi_comm_world
  call init_mpi_data_parameters(comm)
  call read_commandline_options(comm)

  specint = 'no_default_set'
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-specint", specint, lflg, ierr); call CHKERR(ierr)

  atm_filename = share_dir//'tenstream_default.atm'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', &
    & atm_filename, lflg, ierr); call CHKERR(ierr)
  inquire (file=trim(atm_filename), exist=lfile_exists)
  if (.not. lfile_exists) then
    ierr = 1
  else
    ierr = 0
  end if
  call CHKERR(ierr, 'background atmosphere file: `'//trim(atm_filename)//&
    & '` does not exist! Please provide a path with option -atm <atmfile>')

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', outfile, lhave_outfile, ierr); call CHKERR(ierr)

  lsolar = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-solar', &
                     lsolar, lflg, ierr); call CHKERR(ierr)

  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-thermal', &
                     lthermal, lflg, ierr); call CHKERR(ierr)

  Nx = 5; Ny = 5; Nlay = 6
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Nz", Nlay, lflg, ierr); call CHKERR(ierr)
  icollapse = -1
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-icollapse", icollapse, lflg, ierr); call CHKERR(ierr)
  Ntree_height = 4
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ntree", Ntree_height, lflg, ierr); call CHKERR(ierr)

  luse_usgs_db = .false.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-usgs', &
                     luse_usgs_db, lflg, ierr); call CHKERR(ierr)

  dx = 100
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dy = dx
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)

  phi0 = 180._ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 0._ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  Ag_solar = 0.15_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_solar", Ag_solar, lflg, ierr); call CHKERR(ierr)

  Ag_thermal = 0.05_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag_thermal", Ag_thermal, lflg, ierr); call CHKERR(ierr)

  lverbose = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-verbose', &
                     lverbose, lflg, ierr); call CHKERR(ierr)

  lrayli_opts = .false.
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-rayli_opts', &
                     lrayli_opts, lflg, ierr); call CHKERR(ierr)

  if (lrayli_opts) then
    rayli_options = ''
    rayli_options = trim(rayli_options)//' -solver rayli'
    rayli_options = trim(rayli_options)//' -rayli_cyclic_bc'
    rayli_options = trim(rayli_options)//' -show_rayli_dm3d hdf5:dm.h5'

    rayli_options = trim(rayli_options)//' -rayli_snapshot'
    rayli_options = trim(rayli_options)//' -rayli_snap_Nx 256'
    rayli_options = trim(rayli_options)//' -rayli_snap_Ny 256'
    rayli_options = trim(rayli_options)//' -visit_image_zoom .30'
    rayli_options = trim(rayli_options)//' -visit_parallel_scale 291.5'
    rayli_options = trim(rayli_options)//' -visit_focus 500,500,0'
    rayli_options = trim(rayli_options)//' -visit_view_normal -0.2811249083446944,-0.7353472951470268,0.6166304739697339'
    rayli_options = trim(rayli_options)//' -visit_view_up 0.1878717450780742,0.5879401184069877,0.7867849925925738'

    if (lverbose) print *, 'Adding rayli Petsc Options:', trim(rayli_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, trim(rayli_options), ierr); call CHKERR(ierr)
  end if

  call ex_pprts_specint_tree(   &
    & specint,                  &
    & comm, lverbose,           &
    & lthermal, lsolar,         &
    & Nx, Ny, Nlay,             &
    & dx, dy,                   &
    & atm_filename,             &
    & phi0, theta0,             &
    & Ag_solar, Ag_thermal,     &
    & Ntree_height,             &
    & luse_usgs_db,             &
    & gedir, gedn, geup, gabso, &
    & icollapse=icollapse)

  if (lhave_outfile) then
    groups(1) = trim(outfile)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    if (myid .eq. 0_mpiint) then
      dimnames = [character(len=default_str_len) :: 'zlev', 'x', 'y']
      if (lsolar) then
        groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames); call CHKERR(ierr)
      end if
      groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames); call CHKERR(ierr)
      groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames); call CHKERR(ierr)
      dimnames = [character(len=default_str_len) :: 'zlay', 'x', 'y']
      groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames); call CHKERR(ierr)
    end if
  end if

  call finalize_mpi(ierr, .true., .true.)
end program

