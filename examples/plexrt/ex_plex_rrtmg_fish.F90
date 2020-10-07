program main
  use m_examples_plex_rrtmg_fish
  implicit none

  character(len=default_str_len) :: outfile
  logical :: lflg, lthermal, lsolar
  integer(mpiint) :: ierr
  character(len=10*default_str_len) :: default_options
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dx, dz, Ag


  !character(len=*),parameter :: ex_out='plex_ex_dom1_out.h5'
  !character(len=*),parameter :: ex_out='plex_test_out.h5'
  !character(len=*),parameter :: lwcfile='lwc_ex_24_3.nc'
  !character(len=*),parameter :: lwcfile='lwc_ex_dom1.nc'

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
  call init_mpi_data_parameters(PETSC_COMM_WORLD)
  call read_commandline_options(PETSC_COMM_WORLD)

  Nx = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
  Ny = 3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
  Nz = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg,ierr) ; call CHKERR(ierr)
  dx = 1e3_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg,ierr) ; call CHKERR(ierr)
  dz = 1e2_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)
  Ag = .1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)
  lsolar = .True.
  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)


  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

  default_options='-polar_coords no'
  default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
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

  print *,'Adding default Petsc Options:', trim(default_options)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

  call ex_plex_rrtmg_fish(PETSC_COMM_WORLD, Nx, Ny, Nz, dx, dz, Ag, lthermal, lsolar)

  call mpi_barrier(PETSC_COMM_WORLD, ierr)
  call PetscFinalize(ierr)
end program
