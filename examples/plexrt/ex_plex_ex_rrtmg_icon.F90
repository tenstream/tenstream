program main
  use m_examples_plex_ex_rrtmg_icon
  implicit none

  character(len=default_str_len) :: gridfile, icondatafile, outfile
  logical :: lflg
  integer(mpiint) :: myid, ierr
  character(len=10*default_str_len) :: default_options
  real(ireals) :: Ag
  logical :: lthermal, lsolar

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
  call init_mpi_data_parameters(PETSC_COMM_WORLD)
  call read_commandline_options(PETSC_COMM_WORLD)
  call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid', gridfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a grid filename... please call with -grid <fname_of_icon_gridfile.nc>'

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-data', icondatafile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a icondata filename... please call with -data <fname_of_icondatafile.nc>'

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

  Ag = .1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

  lsolar = .True.
  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)

  default_options=''
  default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
  default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

  !default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_abso_direct hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -show_abso_diffuse hdf5:'//trim(outfile)//'::append'

  default_options=trim(default_options)//' -plexrt_dump_thermal_Edn_2_ke1 hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_thermal_Eup_2_ke1 hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_thermal_abso hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_Edir_2_ke1 hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_Edn_2_ke1 hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_Eup_2_ke1 hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_abso hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_lwc hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_iwc hdf5:'//trim(outfile)//'::append'
  default_options=trim(default_options)//' -plexrt_dump_temp hdf5:'//trim(outfile)//'::append'

  !default_options=trim(default_options)//' -dump_optprop_kabs hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -dump_optprop_ksca hdf5:'//trim(outfile)//'::append'
  !default_options=trim(default_options)//' -dump_optprop_g hdf5:'//trim(outfile)//'::append'

  default_options=trim(default_options)//' -show_fV2cV_level_heights_vec hdf5:lvl_'//trim(outfile)//''

  if(myid.eq.0) print *,'Adding default Petsc Options:', trim(default_options)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

  call plex_ex_rrtmg_icon(PETSC_COMM_WORLD, gridfile, icondatafile, Ag, lthermal, lsolar)

  call mpi_barrier(PETSC_COMM_WORLD, ierr)
  call PetscFinalize(ierr)
end program
