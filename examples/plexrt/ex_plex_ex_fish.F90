program main
  use m_examples_plex_ex_fish
  implicit none

  character(len=default_str_len) :: outfile
  logical :: lflg
  integer(mpiint) :: comm, ierr
  character(len=10*default_str_len) :: default_options
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals) :: dz, Ag, sundir(3), dtau, w0, g, B0
  logical :: lverbose, lregular_mesh, lthermal, lsolar
  real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)
  comm = PETSC_COMM_WORLD

  call init_mpi_data_parameters(comm)
  call read_commandline_options(comm)

  Nx = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  Ny = 3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  Nz = 2
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr); call CHKERR(ierr)
  dz = one / real(Nz, ireals)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr); call CHKERR(ierr)
  Ag = .1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr); call CHKERR(ierr)
  dtau = one
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-tau", dtau, lflg, ierr); call CHKERR(ierr)
  w0 = zero
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-w0", w0, lflg, ierr); call CHKERR(ierr)
  g = zero
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-g", g, lflg, ierr); call CHKERR(ierr)
  B0 = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-B0", B0, lflg, ierr); call CHKERR(ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

  lverbose = .true.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-verbose', &
                           lverbose, lflg, ierr); call CHKERR(ierr)

  lregular_mesh = .false.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                           '-use_regular_mesh', lregular_mesh, lflg, ierr); call CHKERR(ierr)

  lthermal = .true.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-thermal', &
                           lthermal, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-solar', &
                           lsolar, lflg, ierr); call CHKERR(ierr)

  call init_sundir()

  default_options = ''
  default_options = trim(default_options)//' -show_plex hdf5:'//trim(outfile)
  default_options = trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_iwc hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
  default_options = trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

  if (lverbose) print *, 'Adding default Petsc Options:', trim(default_options)
  !call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

  call plex_ex_fish(comm, &
    & lverbose, &
    & lthermal, lsolar, &
    & lregular_mesh, &
    & Nx, Ny, Nz, &
    & dz, Ag, &
    & sundir, &
    & dtau, w0, g, B0, &
    & edir, edn, eup, abso)

  call PetscFinalize(ierr)

contains
  subroutine init_sundir()
    use m_helper_functions, only: cross_3d, &
      & rotation_matrix_world_to_local_basis, &
      & rotation_matrix_local_basis_to_world, &
      & rotation_matrix_around_axis_vec, &
      & rotate_angle_x

    logical :: lflg
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: nargs
    real(ireals) :: first_normal(3)
    real(ireals) :: rot_angle, Mrot(3, 3), U(3), rot_sundir(3)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    first_normal(:) = [0, 0, -1]

    if (lverbose .and. myid .eq. 0) print *, myid, 'determine initial sundirection ...'
    nargs = i3
    call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                                  "-sundir", sundir, nargs, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      call CHKERR(int(nargs - i3, mpiint), 'must provide exactly 3 values for -sundir. '// &
                  'Need to be given comma separated without spaces')
    else
      sundir = first_normal
    end if
    sundir = sundir / norm2(sundir)

    rot_angle = zero
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_phi", &
                             rot_angle, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), first_normal)
      rot_sundir = matmul(Mrot, sundir)
      if (lverbose .and. myid .eq. 0) print *, 'rot_sundir', rot_sundir
      if (myid .eq. 0) &
        print *, 'rotated sundirection = ', rot_sundir, &
        & ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)), 'deg'
      sundir = rot_sundir
    end if

    rot_angle = zero
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_theta", &
                             rot_angle, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      U = cross_3d(first_normal, sundir)
      Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), U)
      rot_sundir = matmul(Mrot, sundir)
      if (lverbose .and. myid .eq. 0) print *, 'S', sundir, norm2(sundir)
      if (lverbose .and. myid .eq. 0) print *, 'U', U, norm2(U)
      if (lverbose .and. myid .eq. 0) print *, 'rot_sundir', rot_sundir
      if (lverbose .and. myid .eq. 0) &
        print *, 'rotated sundirection = ', rot_sundir, ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)), 'deg'
      sundir = rot_sundir
    end if

    if (lverbose .and. myid .eq. 0) print *, 'determine initial sundirection ... done'
  end subroutine
end program
