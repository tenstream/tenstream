module m_vegetation_optprop

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & toStr, &
    & imp_bcast, &
    & linspace, &
    & get_arg, &
    & split

  use m_tenstream_interpolation, only: interp_1d

  use m_search, only: find_real_location

  implicit none

  type veg_type
    character(len=default_str_len) :: vname
    real(ireals), allocatable :: lambda(:) ! wavelength of sample pts in [m]
    real(ireals), allocatable :: albedo(:) ! albedo at lambda sample point
  end type

  type(veg_type), allocatable :: veg(:)

contains

  subroutine init_vegetation_types_simple(ierr)
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if (allocated(veg)) return

    allocate (veg(3))

    veg(1)%vname = "bark" ! from usgs WhitebarkPine YNP-WB-1 frst AVIRISb RTGC
    veg(1)%lambda = [real(ireals) :: 0.419, 0.547, 0.676, 0.695, 0.743, 0.772, 0.915, 0.954, 1.069, 1.175,&
      & 1.264, 1.304, 1.424, 1.483, 1.663, 1.812, 1.99, 2.249, 2.478] * 1e-6_ireals
    veg(1)%albedo = [real(ireals) :: 0.013, 0.029, 0.023, 0.033, 0.113, 0.128, 0.144, 0.138, 0.158, 0.129,&
      & 0.146, 0.142, 0.057, 0.056, 0.091, 0.078, 0.043, 0.052, 0.033]

    veg(2)%vname = "grass" ! from usgs LawnGrass GDS91b shifted 3nm  BECKa AREF
    veg(2)%lambda = [real(ireals) :: &
                     0.205, 0.403, 0.499, 0.515, 0.527, 0.543, 0.563, 0.591, 0.679, 0.688, &
                     0.694, 0.708, 0.714, 0.72, 0.736, 0.746, 0.754, 0.76, 0.775, 0.807, &
                     0.914, 0.933, 0.957, 0.976, 1.068, 1.104, 1.128, 1.144, 1.154, 1.179, &
                     1.204, 1.258, 1.299, 1.318, 1.353, 1.374, 1.398, 1.404, 1.408, 1.423, &
                     1.448, 1.474, 1.592, 1.634, 1.676, 1.716, 1.772, 1.835, 1.855, 1.865, &
                     1.885, 1.895, 1.915, 1.935, 2.015, 2.125, 2.215, 2.265, 2.466, 2.656, &
                     2.688, 2.752, 2.784, 2.816, 2.944, 2.976] * 1e-6_ireals
    veg(2)%albedo = [real(ireals) :: &
                     0.02, 0.027, 0.041, 0.06, 0.085, 0.095, 0.089, 0.064, 0.039, 0.049, &
                     0.072, 0.172, 0.224, 0.287, 0.485, 0.588, 0.637, 0.66, 0.686, 0.7, &
                     0.704, 0.695, 0.663, 0.659, 0.699, 0.693, 0.671, 0.619, 0.593, 0.574, &
                     0.572, 0.594, 0.575, 0.537, 0.452, 0.395, 0.224, 0.197, 0.181, 0.154, &
                     0.142, 0.153, 0.304, 0.336, 0.344, 0.327, 0.285, 0.277, 0.242, 0.205, &
                     0.1, 0.066, 0.043, 0.042, 0.079, 0.144, 0.168, 0.155, 0.054, 0.029, &
                     0.013, 0.003, 0.02, 0.0, 0.012, 0.003]

    veg(3)%vname = "leaf" ! from usgs Aspen_Leaf-A DW92-2           BECKa AREF
    veg(3)%lambda = [real(ireals) :: &
                     0.353, 0.499, 0.519, 0.531, 0.553, 0.597, 0.688, 0.694, 0.7, 0.708, &
                     0.728, 0.734, 0.74, 0.751, 0.769, 0.851, 0.949, 1.084, 1.198, 1.244, &
                     1.303, 1.333, 1.378, 1.384, 1.408, 1.423, 1.448, 1.468, 1.534, 1.592, &
                     1.647, 1.7, 1.772, 1.835, 1.855, 1.865, 1.895, 1.905, 1.945, 2.155, &
                     2.235, 2.285, 2.496, 2.56, 2.592] * 1e-6_ireals
    veg(3)%albedo = [real(ireals) :: &
                     0.032, 0.039, 0.053, 0.076, 0.087, 0.053, 0.037, 0.045, 0.073, 0.135, &
                     0.332, 0.379, 0.411, 0.442, 0.458, 0.46, 0.444, 0.438, 0.398, 0.403, &
                     0.398, 0.375, 0.322, 0.305, 0.181, 0.141, 0.124, 0.126, 0.199, 0.247, &
                     0.268, 0.263, 0.233, 0.233, 0.217, 0.192, 0.058, 0.038, 0.029, 0.111, &
                     0.127, 0.093, 0.031, 0.025, 0.031]
  end subroutine

  subroutine init_vegetation_types_from_usgs(comm, ierr, verbose, veg_database_file)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: veg_database_file
    character(len=default_str_len) :: db
    integer(mpiint) :: myid
    ierr = 0

    if (allocated(veg)) return

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      db = get_arg("usgs_material_albedi.txt", veg_database_file)
      call read_usgs_db_file(db, veg, ierr, verbose=verbose); call CHKERR(ierr)
    end if

    call bcast_veg_struct(comm, veg, ierr); call CHKERR(ierr)
  end subroutine

  subroutine bcast_veg_struct(comm, veg, ierr)
    integer(mpiint), intent(in) :: comm
    type(veg_type), intent(inout), allocatable :: veg(:)
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer(iintegers) :: Nmat, i

    ierr = 0

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      if (.not. allocated(veg)) call CHKERR(1_mpiint, 'I expected that rank 0 has database allocated but it isnt')
      Nmat = size(veg)
    end if
    call imp_bcast(comm, Nmat, ierr); call CHKERR(ierr)
    if (.not. allocated(veg)) then
      allocate (veg(Nmat))
    end if

    do i = 1, Nmat
      call imp_bcast(comm, veg(i)%vname, ierr); call CHKERR(ierr)
      call imp_bcast(comm, veg(i)%lambda, ierr); call CHKERR(ierr)
      call imp_bcast(comm, veg(i)%albedo, ierr); call CHKERR(ierr)
    end do
  end subroutine

  subroutine read_usgs_db_file(db, veg, ierr, verbose)
    character(len=*), intent(in) :: db
    type(veg_type), intent(out), allocatable :: veg(:)
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: verbose

    integer :: funit, io
    logical :: file_exists, lverbose
    integer(iintegers) :: nlines, Nmat, i
    character(len=2**14) :: line_str
    !real(ireals), allocatable :: line(:)

    lverbose = get_arg(.false., verbose)
    ierr = 0
    inquire (file=db, exist=file_exists)
    if (.not. file_exists) then
      ierr = 1
      call CHKERR(ierr, 'File '//trim(db)//' does not exist!')
    end if

    open (newunit=funit, file=db)
    read (funit, *) ! skip first line

    ! first count no of lines
    nlines = 0
    do
      read (funit, *, iostat=io)
      if (io .lt. 0) exit ! end of file
      nlines = nlines + 1
    end do
    Nmat = nlines / 3
    if (lverbose) print *, 'Found '//toStr(nlines)//' in usgs database, i.e. '//toStr(Nmat)//' materials.'
    i = modulo(Nlines, 3_iintegers)
    call CHKERR(int(i, mpiint), 'lines in usgs database has to be divideable by 3 !'//&
      & ' However, we have '//toStr(Nmat)//' % 3 = '//toStr(i))

    ! go through it a second time and allocate space
    rewind (funit)
    read (funit, *) ! skip first line

    allocate (veg(Nmat))

    do i = 1, Nmat
      read (funit, '(A)') veg(i)%vname

      read (funit, '(A)', iostat=io) line_str
      call split(line_str, veg(i)%lambda, ' ', ierr); call CHKERR(ierr)

      read (funit, '(A)', iostat=io) line_str
      call split(line_str, veg(i)%albedo, ' ', ierr); call CHKERR(ierr)

      if (lverbose) print *, 'Add '//trim(veg(i)%vname)//' to materials database'
    end do

    close (funit)
  end subroutine

  function get_veg_type_id(vegname) result(id)
    character(len=*), intent(in) :: vegname
    integer(iintegers) :: id
    do id = 1, size(veg)
      if (trim(veg(id)%vname) == trim(vegname)) return
    end do
    call CHKERR(-1_mpiint, 'could not find veg type for veg_name ('//vegname//')')
  end function

  function get_albedo_for_range(veg_id, lambda_min, lambda_max) result(albedo)
    integer(iintegers), intent(in) :: veg_id
    real(ireals), intent(in) :: lambda_min, lambda_max
    real(ireals) :: albedo, albedo1, albedo2

    integer(iintegers) :: Nsample, i
    real(ireals) :: lstart, lend, lambda

    if (lambda_min .gt. lambda_max) &
      & call CHKERR(1_mpiint, 'lambda start has to be smaller than lambda max')
    if (veg_id .lt. 1 .or. veg_id .gt. size(veg)) &
      & call CHKERR(1_mpiint, 'bad veg_id ('//toStr(veg_id)//')')
    lstart = find_real_location(veg(veg_id)%lambda, lambda_min)
    lend = find_real_location(veg(veg_id)%lambda, lambda_max)
    ! ok so we have to integrate over wavelength from [lstart to lend]
    ! which could hold a couple of spectral albedo sample points.
    ! We could try to integrate it thoroughly here
    ! but this will not give the right results anyway because we convolve it with RT results anyway.
    ! So for now we just sample the wavelength spectra in equidistant steps depending on the number of
    ! sample points that we have in the input data because this is the easiest to code and should be fairly robust.
    Nsample = 1_iintegers + int(ceiling(lend - lstart), iintegers)
    albedo = 0
    do i = 1, Nsample
      lambda = linspace(i, [lstart, lend], Nsample)
      albedo = albedo + interp_1d(lambda, veg(veg_id)%albedo)
    end do
    albedo = albedo / real(Nsample, ireals)
    albedo1 = interp_1d(lstart, veg(veg_id)%albedo)
    albedo2 = interp_1d(lend, veg(veg_id)%albedo)
  end function
end module
