program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : &
    & iintegers, mpiint, ireals, default_str_len, pi, &
    & init_mpi_data_parameters, finalize_mpi
  use m_helper_functions, only: CHKERR
  use m_tenstream_options, only: read_commandline_options
  use m_buildings, only: t_pprts_buildings
  use m_examples_pprts_buildings, only: ex_pprts_buildings
  use mpi, only : MPI_COMM_WORLD
  implicit none

  character(len=default_str_len) :: outfile
  integer(iintegers) :: Nx, Ny, Nlay, icollapse
  integer(iintegers) :: glob_box_i, glob_box_j, box_k
  integer(iintegers) :: box_Ni, box_Nj, box_Nk
  real(ireals) :: box_albedo, box_planck
  real(ireals) :: dx, dy, dz
  real(ireals) :: S0, phi0, theta0
  real(ireals) :: Ag, dtau, w0
  real(ireals),allocatable,dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays on rank 0
  type(t_pprts_buildings), allocatable :: buildings

  character(len=10*default_str_len) :: rayli_options
  logical :: lflg, lverbose, lrayli_opts, lsolar, lthermal
  integer(mpiint) :: ierr

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  lsolar = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-solar', &
    lsolar, lflg, ierr) ; call CHKERR(ierr)

  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-thermal', &
    lthermal, lflg, ierr) ; call CHKERR(ierr)

  Nx=5; Ny=5; Nlay=3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nlay, lflg, ierr); call CHKERR(ierr)

  icollapse=1
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-icollapse", icollapse, lflg, ierr); call CHKERR(ierr)

  box_k = Nlay-(icollapse-1) ! touching the surface
  box_k = box_k-1 ! one above, i.e. hovering
  glob_box_i = int(real(Nx+1)/2.)
  glob_box_j = int(real(Nx+1)/2.)
  box_Ni = 1
  box_Nj = 1
  box_Nk = 1
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Bx", glob_box_i, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-By", glob_box_j, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Bz", box_k     , lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-BNx",box_Ni    , lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-BNy",box_Nj    , lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-BNz",box_Nk    , lflg, ierr); call CHKERR(ierr)

  box_albedo = .1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-BAg", box_albedo, lflg, ierr); call CHKERR(ierr)
  box_planck = 100._ireals / pi
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Bplanck", box_planck, lflg, ierr); call CHKERR(ierr)

  dx = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)
  dz = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr); call CHKERR(ierr)

  S0 = 1._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-S0", S0, lflg, ierr); call CHKERR(ierr)

  phi0 = 180._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 0._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  Ag=0.1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr); call CHKERR(ierr)
  dtau=1._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dtau", dtau, lflg, ierr); call CHKERR(ierr)
  w0=0.5_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-w0", w0, lflg, ierr); call CHKERR(ierr)

  lverbose = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-verbose', &
    lverbose, lflg, ierr) ; call CHKERR(ierr)

  lrayli_opts = .False.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-rayli_opts', &
    lrayli_opts, lflg, ierr) ; call CHKERR(ierr)

  if(lrayli_opts) then
    rayli_options=''
    rayli_options=trim(rayli_options)//' -pprts_use_rayli'
    rayli_options=trim(rayli_options)//' -rayli_diff_flx_origin 0,0,-inf'
    rayli_options=trim(rayli_options)//' -rayli_cyclic_bc'
    rayli_options=trim(rayli_options)//' -show_rayli_dm3d hdf5:dm.h5'

    rayli_options=trim(rayli_options)//' -rayli_snapshot'
    rayli_options=trim(rayli_options)//' -rayli_snap_Nx 256'
    rayli_options=trim(rayli_options)//' -rayli_snap_Ny 256'
    rayli_options=trim(rayli_options)//' -visit_image_zoom 5.75'
    rayli_options=trim(rayli_options)//' -visit_parallel_scale 400'
    rayli_options=trim(rayli_options)//' -visit_focus 250,250,0'
    rayli_options=trim(rayli_options)//' -visit_view_normal 0.4,-0.4,0.8'
    rayli_options=trim(rayli_options)//' -visit_view_up -0.5,0.5,0.6'

    if(lverbose) print *,'Adding rayli Petsc Options:', trim(rayli_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, trim(rayli_options), ierr); call CHKERR(ierr)
  endif

  if(lsolar) then
    call ex_pprts_buildings(mpi_comm_world, lverbose, &
      & .False., lsolar, Nx, Ny, Nlay, icollapse, &
      & glob_box_i, glob_box_j, box_k,            &
      & box_Ni, box_Nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings,                                &
      & outfile=outfile)
  endif
  if(lthermal) then
    call ex_pprts_buildings(mpi_comm_world, lverbose, &
      & lthermal, .False., Nx, Ny, Nlay, icollapse, &
      & glob_box_i, glob_box_j, box_k,            &
      & box_Ni, box_Nj, box_Nk,                   &
      & box_albedo, box_planck,                   &
      & dx, dy, dz,                               &
      & S0, phi0, theta0,                         &
      & Ag, dtau, w0,                             &
      & gedir, gedn, geup, gabso,                 &
      & buildings,                                &
      & outfile=outfile)
  endif

  call finalize_mpi(ierr, .True., .True.)
end program
