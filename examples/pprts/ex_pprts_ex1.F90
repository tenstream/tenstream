program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, pi
  use m_helper_functions, only : CHKERR
  use m_examples_pprts_ex1, only: pprts_ex1

  implicit none

  integer(iintegers) :: Nx, Ny, Nlay
  real(ireals) :: dx, dy, dz
  real(ireals) :: phi0, theta0
  real(ireals) :: albedo
  real(ireals) :: incSolar
  real(ireals) :: Bplck, Bplck_srfc
  logical :: lthermal, lsolar

  real(ireals) :: dtau_clearsky, w0_clearsky, g_clearsky
  integer(iintegers) :: cld_layer_idx(2)
  real(ireals) :: dtau_cloud, w0_cloud, g_cloud

  real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

  integer(iintegers) :: k, Ncld_idx
  logical :: lflg, lverbose
  integer(mpiint) :: myid, ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call init_mpi_data_parameters(mpi_comm_world)

  Nx=3; Ny=3; Nlay=10
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nlay, lflg, ierr); call CHKERR(ierr)

  dx=100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dy = dx
  dz = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr); call CHKERR(ierr)

  lthermal=.False.
  lsolar  =.True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg, ierr); call CHKERR(ierr)

  phi0 = 180
  theta0 = 0
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  albedo = 0.1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", albedo, lflg, ierr); call CHKERR(ierr)

  incSolar = 1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-S0", incSolar, lflg, ierr); call CHKERR(ierr)

  Bplck = 100._ireals / pi
  Bplck_srfc = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-B0", Bplck, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Bs", Bplck_srfc, lflg, ierr); call CHKERR(ierr)


  dtau_clearsky=1; w0_clearsky=.5; g_clearsky=0
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dtau_clr", dtau_clearsky, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-w0_clr", w0_clearsky, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-g_clr", g_clearsky, lflg, ierr); call CHKERR(ierr)

  cld_layer_idx = [Nlay/2+1, Nlay/2+1]
  Ncld_idx = 2
  call PetscOptionsGetIntArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-cld_idx", &
    & cld_layer_idx, Ncld_idx, lflg, ierr); call CHKERR(ierr)

  dtau_cloud=1; w0_cloud=.99; g_cloud=.9
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dtau_cld", dtau_cloud, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-w0_cld", w0_cloud, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-g_cld", g_cloud, lflg, ierr); call CHKERR(ierr)

  lverbose=.True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-v", lverbose, lflg, ierr); call CHKERR(ierr)


  call pprts_ex1( &
      & mpi_comm_world, &
      & lthermal, &
      & lsolar, &
      & Nx, Ny, Nlay, &
      & dx, dy, &
      & phi0, theta0, &
      & albedo, dz, &
      & incSolar, &
      & Bplck, &
      & Bplck_srfc, &
      & dtau_clearsky, w0_clearsky, g_clearsky, &
      & cld_layer_idx, &
      & dtau_cloud, w0_cloud, g_cloud, &
      & fdir,fdn,fup,fdiv, &
      & lverbose )

  if(myid.eq.0 .and. lverbose) then
    print *,''
    print *,'Fluxes in first column on rank 0'
    do k = lbound(fdiv,1), ubound(fdiv,1)
      print *, 'k=', k, &
        & 'fdir:', fdir(k,1,1), &
        & 'fdn:',  fdn (k,1,1), &
        & 'fup:',  fup (k,1,1), &
        & 'fdiv:', fdiv(k,1,1)
    enddo
    k = ubound(fdir,1)
    print *,'k=', k, &
      & 'fdir:', fdir(k,1,1), &
      & 'fdn:',  fdn (k,1,1), &
      & 'fup:',  fup (k,1,1)
  endif

  if(myid.eq.0) then
    print *,''
    print *,''
    print *,'Call this example e.g. with options: -show_edir hdf5:edir.h5'
    print *,'and plot results with python:'
    print *,'import h5py as H; h=H.File("edir.h5","r"); edir = h["edir0"][:]'
    print *,'imshow(edir[0,:,:,0].T,interpolation="nearest");' ! has dimension nyp,nxp,nzp,8streams
    print *,'colorbar(); savefig("edir_x0.pdf")'
  endif
  call mpi_finalize(ierr)
end program
