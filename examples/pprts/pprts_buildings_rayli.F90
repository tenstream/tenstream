module m_pprts_buildings
  use m_data_parameters, only : init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & zero, pi, i1, i2, default_str_len

  use m_helper_functions, only : CHKERR, spherical_2_cartesian

  use m_pprts, only : init_pprts, &
    & set_optical_properties, solve_pprts, &
    & pprts_get_result, set_angles

  use m_pprts_base, only: t_solver, &
    & allocate_pprts_solver_from_commandline, destroy_pprts, &
    & t_solver_1_2, t_solver_3_6, t_solver_3_10, t_solver_3_16, &
    & t_solver_8_10, t_solver_8_16, t_solver_8_18

  use m_tenstream_options, only: read_commandline_options

  use m_buildings, only: t_pprts_buildings, &
    & faceidx_by_cell_plus_offset, &
    & PPRTS_TOP_FACE, PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, PPRTS_RIGHT_FACE, &
    & PPRTS_REAR_FACE, PPRTS_FRONT_FACE, &
    & init_buildings

  implicit none

contains
  subroutine pprts_buildings(comm, Nx, Ny, Nlay, dx, dy, dz, phi0, theta0, albedo)
    integer(iintegers), intent(in) :: Nx,Ny,Nlay
    real(ireals), intent(in) :: dx, dy, dz
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo
    integer(mpiint), intent(in) :: comm
    real(ireals),parameter :: incSolar = 1
    real(ireals) :: dz1d(Nlay)

    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver
    type(t_pprts_buildings), allocatable :: buildings

    integer(iintegers) :: i, box_k, box_i, box_j
    integer(mpiint) :: ierr

    dz1d = dz

    call init_mpi_data_parameters(comm)
    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(comm, Nlay, Nx, Ny, dx,dy, sundir, solver, dz1d)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = .1_ireals/(dz*Nlay)
    ksca = 1e-3_ireals/(dz*Nlay)
    g    = zero

    call init_buildings(buildings, &
      & [integer(iintegers) :: 6, solver%C_one%zm, solver%C_one%xm,  solver%C_one%ym], &
      & ierr); call CHKERR(ierr)
    allocate(buildings%albedo(6), buildings%iface(6))

    box_k = int((1+Nlay )/ 2.)
    box_i = int((1+Nx)/ 2.)
    box_j = int((1+Ny)/ 2.)

    do i=1,6
      buildings%iface(i) = faceidx_by_cell_plus_offset( &
        & buildings%da_offsets, box_k, box_i, box_j, i )
      buildings%albedo(i) = .2_ireals + i/10._ireals
    enddo

    !buildings%iface(5) = buildings%iface(2)
    !buildings%iface(6) = buildings%iface(3)


    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call set_angles(solver, sundir)

    call solve_pprts(solver, incSolar, opt_buildings=buildings)

    call pprts_get_result(solver, fdn,fup,fdiv,fdir)

    call destroy_pprts(solver, .True.)
  end subroutine

end module


program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_pprts_buildings
  use mpi, only : MPI_COMM_WORLD
  implicit none

  integer(iintegers) :: Nx,Ny,Nlay
  real(ireals) :: dx, dy, dz
  real(ireals) :: phi0, theta0
  real(ireals) :: Ag

  character(len=10*default_str_len) :: default_options
  logical :: lflg, lverbose
  integer(mpiint) :: ierr

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)

  Nx=5; Ny=5; Nlay=3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nlay, lflg, ierr)

  dx = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr)
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr)
  dz = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr)

  phi0 = 180._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr)
  theta0 = 0._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr)

  Ag=0.1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr)

  lverbose = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-verbose', &
    lverbose, lflg, ierr) ; call CHKERR(ierr)


  default_options=''
  default_options=trim(default_options)//' -pprts_use_rayli'
  default_options=trim(default_options)//' -rayli_diff_flx_origin 0,0,-inf'
  default_options=trim(default_options)//' -rayli_cyclic_bc'
  default_options=trim(default_options)//' -show_rayli_dm3d hdf5:dm.h5'

  default_options=trim(default_options)//' -rayli_snapshot'
  default_options=trim(default_options)//' -rayli_snap_photons 1e6'
  default_options=trim(default_options)//' -rayli_snap_Nx 256'
  default_options=trim(default_options)//' -rayli_snap_Ny 256'
  default_options=trim(default_options)//' -visit_image_zoom .75'
  default_options=trim(default_options)//' -visit_parallel_scale 291.5'
  default_options=trim(default_options)//' -visit_focus 300,300,0'
  default_options=trim(default_options)//' -visit_view_normal -0.2811249083446944,-0.7353472951470268,0.6166304739697339'
  default_options=trim(default_options)//' -visit_view_up 0.1878717450780742,0.5879401184069877,0.7867849925925738'

  if(lverbose) print *,'Adding default Petsc Options:', trim(default_options)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, trim(default_options), ierr); call CHKERR(ierr)

  call pprts_buildings(mpi_comm_world, Nx, Ny, Nlay, dx, dy, dz, phi0, theta0, Ag)

  call mpi_finalize(ierr)
end program
