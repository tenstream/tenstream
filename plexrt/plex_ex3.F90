module m_mpi_plex_ex3

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg

use m_data_parameters, only : ireals, iintegers, mpiint, &
  default_str_len, &
  i0, i1, i2, i3, i4, i5,  &
  zero, one,       &
  init_mpi_data_parameters

use m_icon_grid, only: t_icongrid, read_icon_grid_file, &
  bcast_icongrid, distribute_icon_grid

use m_plex_grid, only: t_plexgrid, create_plex_from_icongrid, &
  setup_edir_dmplex, setup_abso_dmplex, compute_face_geometry, &
  distribute_plexgrid_dm, ncvar2d_to_globalvec, setup_plexgrid, &
  gen_test_mat

use m_plex_rt, only: get_normal_of_first_toa_face, compute_face_geometry, &
  t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
  destroy_plexrt_solver

use m_netcdfio, only : ncload

use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D, dump_ownership

implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex3(comm, Nx, Ny, Nz, dz)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nx, Ny, Nz
      real(ireals), intent(in) :: dz

      type(tDM) :: dm2d, dm3d
      real(ireals) :: hhl(Nz)

      integer(mpiint) :: myid, numnodes, ierr
      integer(iintegers) :: k
      !type(tVec) :: lwcvec, iwcvec

      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system
      !real(ireals),allocatable :: hhl(:)
      !character(len=default_str_len) :: ncgroups(2)

      type(t_plexgrid), allocatable :: plex
      integer(iintegers), allocatable :: zindex(:)
      type(t_plex_solver), allocatable :: solver

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call create_2d_fish_plex(dm2d, Nx, Ny)

      hhl(1) = zero
      do k=2,Nz
        hhl(k) = hhl(k-1) - dz
      enddo

      call dmplex_2D_to_3D(dm2d, hhl, dm3d, zindex)

      call setup_plexgrid(dm3d, zindex, hhl, plex)
      deallocate(zindex)

      if(ldebug.and.myid.eq.0) print *,'create_plex_from_icongrid :: Setup Connections : show plex'
      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)
      call dump_ownership(plex%dm, '-show_ownership')

      !call create_plex_from_icongrid(comm, Nz, hhl, cell_ao, local_icongrid, plexgrid)
      !deallocate(local_icongrid)

      !call ncvar2d_to_globalvec(plexgrid, icondatafile, 'clw', lwcvec)
      !call PetscObjectViewFromOptions(lwcvec, PETSC_NULL_VEC, '-show_lwc', ierr); call CHKERR(ierr)

      !call ncvar2d_to_globalvec(plexgrid, icondatafile, 'cli', iwcvec)
      !call PetscObjectViewFromOptions(iwcvec, PETSC_NULL_VEC, '-show_iwc', ierr); call CHKERR(ierr)

      call init_plex_rt_solver(plex, solver)

      call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal + [zero, .001_ireals, zero]
      !sundir = -[0.688915, -0.422213, 0.589179]
      !sundir = -[-0.003132140126484546, -0.9198186384401722, 0.3923313167162367]
      !!sundir = -[0.677688, 0.0758756, 0.731425]
      !!sundir = -[0.826811, 0.0269913, 0.561832]
      !!sundir = -[0.775165, 0.0535335, 0.629487] ! sza 10deg
      !!sundir = -[0.72362, 0.0793532, 0.685621]
      !!sundir = -[0.717145, 0.0262298, 0.696431] !further left
      !sundir = -[0.717607, -0.690197, -0.0930992]
      !sundir = -[0.88119, -0.0874145, 0.46461]
      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      !!call set_plex_rt_optprop(solver, vlwc=lwcvec, viwc=iwcvec)
      call set_plex_rt_optprop(solver)

      call run_plex_rt_solver(solver, sundir)
      call destroy_plexrt_solver(solver, lfinalizepetsc=.False.)
    end subroutine

  end module

  program main
    use m_mpi_plex_ex3
    implicit none

    character(len=default_str_len) :: outfile
    logical :: lflg
    integer(mpiint) :: ierr
    character(len=10*default_str_len) :: default_options
    integer(iintegers) :: Nx, Ny, Nz
    real(ireals) :: dz

    !character(len=*),parameter :: ex_out='plex_ex_dom1_out.h5'
    !character(len=*),parameter :: ex_out='plex_test_out.h5'
    !character(len=*),parameter :: lwcfile='lwc_ex_24_3.nc'
    !character(len=*),parameter :: lwcfile='lwc_ex_dom1.nc'

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)
    call read_commandline_options(PETSC_COMM_WORLD)

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) Nx = 2
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) Ny = 3
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) Nz = 2
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) dz = one/Nz

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

    default_options='-polar_coords no'
    default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
    default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_iwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex3(PETSC_COMM_WORLD, Nx, Ny, Nz, dz)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
