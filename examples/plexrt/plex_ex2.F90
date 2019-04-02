module m_mpi_plex_ex2

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg

use m_data_parameters, only : ireals, iintegers, mpiint, &
  default_str_len, &
  i0, i1, i2, i3, i4, i5,  &
  zero, one,       &
  init_mpi_data_parameters

use m_icon_plex_utils, only: gen_2d_plex_from_icongridfile, icon_hdcp2_default_hhl, &
  dump_ownership, dmplex_2D_to_3D, icon_ncvec_to_plex

use m_plex_grid, only: t_plexgrid, setup_plexgrid, get_normal_of_first_toa_face

use m_plex_rt, only: compute_face_geometry, allocate_plexrt_solver_from_commandline, &
  t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
  destroy_plexrt_solver


use m_netcdfio, only : ncload


implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex2(comm, gridfile, icondatafile)
      MPI_Comm, intent(in) :: comm
      character(len=default_str_len), intent(in) :: gridfile, icondatafile

      integer(mpiint) :: myid, numnodes, ierr
      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: migration_sf
      AO, allocatable :: cell_ao_2d
      type(t_plexgrid), allocatable :: plex
      type(tVec), allocatable :: lwcvec, iwcvec
      real(ireals), parameter :: Ag=.15

      integer(iintegers) :: Nlev
      integer(iintegers), allocatable :: zindex(:)

      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system

      class(t_plex_solver), allocatable :: solver

      call init_mpi_data_parameters(comm)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call gen_2d_plex_from_icongridfile(comm, gridfile, dm2d, dm2d_dist, &
        migration_sf, cell_ao_2d)
      Nlev = size(icon_hdcp2_default_hhl, kind=iintegers)
      call dmplex_2D_to_3D(dm2d_dist, Nlev, icon_hdcp2_default_hhl, dm3d, zindex)

      call dump_ownership(dm3d, '-dump_ownership', '-show_plex')
      call setup_plexgrid(dm3d, Nlev-1, zindex, plex)

      call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'clw', lwcvec, dm3d=dm3d)
      call PetscObjectViewFromOptions(lwcvec, PETSC_NULL_VEC, '-show_lwc', ierr); call CHKERR(ierr)

      call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'cli', iwcvec, dm3d=dm3d)
      call PetscObjectViewFromOptions(iwcvec, PETSC_NULL_VEC, '-show_iwc', ierr); call CHKERR(ierr)

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

      call allocate_plexrt_solver_from_commandline(solver, '5_8')
      call init_plex_rt_solver(plex, solver)

      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal
      !sundir = -[0.677688, 0.0758756, 0.731425]
      !sundir = -[0.826811, 0.0269913, 0.561832]
      !sundir = -[0.775165, 0.0535335, 0.629487] ! sza 10deg
      !sundir = -[0.72362, 0.0793532, 0.685621]
      !sundir = -[0.717145, 0.0262298, 0.696431] !further left
      sundir = -[0.717607, -0.690197, -0.0930992]
      sundir = -[0.88119, -0.0874145, 0.46461]
      sundir = sundir/norm2(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      call set_plex_rt_optprop(solver, vlwc=lwcvec, viwc=iwcvec)

      if(.not.allocated(solver%albedo)) then
        allocate(solver%albedo)
        call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
      endif
      call VecSet(solver%albedo, Ag, ierr); call CHKERR(ierr)

      call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., sundir=sundir)
      call destroy_plexrt_solver(solver, lfinalizepetsc=.False.)
    end subroutine

  end module

  program main
    use m_mpi_plex_ex2
    implicit none

    character(len=default_str_len) :: gridfile, icondatafile, outfile
    logical :: lflg
    integer(mpiint) :: ierr
    character(len=10*default_str_len) :: default_options
    !character(len=*),parameter :: ex_out='plex_ex_dom1_out.h5'
    !character(len=*),parameter :: ex_out='plex_test_out.h5'
    !character(len=*),parameter :: lwcfile='lwc_ex_24_3.nc'
    !character(len=*),parameter :: lwcfile='lwc_ex_dom1.nc'

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)
    call read_commandline_options(PETSC_COMM_WORLD)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid', gridfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a grid filename... please call with -grid <fname_of_icon_gridfile.nc>'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-data', icondatafile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a icondata filename... please call with -data <fname_of_icondatafile.nc>'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

    default_options=''
    default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
    default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_iwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex2(PETSC_COMM_WORLD, gridfile, icondatafile)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
