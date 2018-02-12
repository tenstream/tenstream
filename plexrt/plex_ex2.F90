module m_mpi_plex_ex2

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
  distribute_plexgrid_dm, ncvar2d_to_globalvec

use m_plex_rt, only: get_normal_of_first_toa_face, compute_face_geometry, &
  t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop

use m_netcdfio, only : ncload


implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex2(comm, gridfile, hhlfile, icondatafile)
      MPI_Comm, intent(in) :: comm
      character(len=default_str_len), intent(in) :: gridfile, hhlfile, icondatafile

      type(t_icongrid),allocatable :: icongrid, local_icongrid

      integer(mpiint) :: myid, numnodes, ierr
      type(t_plexgrid), allocatable :: plexgrid
      type(tVec) :: lwcvec, iwcvec
      AO :: cell_ao

      integer(iintegers) :: Nz
      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system
      real(ireals),allocatable :: hhl(:)
      character(len=default_str_len) :: ncgroups(2)

      type(t_plex_solver), allocatable :: solver

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
          call read_icon_grid_file(gridfile, icongrid)
      endif
      call bcast_icongrid(comm, icongrid)

      call distribute_icon_grid(comm, icongrid, local_icongrid, cell_ao)
      deallocate(icongrid)

      if(myid.eq.0) then
        ncgroups(1) = trim(hhlfile)
        ncgroups(2) = 'height_level'; call ncload(ncgroups, hhl, ierr)
      endif
      call imp_bcast(comm, hhl, 0)
      Nz = size(hhl)-1

      call create_plex_from_icongrid(comm, Nz, hhl, cell_ao, local_icongrid, plexgrid)
      deallocate(local_icongrid)

      call ncvar2d_to_globalvec(plexgrid, icondatafile, 'clw', lwcvec)
      call PetscObjectViewFromOptions(lwcvec, PETSC_NULL_VEC, '-show_lwc', ierr); call CHKERR(ierr)

      call ncvar2d_to_globalvec(plexgrid, icondatafile, 'cli', iwcvec)
      call PetscObjectViewFromOptions(iwcvec, PETSC_NULL_VEC, '-show_iwc', ierr); call CHKERR(ierr)

      call init_plex_rt_solver(plexgrid, solver)

      call compute_face_geometry(solver%plex, solver%plex%geom_dm)
      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal
      !sundir = -[0.677688, 0.0758756, 0.731425]
      !sundir = -[0.826811, 0.0269913, 0.561832]
      !sundir = -[0.775165, 0.0535335, 0.629487] ! sza 10deg
      !sundir = -[0.72362, 0.0793532, 0.685621]
      !sundir = -[0.717145, 0.0262298, 0.696431] !further left
      sundir = -[0.717607, -0.690197, -0.0930992]
      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      call set_plex_rt_optprop(solver, vlwc=lwcvec, viwc=iwcvec)

      call run_plex_rt_solver(solver, sundir)
    end subroutine

  end module

  program main
    use m_mpi_plex_ex2
    implicit none

    character(len=default_str_len) :: gridfile, hhlfile, icondatafile, outfile
    logical :: lflg
    integer(mpiint) :: ierr
    character(len=10*default_str_len) :: default_options
    !character(len=*),parameter :: ex_out='plex_ex_dom1_out.h5'
    !character(len=*),parameter :: ex_out='plex_test_out.h5'
    !character(len=*),parameter :: lwcfile='lwc_ex_24_3.nc'
    !character(len=*),parameter :: lwcfile='lwc_ex_dom1.nc'

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)
    call read_commandline_options()

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid', gridfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a grid filename... please call with -grid <fname_of_icon_gridfile.nc>'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-hhl', hhlfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a hhl filename... please call with -hhl <fname_of_hhlfile.nc>'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-icondata', icondatafile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a icondata filename... please call with -icondata <fname_of_icondatafile.nc>'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

    default_options=''
    default_options=trim(default_options)//' -show_plex_dump hdf5:'//trim(outfile)
    default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex2(PETSC_COMM_WORLD, gridfile, hhlfile, icondatafile)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
