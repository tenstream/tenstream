module m_mpi_plex_ex2

#include "petsc/finclude/petsc.h"
use petsc
use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, spherical_2_cartesian
use m_data_parameters, only : ireals, iintegers, mpiint, &
  default_str_len, &
  i0, i1, i2, i3, i4, i5,  &
  zero, one,       &
  init_mpi_data_parameters, myid

use m_icon_grid, only: t_icongrid, read_icon_grid_file, &
  decompose_icon_grid, bcast_icongrid, distribute_icon_grid

use m_plex_grid, only: t_plexgrid, create_plex_from_icongrid, &
  setup_edir_dmplex, setup_abso_dmplex, compute_face_geometry, &
  distribute_plexgrid_dm

use m_plex_rt, only: get_normal_of_first_toa_face, create_src_vec, &
  create_edir_mat, solve_plex_rt, compute_edir_absorption


implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex2(comm, gridfile)
      MPI_Comm, intent(in) :: comm
      character(len=default_str_len) :: gridfile

      type(t_icongrid),allocatable :: icongrid, local_icongrid

      integer(mpiint) :: myid, numnodes, ierr
      type(t_plexgrid), allocatable :: plexgrid
      type(tVec) :: b, edir, abso
      type(tMat) :: A

      integer(iintegers), parameter :: Nz=10
      real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
          call read_icon_grid_file(gridfile, icongrid)
      endif
      call bcast_icongrid(comm, icongrid)

      call distribute_icon_grid(comm, icongrid, local_icongrid)
      deallocate(icongrid)

      call create_plex_from_icongrid(comm, Nz, local_icongrid, plexgrid)
      deallocate(local_icongrid)

      call compute_face_geometry(plexgrid) ! setup the geometry info for the plexgrid object (sets up the plexgrid%geom_dm)

      call setup_edir_dmplex(plexgrid, plexgrid%edir_dm)

      sundir = get_normal_of_first_TOA_face(plexgrid) + [0.,0.,-.1]
      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir

      call create_src_vec(plexgrid%edir_dm, b)

      call setup_abso_dmplex(plexgrid, plexgrid%abso_dm)

      call create_edir_mat(plexgrid, sundir, A)

      call solve_plex_rt(plexgrid, b, A, edir)

      call PetscObjectViewFromOptions(edir, PETSC_NULL_VEC, '-show_edir', ierr); call CHKERR(ierr)

      call compute_edir_absorption(plexgrid, edir, sundir, abso)

    end subroutine
  end module

  program main
    use m_mpi_plex_ex2
    implicit none

    character(len=default_str_len) :: gridfile
    logical :: lflg
    integer(mpiint) :: ierr
    character(len=*),parameter :: default_options='-show_plex hdf5:plex_ex2_out.h5 -show_abso hdf5:plex_ex2_out.h5::append'

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)

    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid', gridfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a plex filename... please call with -grid <fname_of_icon_gridfile.nc>'

    call plex_ex2(PETSC_COMM_WORLD, gridfile)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
