module m_mpi_plex_ex1

#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, spherical_2_cartesian
  use m_data_parameters, only : ireals, iintegers, mpiint, &
                 default_str_len, i0, i1, i2, i3, i4, i5,  &
                 zero, one, init_mpi_data_parameters

  use m_plex_grid, only: t_plexgrid, load_plex_from_file, &
                       compute_face_geometry, print_dmplex,   &
                       setup_edir_dmplex, setup_abso_dmplex

  use m_plex_rt, only: get_normal_of_first_toa_face, &
   t_plex_solver, init_plex_rt_solver, run_plex_rt_solver

  implicit none

  logical, parameter :: ldebug=.True.

  integer(mpiint) :: ierr

  contains

    subroutine plex_ex1(plex)
      type(t_plexgrid) :: plex

      real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

      type(t_plex_solver), allocatable :: solver

      call init_plex_rt_solver(plex, solver)

      sundir = get_normal_of_first_TOA_face(plex)
      print *,'get_normal_of_first_TOA_face', sundir
      !sundir = [-0.71184089224108049, -3.7794622710146053E-002, -0.70132311428301020] ! original zenith 0
      !sundir = [-0.71184089224108049, -3.7794622710146053E-002, -0.60132311428301020] ! zenith 10 azi 1
      !sundir = [-0.51184089224108049, +0.37794622710146053, -0.00132311428301020] ! zenith 63 azi 17

      !sundir = [-one/100,-one/10,-one]
      !sundir = [one,one,-one]
      !sundir = spherical_2_cartesian(180*one,60*one)
      !sundir = spherical_2_cartesian(90*one,60*one)

      sundir = sundir / norm(sundir)
      print *,'sundir',sundir

      call print_dmplex(plex%comm, plex%edir_dm)

      call run_plex_rt_solver(solver, sundir)

      !call PetscObjectViewFromOptions(edir, PETSC_NULL_VEC, '-show_edir', ierr); call CHKERR(ierr)

    end subroutine
end module

program main
  use m_mpi_plex_ex1
  implicit none

  character(len=default_str_len) :: gridfile
  logical :: lflg

  type(t_plexgrid) :: plex


  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
  call init_mpi_data_parameters(PETSC_COMM_WORLD)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-plex', gridfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a plex filename... please call with -plex <fname_of_plexfile.h5>'

  call load_plex_from_file(PETSC_COMM_WORLD, gridfile, plex)

  call plex_ex1(plex)

  call PetscFinalize(ierr)
end program
