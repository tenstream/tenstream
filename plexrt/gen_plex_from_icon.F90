module m_gen_plex_from_icon

#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR

  use m_icon_plex_utils, only: gen_2d_plex_from_icongridfile, icon_hdcp2_default_hhl, &
    dump_ownership, dmplex_2D_to_3D
  use m_plex_grid, only: t_plexgrid, setup_plexgrid, ncvar2d_to_globalvec
  use m_data_parameters, only: ireals, iintegers, mpiint, &
    default_str_len, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one,       &
    init_mpi_data_parameters

  implicit none

  logical, parameter :: ldebug=.True.

  contains

    subroutine gen_distributed_3d_plex()
      logical :: lflg
      character(len=default_str_len) :: gridfile, datafile

      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: distribute_point_sf
      type(AO), allocatable :: cell_ao_2d
      type(t_plexgrid), allocatable :: plex
      integer(mpiint) :: comm, myid, ierr

      type(tVec), allocatable :: clw

      integer(iintegers) :: Nlev
      integer(iintegers), allocatable :: zindex(:)

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
      comm = MPI_COMM_WORLD

      call init_mpi_data_parameters(comm)
      call mpi_comm_rank(comm, myid, ierr)

      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-grid', gridfile, lflg, ierr); call CHKERR(ierr)
      if(.not.lflg) then
        print *,'Please specify a grid file with option:'
        print *,'-grid <path_to_icon_grid_file.nc>'
        call CHKERR(1_mpiint, 'Required Option missing')
      endif

      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-data', datafile, lflg, ierr); call CHKERR(ierr)
      if(.not.lflg) then
        print *,'Please specify a data file with option:'
        print *,'-data <path_to_icon_data_file.nc>'
        call CHKERR(1_mpiint, 'Required Option missing')
      endif

      call gen_2d_plex_from_icongridfile(comm, gridfile, dm2d, dm2d_dist, distribute_point_sf, cell_ao_2d)
      Nlev = size(icon_hdcp2_default_hhl, kind=iintegers)
      call dmplex_2D_to_3D(dm2d_dist, Nlev, icon_hdcp2_default_hhl, dm3d, zindex)

      call dump_ownership(dm3d, '-dump_ownership', '-show_plex')
      call setup_plexgrid(dm3d, Nlev-1, zindex, plex)

      call ncvar2d_to_globalvec(plex, datafile, 'clw', clw, cell_ao_2d=cell_ao_2d)
      call PetscObjectViewFromOptions(clw, PETSC_NULL_VEC, '-show_clw', ierr); call CHKERR(ierr)

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)
      call DMDestroy(dm3d, ierr); call CHKERR(ierr)
      call PetscFinalize(ierr)
      print *,'m_gen_plex_from_icon... done'
    end subroutine

end module

program main
  use m_gen_plex_from_icon
  implicit none

  call gen_distributed_3d_plex()
end program
