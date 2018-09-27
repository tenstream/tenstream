module m_mpi_plex_ex3

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg, reverse

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
  gen_test_mat, get_normal_of_first_toa_face

use m_plex_rt, only: compute_face_geometry, &
  t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
  plexrt_get_result, destroy_plexrt_solver

use m_netcdfio, only : ncload

use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D, dump_ownership

implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex3(comm, Nx, Ny, Nz, dz, Ag)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nx, Ny, Nz
      real(ireals), intent(in) :: dz, Ag

      type(tDM) :: dm2d, dm3d
      real(ireals) :: hhl(Nz)
      real(ireals), allocatable :: hhl_3d(:)

      integer(mpiint) :: myid, numnodes, ierr
      integer(iintegers) :: i,k, v2dStart, v2dEnd
      !type(tVec) :: lwcvec, iwcvec

      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system
      !real(ireals),allocatable :: hhl(:)
      !character(len=default_str_len) :: ncgroups(2)

      type(t_plexgrid), allocatable :: plex
      integer(iintegers), allocatable :: zindex(:)
      class(t_plex_solver), allocatable :: solver

      real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

      logical, parameter :: lthermal = .False.

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call create_2d_fish_plex(Nx, Ny, dm2d)

      hhl(1) = zero
      do k=2,Nz
        hhl(k) = hhl(k-1) + dz
      enddo
      hhl = reverse(hhl)

      call DMPlexGetDepthStratum(dm2d, i0, v2dStart, v2dEnd, ierr); call CHKERR(ierr) ! 2D vertices

      allocate(hhl_3d(Nz*(v2dEnd-v2dStart)))
      do i=v2dStart, v2dEnd-1
          do k=1,Nz
            hhl_3d( (i-v2dStart)*Nz + k) = hhl(k) + real(i, ireals)/Nx/Ny
        enddo
      enddo
      print *,'hhl_3d', hhl_3d

      call dmplex_2D_to_3D(dm2d, Nz, hhl_3d, dm3d, zindex)

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

      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal !+ [zero, -.00001_ireals, zero]

      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      !!call set_plex_rt_optprop(solver, vlwc=lwcvec, viwc=iwcvec)
      call set_plex_rt_optprop(solver)


      if(.not.allocated(solver%albedo)) then
        allocate(solver%albedo)
        call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
      endif
      call VecSet(solver%albedo, Ag, ierr); call CHKERR(ierr)

      if(lthermal) then
        if(.not.allocated(solver%plck)) then
          allocate(solver%plck)
          call DMCreateGlobalVector(solver%plex%cell1_dm, solver%plck, ierr); call CHKERR(ierr)
        endif
        call VecSet(solver%plck, 100._ireals, ierr); call CHKERR(ierr)

        if(.not.allocated(solver%srfc_emission)) then
          allocate(solver%srfc_emission)
          call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%srfc_emission, ierr); call CHKERR(ierr)
        endif
        call VecSet(solver%srfc_emission, 400._ireals/3.1415_ireals, ierr); call CHKERR(ierr)
        call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=sundir)
      else
        call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., sundir=sundir)
      endif

      call plexrt_get_result(solver, edn, eup, abso, redir=edir)
      call destroy_plexrt_solver(solver, lfinalizepetsc=.False.)

      print *,'k                Edir                Edn                          Eup                          abso'
      do k = 1, ubound(abso,1)
        print *,k, edir(k,1), edn(k,1), eup(k,1), abso(k,1)
      enddo
      print *,k, edir(k,1), edn(k,1), eup(k,1)
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
    real(ireals) :: dz, Ag

    !character(len=*),parameter :: ex_out='plex_ex_dom1_out.h5'
    !character(len=*),parameter :: ex_out='plex_test_out.h5'
    !character(len=*),parameter :: lwcfile='lwc_ex_24_3.nc'
    !character(len=*),parameter :: lwcfile='lwc_ex_dom1.nc'

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)
    call read_commandline_options(PETSC_COMM_WORLD)

    Nx = 2
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
    Ny = 3
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
    Nz = 2
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg,ierr) ; call CHKERR(ierr)
    dz = one/Nz
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)
    Ag = .1_ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

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
    default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex3(PETSC_COMM_WORLD, Nx, Ny, Nz, dz, Ag)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
