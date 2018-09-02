module m_mpi_plex_ex4

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg, meanvec, reverse

use m_data_parameters, only : ireals, iintegers, mpiint, &
  default_str_len, &
  i0, i1, i2, i3, i4, i5,  &
  zero, one,       &
  init_mpi_data_parameters

use m_icon_plex_utils, only: gen_2d_plex_from_icongridfile, icon_hdcp2_default_hhl, &
  dump_ownership, dmplex_2D_to_3D, icon_ncvec_to_plex

use m_plex_grid, only: t_plexgrid, setup_plexgrid, get_normal_of_first_toa_face

use m_plex_rt, only: compute_face_geometry, &
  t_plex_solver, init_plex_rt_solver

use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, hydrostat_plev, print_tenstr_atm

use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

use m_netcdfio, only : ncload


implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex4(comm, gridfile, icondatafile)
      MPI_Comm, intent(in) :: comm
      character(len=default_str_len), intent(in) :: gridfile, icondatafile

      integer(mpiint) :: myid, numnodes, ierr
      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: migration_sf
      type(AO), allocatable :: cell_ao_2d
      type(t_plexgrid), allocatable :: plex
      !type(tVec), allocatable :: lwcvec2d, iwcvec2d
      real(ireals), allocatable :: col_plev(:,:), col_tlev(:,:), dp(:)
      real(ireals), allocatable :: col_lwc(:,:), col_reliq(:,:)
      real(ireals), allocatable :: col_iwc(:,:), col_reice(:,:)
      real(ireals), parameter :: Ag=.15, lapse_rate=6.5e-3, Tsrfc=288.3, minTemp=Tsrfc-50._ireals

      character(len=default_str_len),parameter :: atm_filename='afglus_100m.dat'
      type(t_tenstr_atm) :: atm
      integer(iintegers), allocatable :: zindex(:)
      integer(iintegers) :: k, fStart, fEnd, Ncol, Nlev, Nlay
      real(ireals), allocatable, dimension(:,:) :: edir,edn,eup,abso

      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system

      class(t_plex_solver), allocatable :: solver

      call init_mpi_data_parameters(comm)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call gen_2d_plex_from_icongridfile(comm, gridfile, dm2d, dm2d_dist, &
        migration_sf, cell_ao_2d)

      if(myid.eq.0) print *,'Could read data from icondatafile', trim(icondatafile)
      !call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'clw', lwcvec2d)
      !call PetscObjectViewFromOptions(lwcvec2d, PETSC_NULL_VEC, '-show_lwc', ierr); call CHKERR(ierr)

      !call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'cli', iwcvec2d)
      !call PetscObjectViewFromOptions(iwcvec2d, PETSC_NULL_VEC, '-show_iwc', ierr); call CHKERR(ierr)

      call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
      Ncol = fEnd - fStart
      Nlev = size(icon_hdcp2_default_hhl); Nlay = Nlev-1

      if(myid.eq.0) print *,'Dynamics Grid has Size Nlev, Ncol:', Nlev, Ncol

      ! prepare atmosphere
      allocate(col_tlev(Nlev, Ncol))
      allocate(col_plev(Nlev, Ncol))

      do k = 1, Nlev
        col_tlev(k, i1) = max(minTemp, Tsrfc - (lapse_rate * icon_hdcp2_default_hhl(Nlev-(k-1))))
      enddo

      allocate(dp(Nlay))
      call hydrostat_plev(1013._ireals, meanvec(col_tlev(:,i1)), reverse(icon_hdcp2_default_hhl), &
        col_plev(:,i1), dp)
      deallocate(dp)

      do k = 2, Ncol
        col_plev(:, k) =  col_plev(:, i1)
        col_tlev(:, k) =  col_tlev(:, i1)
      enddo

      allocate(col_lwc(Nlay, Ncol), col_reliq(Nlay, Ncol), &
               col_iwc(Nlay, Ncol), col_reice(Nlay, Ncol), source=zero)

      call setup_tenstr_atm(comm, .False., atm_filename, &
        col_plev, col_tlev, atm, &
        d_lwc=col_lwc, d_reliq=col_reliq, &
        d_iwc=col_iwc, d_reice=col_reice)

      !call print_tenstr_atm(atm)

      call dmplex_2D_to_3D(dm2d_dist, reverse(atm%zt(:, i1)), dm3d, zindex)

      call dump_ownership(dm3d, '-dump_ownership', '-show_plex')
      call setup_plexgrid(dm3d, zindex, reverse(atm%zt(:, i1)), plex)

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

      call init_plex_rt_solver(plex, solver)

      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal
      sundir = -[0.717607, -0.690197, -0.0930992]
      sundir = -[0.88119, -0.0874145, 0.46461]
      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      call plexrt_rrtmg(solver, atm, sundir, &
        albedo_thermal=zero, albedo_solar=Ag, &
        lthermal=.True., lsolar=.False., &
        edir=edir, edn=edn, eup=eup, abso=abso)

      call destroy_plexrt_rrtmg(solver, lfinalizepetsc=.False.)
    end subroutine

  end module

  program main
    use m_mpi_plex_ex4
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

    call plex_ex4(PETSC_COMM_WORLD, gridfile, icondatafile)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
