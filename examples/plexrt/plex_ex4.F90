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
  dump_ownership, dmplex_2D_to_3D, icon_ncvec_to_plex, celldm_veccopy, celldm1_vec_to_nz_ncol, &
  dm2d_vec_to_Nz_Ncol

use m_plex_grid, only: t_plexgrid, setup_plexgrid, get_normal_of_first_toa_face

use m_plex_rt, only: compute_face_geometry, &
  t_plex_solver, init_plex_rt_solver

use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, hydrostat_plev, print_tenstr_atm

use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

use m_netcdfio, only : ncload


implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine plex_ex4(comm, gridfile, icondatafile, Ag, lthermal, lsolar)
      MPI_Comm, intent(in) :: comm
      character(len=default_str_len), intent(in) :: gridfile, icondatafile
      real(ireals), intent(in) :: Ag
      logical, intent(in) :: lthermal, lsolar

      integer(mpiint) :: myid, numnodes, ierr
      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: migration_sf
      type(AO), allocatable :: cell_ao_2d
      type(t_plexgrid), allocatable :: plex
      !type(tVec), allocatable :: lwcvec2d, iwcvec2d
      real(ireals), allocatable :: col_plev(:,:), col_tlev(:,:), dp(:)
      real(ireals), allocatable :: col_lwc(:,:), col_reliq(:,:)
      real(ireals), allocatable :: col_iwc(:,:), col_reice(:,:)
      real(ireals), parameter :: lapse_rate=6.5e-3, Tsrfc=288.3, minTemp=Tsrfc-50._ireals

      character(len=default_str_len),parameter :: atm_filename='afglus_100m.dat'
      type(t_tenstr_atm) :: atm
      integer(iintegers), allocatable :: zindex(:)
      integer(iintegers) :: k, fStart, fEnd, Ncol, Nlev, Nlay
      real(ireals), allocatable, dimension(:,:) :: edir,edn,eup,abso

      real(ireals) :: first_normal(3), sundir(3) ! cartesian direction of sun rays in a global reference system

      class(t_plex_solver), allocatable :: solver

      type(tVec), allocatable :: lwcvec, iwcvec

      call init_mpi_data_parameters(comm)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call gen_2d_plex_from_icongridfile(comm, gridfile, dm2d, dm2d_dist, &
        migration_sf, cell_ao_2d)

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

      call setup_tenstr_atm(comm, .False., atm_filename, col_plev, col_tlev, atm)

      !call print_tenstr_atm(atm)

      ! Setup 3D DMPLEX grid
      call dmplex_2D_to_3D(dm2d_dist, reverse(atm%zt(:, i1)), dm3d, zindex)

      call dump_ownership(dm3d, '-dump_ownership', '-show_plex')
      call setup_plexgrid(dm3d, zindex, reverse(atm%zt(:, i1)), plex)

      !Load Data from iconfile and distribute it
      if(myid.eq.0) print *,'Read data from icondatafile', trim(icondatafile)
      call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'clw', lwcvec)
      call PetscObjectViewFromOptions(lwcvec, PETSC_NULL_VEC, '-show_lwc', ierr); call CHKERR(ierr)

      call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, 'cli', iwcvec)
      call PetscObjectViewFromOptions(iwcvec, PETSC_NULL_VEC, '-show_iwc', ierr); call CHKERR(ierr)

      allocate(col_lwc(Nlay, Ncol), col_reliq(Nlay, Ncol), &
               col_iwc(Nlay, Ncol), col_reice(Nlay, Ncol), source=zero)

      call dm2d_vec_to_Nz_Ncol(dm2d_dist, lwcvec, col_lwc)
      call dm2d_vec_to_Nz_Ncol(dm2d_dist, iwcvec, col_iwc)
      col_lwc = col_lwc * 1e3_ireals ! kg to g
      col_iwc = col_iwc * 1e3_ireals
      col_reliq(:,:) = 2.5_ireals
      col_reice(:,:) = 10._ireals

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

      call setup_tenstr_atm(comm, .False., atm_filename, &
        col_plev, col_tlev, atm, &
        d_lwc=col_lwc, d_reliq=col_reliq, &
        d_iwc=col_iwc, d_reice=col_reice)

      call init_plex_rt_solver(plex, solver)

      first_normal = get_normal_of_first_TOA_face(solver%plex)
      sundir = first_normal + [zero, zero, 1e-2_ireals]
      !sundir = -[0.717607, -0.690197, -0.0930992]
      !sundir = -[0.88119, -0.0874145, 0.46461]
      sundir = sundir/norm(sundir)
      print *,myid,'Initial sundirection = ', sundir, rad2deg(angle_between_two_vec(sundir, first_normal))

      if(lthermal) then
        call plexrt_rrtmg(solver, atm, sundir, &
          albedo_thermal=zero, albedo_solar=Ag, &
          lthermal=.True., lsolar=.False., &
          edir=edir, edn=edn, eup=eup, abso=abso)
      endif

      if(lsolar) then
        call plexrt_rrtmg(solver, atm, sundir, &
          albedo_thermal=zero, albedo_solar=Ag, &
          lthermal=.False., lsolar=.True., &
          edir=edir, edn=edn, eup=eup, abso=abso)
      endif

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
    real(ireals) :: Ag
    logical :: lthermal, lsolar
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

    Ag = .1
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

    lsolar = .True.
    lthermal = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)

    default_options=''
    default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
    !default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

    !default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_abso_direct hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_abso_diffuse hdf5:'//trim(outfile)//'::append'

    default_options=trim(default_options)//' -plexrt_dump_thermal_Edn_2_ke1 hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_thermal_Eup_2_ke1 hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_thermal_abso hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_Edir_2_ke1 hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_Edn_2_ke1 hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_Eup_2_ke1 hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_abso hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_lwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_iwc hdf5:'//trim(outfile)//'::append'
    default_options=trim(default_options)//' -plexrt_dump_temp hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex4(PETSC_COMM_WORLD, gridfile, icondatafile, Ag, lthermal, lsolar)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
