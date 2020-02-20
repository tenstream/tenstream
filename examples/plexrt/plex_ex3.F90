module m_mpi_plex_ex3

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg, deg2rad, reverse, itoa, cstr, &
  meanval

use m_data_parameters, only : ireals, iintegers, mpiint, &
  default_str_len, &
  i0, i1, i2, i3, i4, i5,  &
  zero, one,       &
  init_mpi_data_parameters

use m_icon_grid, only: t_icongrid, read_icon_grid_file, &
  bcast_icongrid, distribute_icon_grid

use m_plex_grid, only: t_plexgrid, &
  setup_edir_dmplex, setup_abso_dmplex, compute_face_geometry, &
  ncvar2d_to_globalvec, setup_plexgrid, setup_cell1_dmplex, &
  gen_test_mat, get_normal_of_first_toa_face, get_horizontal_faces_around_vertex, &
  atm_dz_to_vertex_heights

use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

use m_plex_rt, only: compute_face_geometry, &
  init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
  plexrt_get_result, destroy_plexrt_solver

use m_netcdfio, only : ncload, ncwrite

use m_icon_plex_utils, only: create_2d_fish_plex, create_2d_regular_plex, &
  dmplex_2D_to_3D, dump_ownership, Nz_Ncol_vec_to_celldm1

implicit none

logical, parameter :: ldebug=.False.

  contains

    subroutine plex_ex3(comm, Nx, Ny, Nz, dz, Ag, lverbose)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nx, Ny, Nz
      real(ireals), intent(in) :: dz, Ag
      logical, intent(in) :: lverbose

      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: migration_sf
      real(ireals) :: hhl(Nz)

      integer(mpiint) :: myid, numnodes, ierr
      integer(iintegers) :: icol, k
      !type(tVec) :: lwcvec, iwcvec

      real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

      type(t_plexgrid), allocatable :: plex
      integer(iintegers), allocatable :: zindex(:)
      class(t_plex_solver), allocatable :: solver

      real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

      logical, parameter :: lthermal = .False.
      logical :: lregular_mesh, lflg

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      lregular_mesh = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        '-use_regular_mesh', lregular_mesh, lflg, ierr) ; call CHKERR(ierr)
      if(lregular_mesh) then
        call create_2d_regular_plex(comm, Nx, Ny, dm2d, dm2d_dist, migration_sf, lverbose=lverbose)
      else
        call create_2d_fish_plex(comm, Nx, Ny, dm2d, dm2d_dist, migration_sf, lverbose=lverbose)
      endif

      hhl(1) = zero
      do k=2,Nz
        hhl(k) = hhl(k-1) + dz
      enddo
      hhl = reverse(hhl)

      call dmplex_2D_to_3D(dm2d_dist, Nz, hhl, dm3d, zindex, lverbose=lverbose)

      call setup_plexgrid(dm2d_dist, dm3d, Nz-1, zindex, plex, hhl)
      deallocate(zindex)

      call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

      if(lregular_mesh) then
        call allocate_plexrt_solver_from_commandline(solver, 'rectilinear_5_8')
      else
        call allocate_plexrt_solver_from_commandline(solver, '5_8')
      endif
      call init_plex_rt_solver(plex, solver)

      call init_sundir()

      call set_plex_rt_optprop(solver, vert_integrated_kabs=1e-0_ireals, vert_integrated_ksca=one, lverbose=lverbose)

      if(.not.allocated(solver%albedo)) then
        allocate(solver%albedo)
        call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
      endif
      call VecSet(solver%albedo, Ag, ierr); call CHKERR(ierr)

      if(lthermal) then
        if(.not.allocated(solver%plck)) then
          allocate(solver%plck)
          call DMCreateGlobalVector(solver%plex%horizface1_dm, solver%plck, ierr); call CHKERR(ierr)
        endif
        call VecSet(solver%plck, 100._ireals, ierr); call CHKERR(ierr)

        call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=sundir)
      else
        call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., sundir=sundir)
      endif

      call plexrt_get_result(solver, edn, eup, abso, redir=edir)

      if(lverbose) then
        do icol = 1, Nx
        print *,''
        print *,cstr('Column ','blue')//cstr(itoa(icol),'red')//&
          cstr('  k  Edir              Edn              Eup              abso', 'blue')
        do k = 1, ubound(abso,1)
        print *,k, edir(k,icol), edn(k,icol), eup(k,icol), abso(k,icol)
        enddo
        print *,k, edir(k,icol), edn(k,icol), eup(k,icol)
        enddo

        print *, ''
        print *, 'Averages'
        do k = 1, ubound(abso,1)
        print *,k, meanval(edir(k,:)), meanval(edn(k,:)), meanval(eup(k,:)), meanval(abso(k,:))
        enddo
        print *,k, meanval(edir(k,:)), meanval(edn(k,:)), meanval(eup(k,:))
      endif

      call dump_result()

      call destroy_plexrt_solver(solver, lfinalizepetsc=.False.)

      contains
        subroutine dump_result()
          logical :: lflg, ldump
          type(tVec) :: vec=PETSC_NULL_VEC
          call DMGetGlobalVector(solver%plex%cell1_dm, vec, ierr); call CHKERR(ierr)

          ldump = .False.
          call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-ex_dump_result', &
            ldump, lflg, ierr) ; call CHKERR(ierr)
          if(ldump) then
            call Nz_Ncol_vec_to_celldm1(solver%plex, edir(2:,:), vec)
            call PetscObjectSetName(vec, 'ex_dump_result_edir', ierr);call CHKERR(ierr)
            call PetscObjectViewFromOptions(solver%plex%cell1_dm, PETSC_NULL_DM, "-ex_dump_result_edir_dm", ierr); call CHKERR(ierr)
            call PetscObjectViewFromOptions(vec, PETSC_NULL_DM, "-ex_dump_result_edir", ierr); call CHKERR(ierr)
          endif
          call DMRestoreGlobalVector(solver%plex%cell1_dm, vec, ierr); call CHKERR(ierr)
        end subroutine
        subroutine init_sundir()
          use m_helper_functions, only: cross_3d, rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world, &
            rotate_angle_x, rotation_matrix_around_axis_vec
          logical :: lflg
          real(ireals) :: first_normal(3)
          integer(mpiint) :: myid, ierr
          integer(iintegers) :: nargs
          real(ireals) :: rot_angle, Mrot(3,3), U(3), rot_sundir(3)

          call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
          first_normal = get_normal_of_first_TOA_face(solver%plex)

          if(ldebug.and.myid.eq.0) print *,myid, 'determine initial sundirection ...'
          nargs = i3
          call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
            "-sundir", sundir, nargs, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            call CHKERR(int(nargs-i3, mpiint), 'must provide exactly 3 values for -sundir. '// &
              'Need to be given comma separated without spaces')
          else
            sundir = first_normal + [zero, -.5_ireals, zero]
            sundir = sundir/norm2(sundir)
          endif
          sundir = sundir/norm2(sundir)

          if(lverbose.and.myid.eq.0) then
            print *,'Initial sundirection = ', sundir, ': sza', angle_between_two_vec(sundir, first_normal), 'rad'
            print *,'Initial sundirection = ', sundir, ': sza', rad2deg(angle_between_two_vec(sundir, first_normal)),'deg'
          endif


          rot_angle = zero
          call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_phi", &
            rot_angle, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), first_normal)
            rot_sundir = matmul(Mrot, sundir)
            if(ldebug.and.myid.eq.0) print *,'rot_sundir', rot_sundir
            if(myid.eq.0) &
              print *,'rotated sundirection = ', rot_sundir, ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)),'deg'
            sundir = rot_sundir
          endif

          rot_angle = zero
          call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_theta", &
            rot_angle, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            U = cross_3d(first_normal, sundir)
            Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), U)
            rot_sundir = matmul(Mrot, sundir)
            if(ldebug.and.myid.eq.0) print *,'S', sundir, norm2(sundir)
            if(ldebug.and.myid.eq.0) print *,'U', U, norm2(U)
            if(ldebug.and.myid.eq.0) print *,'rot_sundir', rot_sundir
            if(lverbose.and.myid.eq.0) &
              print *,'rotated sundirection = ', rot_sundir, ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)),'deg'
            sundir = rot_sundir
          endif

          if(ldebug.and.myid.eq.0) print *,'determine initial sundirection ... done'
        end subroutine
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
    logical :: lverbose

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
    call init_mpi_data_parameters(PETSC_COMM_WORLD)
    call read_commandline_options(PETSC_COMM_WORLD)

    Nx = 2
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg,ierr) ; call CHKERR(ierr)
    Ny = 3
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg,ierr) ; call CHKERR(ierr)
    Nz = 2
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg,ierr) ; call CHKERR(ierr)
    dz = one/real(Nz, ireals)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)
    Ag = .1_ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

    lverbose = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-verbose', &
      lverbose, lflg, ierr) ; call CHKERR(ierr)

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

    if(lverbose) print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call plex_ex3(PETSC_COMM_WORLD, Nx, Ny, Nz, dz, Ag, lverbose)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
