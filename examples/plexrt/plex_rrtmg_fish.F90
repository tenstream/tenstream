module m_plex_rrtmg_fish

#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options

use m_helper_functions, only: CHKERR, CHKWARN, imp_bcast, determine_normal_direction, &
  spherical_2_cartesian, angle_between_two_vec, rad2deg, deg2rad, reverse, itoa, cstr, &
  meanval, meanvec

use m_helper_functions, only: cross_3d, rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world, &
  rotate_angle_x, rotation_matrix_around_axis_vec

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

use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

use m_netcdfio, only : ncload, ncwrite

use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, print_tenstr_atm, &
  reff_from_lwc_and_N, hydrostat_plev

use m_icon_plex_utils, only: create_2d_fish_plex, create_2d_regular_plex, &
  dmplex_2D_to_3D, dump_ownership, Nz_Ncol_vec_to_celldm1

implicit none

logical, parameter :: ldebug=.True.

  contains

    subroutine ex_plex_rrtmg_fish(comm, Nx, Ny, Nz, dx, dz, Ag, lthermal, lsolar)
      MPI_Comm, intent(in) :: comm
      integer(iintegers), intent(in) :: Nx, Ny, Nz
      real(ireals), intent(in) :: dx, dz, Ag
      logical, intent(in) :: lthermal, lsolar

      type(tDM) :: dm2d, dm2d_dist, dm3d
      type(tPetscSF) :: migration_sf
      real(ireals) :: hhl(Nz)

      integer(mpiint) :: myid, numnodes, ierr

      real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

      type(t_plexgrid), allocatable :: plex
      integer(iintegers), allocatable :: zindex(:)
    class(t_plex_solver), allocatable :: solver

      type(t_tenstr_atm) :: atm
      character(len=default_str_len) :: atm_filename
      integer(iintegers) :: k, fStart, fEnd, Ncol, dNlev, dNlay, Nlev
      real(ireals), allocatable :: col_plev(:,:), col_tlev(:,:), dp(:)
      real(ireals), allocatable :: col_lwc(:,:), col_reff(:,:)
      real(ireals), parameter :: lapse_rate=9.5e-3, Tsrfc=288.3, minTemp=Tsrfc-60._ireals

      integer(iintegers) :: solve_iterations, iter
      real(ireals) :: solve_iterations_scale, lwcval
      logical :: lregular_mesh, lflg

      real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      lregular_mesh = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        '-use_regular_mesh', lregular_mesh, lflg, ierr) ; call CHKERR(ierr)
      if(lregular_mesh) then
        call create_2d_regular_plex(comm, Nx, Ny, dm2d, dm2d_dist, migration_sf, opt_dx=dx)
      else
        call create_2d_fish_plex(comm, Nx, Ny, dm2d, dm2d_dist, migration_sf, opt_dx=dx)
      endif

      hhl(1) = zero
      do k=2,Nz
        hhl(k) = hhl(k-1) + dz
      enddo

      call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
      Ncol = fEnd - fStart
      dNlev = Nz; dNlay = dNlev-1

      if(myid.eq.0) print *,'Dynamics Grid has Size Nlev, Ncol:', dNlev, Ncol
      if(Ncol.eq.0) call CHKERR(1_mpiint, 'We have a process that has nothing to do? '// &
        'Maybe decrease the number of processes or increase the problem size')

      ! prepare atmosphere
      allocate(col_tlev(dNlev, Ncol))
      allocate(col_plev(dNlev, Ncol))

      do k = 1, dNlev
        col_tlev(k, i1) = max(minTemp, Tsrfc - (lapse_rate * hhl(k)))
      enddo

      allocate(dp(dNlay))
      call hydrostat_plev(1013._ireals, meanvec(col_tlev(:,i1)), hhl, &
        col_plev(:,i1), dp)
      deallocate(dp)

      if(myid.eq.0) then
        print *,'Dynamics Grid Pressure and Temperature'
        do k = 1, dNlev
          print *,k, col_plev(k, i1), col_tlev(k, i1)
        enddo
      endif

      do k = 2, Ncol
        col_plev(:, k) =  col_plev(:, i1)
        col_tlev(:, k) =  col_tlev(:, i1)
      enddo

      allocate(col_lwc(dNlay, Ncol), col_reff(dNlay, Ncol))
      col_lwc = 0
      col_reff = 10

      lwcval = .1
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-lwc", lwcval, lflg,ierr) ; call CHKERR(ierr)

      col_lwc(1+dNlay/2:dNlay, 1:Ncol:10) = lwcval
      col_lwc(1+dNlay/2:dNlay, 2:Ncol:10) = lwcval
      col_lwc(1+dNlay/2:dNlay, 3:Ncol:10) = lwcval
      col_lwc(1+dNlay/2:dNlay, 4:Ncol:10) = lwcval
      col_lwc(1+dNlay/2:dNlay, 5:Ncol:10) = lwcval

      call init_data_strings()

      call setup_tenstr_atm(comm, .False., atm_filename, &
        col_plev, col_tlev, atm, &
        d_lwc=col_lwc, d_reliq=col_reff)

      Nlev = size(atm%plev,1,kind=iintegers)
      call dmplex_2D_to_3D(dm2d_dist, Nlev, reverse(atm%zt(:, i1)), dm3d, zindex)

      call dump_ownership(dm3d, '-dump_ownership', '-show_plex')

      call setup_plexgrid(dm2d_dist, dm3d, Nlev-1, zindex, plex, hhl=reverse(atm%zt(:, i1)))

      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)
      deallocate(zindex)

      if(lregular_mesh) then
        call allocate_plexrt_solver_from_commandline(solver, 'rectilinear_5_8')
      else
        call allocate_plexrt_solver_from_commandline(solver, '5_8')
      endif
      call init_plex_rt_solver(plex, solver)

      call init_sundir()

      solve_iterations = 1
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solve_iterations", &
        solve_iterations, lflg,ierr) ; call CHKERR(ierr)
      solve_iterations_scale = 0
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solve_iterations_scale", &
        solve_iterations_scale, lflg,ierr) ; call CHKERR(ierr)

      do iter = 1, solve_iterations
        if(lthermal) then

          ! Perform a circular shift of cloud field by one column
          col_lwc = cshift(col_lwc, shift=1, dim=2)
          call setup_tenstr_atm(comm, .False., atm_filename, &
            col_plev, col_tlev, atm, &
            d_lwc=col_lwc, d_reliq=col_reff)

          call plexrt_rrtmg(solver, atm, sundir, &
            albedo_thermal=zero, albedo_solar=Ag, &
            lthermal=.True., lsolar=.False., &
            edir=edir, edn=edn, eup=eup, abso=abso)
          if(myid.eq.0) then
            print *, ''
            print *, cstr('Thermal Radiation', 'green')
            print *,cstr('Avg.horiz','blue')//&
              cstr(' k  Edn             Eup           abso', 'blue')
            do k = 1, ubound(abso,1)
              print *,k, meanval(edn(k,:)), meanval(eup(k,:)), meanval(abso(k,:))
            enddo
            print *,k, meanval(edn(k,:)), meanval(eup(k,:))
          endif
          if(meanval(edn).lt.epsilon(edn)) &
            call CHKERR(1_mpiint, 'Mean Edn is really small,'// &
            'the solver probably had a problem but did not fail.'// &
            'Caution! Your results may be garbage!')
        endif

        if(lsolar) then
          call init_sundir()
          call sundir_rot_theta(sundir, real(iter-1, ireals) * solve_iterations_scale)
          call plexrt_rrtmg(solver, atm, sundir, &
            albedo_thermal=zero, albedo_solar=Ag, &
            lthermal=.False., lsolar=.True., &
            edir=edir, edn=edn, eup=eup, abso=abso)

          if(myid.eq.0) then
            print *, ''
            print *, cstr('Solar Radiation', 'green')
            print *, cstr('Avg.horiz','blue')//&
              cstr(' k   Edir          Edn            Eup            abso', 'blue')
            do k = 1, ubound(abso,1)
              print *,k, meanval(edir(k,:)), meanval(edn(k,:)), meanval(eup(k,:)), meanval(abso(k,:))
            enddo
            print *,k, meanval(edir(k,:)), meanval(edn(k,:)), meanval(eup(k,:))
          endif
          if(meanval(edn).lt.epsilon(edn)) &
            call CHKERR(1_mpiint, 'Mean Edn is really small,'// &
            'the solver probably had a problem but did not fail.'// &
            'Caution! Your results may be garbage!')
        endif
      enddo ! solve_iterations

      call dump_result()

      call destroy_plexrt_rrtmg(solver, lfinalizepetsc=.False.)

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
          logical :: lflg
          real(ireals) :: first_normal(3)
          integer(mpiint) :: myid, ierr
          real(ireals) :: rot_angle, Mrot(3,3), rot_sundir(3)
          integer(iintegers) :: argcnt

          call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
          first_normal = get_normal_of_first_TOA_face(solver%plex)

          sundir = zero
          if(ldebug.and.myid.eq.0) print *,myid, 'determine initial sundirection ...'
          argcnt = 3
          call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir", &
            sundir, argcnt, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            call CHKERR(int(argcnt-i3, mpiint), "must provide 3 values for sundir, comma separated, no spaces")
          else
            sundir = first_normal + [zero, +.5_ireals, zero]
          endif
          sundir = sundir/norm2(sundir)

          if(ldebug.and.myid.eq.0) print *,'Initial sundirection = ', sundir, &
            ': sza', angle_between_two_vec(sundir, first_normal), 'rad'
          if(myid.eq.0) print *,'Initial sundirection = ', sundir, &
            ': sza', rad2deg(angle_between_two_vec(sundir, first_normal)),'deg'

          rot_angle = zero
          call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_phi", &
            rot_angle, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), first_normal)
            rot_sundir = matmul(Mrot, sundir)
            if(myid.eq.0) print *,'rotated sundirection = ', rot_sundir, &
              ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)),'deg'
            sundir = rot_sundir
          endif

          rot_angle = zero
          call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-sundir_rot_theta", &
            rot_angle, lflg, ierr) ; call CHKERR(ierr)
          if(lflg) then
            call sundir_rot_theta(sundir, rot_angle)
          endif

          if(ldebug.and.myid.eq.0) print *,'determine initial sundirection ... done'
        end subroutine
        subroutine sundir_rot_theta(sundir, rot_angle)
          real(ireals) :: sundir(3), rot_angle, Mrot(3,3), U(3), rot_sundir(3)
          integer(mpiint) :: myid, ierr
          real(ireals) :: first_normal(3)

          call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
          first_normal = get_normal_of_first_TOA_face(solver%plex)
          U = cross_3d(first_normal, sundir)
          Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), U)
          rot_sundir = matmul(Mrot, sundir)
          if(myid.eq.0) print *,'rotated sundirection = ', rot_sundir, ': sza', &
            rad2deg(angle_between_two_vec(rot_sundir, first_normal)),'deg'
          sundir = rot_sundir
        end subroutine
        subroutine init_data_strings()
          logical :: lflg
          atm_filename='afglus_100m.dat'
          call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atm_filename', &
            atm_filename, lflg, ierr); call CHKERR(ierr)
        end subroutine
    end subroutine

  end module

  program main
    use m_plex_rrtmg_fish
    implicit none

    character(len=default_str_len) :: outfile
    logical :: lflg, lthermal, lsolar
    integer(mpiint) :: ierr
    character(len=10*default_str_len) :: default_options
    integer(iintegers) :: Nx, Ny, Nz
    real(ireals) :: dx, dz, Ag


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
    dx = 1e3_ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg,ierr) ; call CHKERR(ierr)
    dz = 1e2_ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg,ierr) ; call CHKERR(ierr)
    Ag = .1_ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)
    lsolar = .True.
    lthermal = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)


    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
    if(.not.lflg) stop 'need to supply a output filename... please call with -out <fname_of_output_file.h5>'

    default_options='-polar_coords no'
    default_options=trim(default_options)//' -show_plex hdf5:'//trim(outfile)
    !default_options=trim(default_options)//' -show_ownership hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_abso hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_iconindex hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_zindex hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_domainboundary hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_lwc hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_iwc hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_edir hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_srcVec hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_DiffSrcVec hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_fV2cV_ediff hdf5:'//trim(outfile)//'::append'
    !default_options=trim(default_options)//' -show_WedgeOrient hdf5:'//trim(outfile)//'::append'

    print *,'Adding default Petsc Options:', trim(default_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr)

    call ex_plex_rrtmg_fish(PETSC_COMM_WORLD, Nx, Ny, Nz, dx, dz, Ag, lthermal, lsolar)

    call mpi_barrier(PETSC_COMM_WORLD, ierr)
    call PetscFinalize(ierr)
  end program
