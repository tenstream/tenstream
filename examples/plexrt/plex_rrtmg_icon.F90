module m_examples_plex_ex_rrtmg_icon

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only: &
    & angle_between_two_vec, &
    & CHKERR, &
    & deg2rad, &
    & get_petsc_opt, &
    & meanvec, &
    & rad2deg, &
    & reverse

  use m_data_parameters, only: &
    & i0, i1, i2, i3, i4, i5,  &
    & default_str_len, &
    & iintegers, &
    & imp_ireals, &
    & init_mpi_data_parameters, &
    & ireal_dp, &
    & ireals, &
    & mpiint, &
    & one, &
    & zero

  use m_icon_plex_utils, only: gen_2d_plex_from_icongridfile, icon_hdcp2_default_hhl, &
                               dump_ownership, dmplex_2D_to_3D, icon_ncvec_to_plex, celldm_veccopy, celldm1_vec_to_nz_ncol, &
                               dm2d_vec_to_Nz_Ncol

  use m_plex_grid, only: t_plexgrid, setup_plexgrid, get_normal_of_first_toa_face, atm_dz_to_vertex_heights

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: compute_face_geometry, init_plex_rt_solver

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, print_tenstr_atm, &
                                reff_from_lwc_and_N, hydrostat_plev

  use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

  use m_netcdfio, only: ncload

  implicit none

  logical, parameter :: ldebug = .true.

contains

  subroutine plex_ex_rrtmg_icon(comm, gridfile, icondatafile, Ag, lthermal, lsolar)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: gridfile, icondatafile
    real(ireals), intent(in) :: Ag
    logical, intent(in) :: lthermal, lsolar

    integer(mpiint) :: myid, numnodes, ierr
    type(tDM) :: dm2d, dm2d_dist, dm3d
    type(tPetscSF) :: migration_sf
    type(tPetscSection) :: par_cell_Section
    AO, allocatable :: cell_ao_2d
    type(t_plexgrid), allocatable :: plex
    !type(tVec), allocatable :: lwcvec2d, iwcvec2d
    real(ireals), allocatable :: col_plev(:, :), col_tlev(:, :), col_qv(:, :), dp(:)
    real(ireals), allocatable :: col_lwc(:, :), col_reliq(:, :), col_qnc(:, :)
    real(ireals), allocatable :: col_iwc(:, :), col_reice(:, :), col_qni(:, :)
    real(ireals), parameter :: lapse_rate = 6.5e-3, Tsrfc = 288.3, minTemp = Tsrfc - 50._ireals
    character(len=default_str_len) :: lwc_data_string, qnc_data_string, iwc_data_string, qni_data_string
    character(len=default_str_len) :: qv_data_string, atm_filename

    type(t_tenstr_atm) :: atm
    integer(iintegers), allocatable :: zindex(:)
    integer(iintegers) :: k, fStart, fEnd, Ncol, dNlev, dNlay, Nlev
    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso

    real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system
    real(ireals) :: coord_displacement(3) ! cartesian displacement of origin

    class(t_plex_solver), allocatable :: solver

    type(tVec), allocatable :: lwcvec, iwcvec, qncvec, qnivec, qvvec

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    call gen_2d_plex_from_icongridfile(comm, gridfile, dm2d, dm2d_dist, &
                                       migration_sf, cell_ao_2d, coord_displacement)

    call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart
    dNlev = size(icon_hdcp2_default_hhl); dNlay = dNlev - 1

    if (myid .eq. 0) print *, 'Dynamics Grid has Size Nlev, Ncol:', dNlev, Ncol

    ! prepare atmosphere
    allocate (col_tlev(dNlev, Ncol))
    allocate (col_plev(dNlev, Ncol))

    do k = 1, dNlev
      col_tlev(k, i1) = max(minTemp, Tsrfc - (lapse_rate * icon_hdcp2_default_hhl(dNlev - (k - 1))))
    end do

    allocate (dp(dNlay))
    call hydrostat_plev(1013._ireals, meanvec(col_tlev(:, i1)), reverse(icon_hdcp2_default_hhl), &
                        col_plev(:, i1), dp)
    deallocate (dp)

    if (myid .eq. 0) then
      print *, 'Dynamics Grid Pressure and Temperature'
      do k = 1, dNlev
        print *, k, col_plev(k, i1), col_tlev(k, i1)
      end do
    end if

    do k = 2, Ncol
      col_plev(:, k) = col_plev(:, i1)
      col_tlev(:, k) = col_tlev(:, i1)
    end do

    call init_data_strings()

    call setup_tenstr_atm(comm, .false., atm_filename, col_plev, col_tlev, atm)
    Nlev = size(atm%plev, 1, kind=iintegers)

    !call print_tenstr_atm(atm)

    ! Setup 3D DMPLEX grid
    call dmplex_2D_to_3D(dm2d_dist, Nlev, reverse(atm%zt(:, i1)), -coord_displacement, dm3d, zindex)
    ! call atm_dz_to_vertex_heights(atm%dz, dm3d)

    call dump_ownership(dm3d, '-dump_ownership', '-show_plex')
    call setup_plexgrid(dm2d_dist, dm3d, Nlev - 1, zindex, plex, reverse(atm%zt(:, i1)))

    !Load Data from iconfile and distribute it
    if (myid .eq. 0) print *, 'Read data from icondatafile ', trim(icondatafile)

    call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, qnc_data_string, par_cell_Section, qncvec)
    call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, lwc_data_string, par_cell_Section, lwcvec)
    call PetscObjectViewFromOptions(PetscObjectCast(lwcvec), PETSC_NULL_OBJECT, '-show_lwc', ierr); call CHKERR(ierr)

    call dm2d_vec_to_Nz_Ncol(par_cell_Section, lwcvec, col_lwc); col_lwc = col_lwc * 1e3
    call VecDestroy(lwcvec, ierr); call CHKERR(ierr)
    call dm2d_vec_to_Nz_Ncol(par_cell_Section, qncvec, col_qnc); col_qnc = col_qnc * 1e-3
    call VecDestroy(qncvec, ierr); call CHKERR(ierr)
    allocate (col_reliq(dNlay, Ncol))
    col_reliq = min(40._ireals, max(2.5_ireals, reff_from_lwc_and_N(col_lwc, col_qnc)))
    deallocate (col_qnc)
    if (myid .eq. 0) print *, 'Min/Max Liquid effective Radius', minval(col_reliq), maxval(col_reliq)

    call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, qni_data_string, par_cell_Section, qnivec)
    call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, iwc_data_string, par_cell_Section, iwcvec)
    call PetscObjectViewFromOptions(PetscObjectCast(iwcvec), PETSC_NULL_OBJECT, '-show_iwc', ierr); call CHKERR(ierr)

    call dm2d_vec_to_Nz_Ncol(par_cell_Section, iwcvec, col_iwc); col_iwc = col_iwc * 1e3
    call VecDestroy(iwcvec, ierr); call CHKERR(ierr)
    call dm2d_vec_to_Nz_Ncol(par_cell_Section, qnivec, col_qni); col_qni = col_qni * 1e-3
    call VecDestroy(qnivec, ierr); call CHKERR(ierr)
    allocate (col_reice(dNlay, Ncol))
    ! k == .8 is for clean air
    col_reice = min(120._ireals, max(5._ireals, reff_from_lwc_and_N(col_iwc, col_qni, k=.8_ireals)))
    deallocate (col_qni)
    if (myid .eq. 0) print *, 'Min/Max Ice effective Radius', minval(col_reice), maxval(col_reice)

    call icon_ncvec_to_plex(dm2d, dm2d_dist, migration_sf, icondatafile, qv_data_string, par_cell_Section, qvvec)
    call dm2d_vec_to_Nz_Ncol(par_cell_Section, qvvec, col_qv)
    call VecDestroy(qvvec, ierr); call CHKERR(ierr)

    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          col_plev, col_tlev, atm, &
                          d_h2ovmr=col_qv, &
                          d_lwc=col_lwc, d_reliq=col_reliq, &
                          d_iwc=col_iwc, d_reice=col_reice)

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)

    call init_sundir()

    if (lthermal) then
      call plexrt_rrtmg(solver, atm, sundir, &
                        albedo_thermal=zero, albedo_solar=Ag, &
                        lthermal=.true., lsolar=.false., &
                        edir=edir, edn=edn, eup=eup, abso=abso)
    end if

    if (lsolar) then
      call plexrt_rrtmg(solver, atm, sundir, &
                        albedo_thermal=zero, albedo_solar=Ag, &
                        lthermal=.false., lsolar=.true., &
                        edir=edir, edn=edn, eup=eup, abso=abso)
    end if

    call destroy_plexrt_rrtmg(solver, lfinalizepetsc=.false.)

  contains
    subroutine init_data_strings()
      logical :: lflg
      lwc_data_string = 'clw'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-lwc_data_string', lwc_data_string, lflg, ierr); call CHKERR(ierr)
      iwc_data_string = 'cli'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-iwc_data_string', iwc_data_string, lflg, ierr); call CHKERR(ierr)
      qnc_data_string = 'qnc'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-qnc_data_string', qnc_data_string, lflg, ierr); call CHKERR(ierr)
      qni_data_string = 'qni'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-qni_data_string', qni_data_string, lflg, ierr); call CHKERR(ierr)
      qv_data_string = 'qv'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-qv_data_string', qv_data_string, lflg, ierr); call CHKERR(ierr)
      atm_filename = 'atm.dat'
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)
    end subroutine
    subroutine init_sundir()
      use m_helper_functions, only: cross_3d, rotation_matrix_world_to_local_basis, rotation_matrix_local_basis_to_world, &
                                    rotate_angle_x, rotation_matrix_around_axis_vec
      logical :: lflg
      PetscBool :: lflg_p
      real(ireals) :: first_normal(3)
      integer(mpiint) :: myid, ierr
      integer(iintegers) :: nargs
      real(ireals) :: rot_angle, Mrot(3, 3), U(3), V(3), rot_sundir(3)

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      first_normal = get_normal_of_first_TOA_face(solver%plex)

      if (ldebug .and. myid .eq. 0) print *, myid, 'determine initial sundirection ...'
      nargs = i3
      call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                                    "-sundir", sundir, nargs, lflg_p, ierr); call CHKERR(ierr)
      if (lflg_p) then
        call CHKERR(int(nargs - i3, mpiint), 'must provide exactly 3 values for -sundir. '// &
                    'Need to be given comma separated without spaces')
      else
        sundir = first_normal + [zero, zero, 1e-2_ireals]
        sundir = -[0.5368672026070134, 0.7320371297951827, 0.419398673538855] ! example on cut50
        sundir = [-0.67442513891627376, 0.68211997954408909, -0.28260054052414002] ! example from ifc_icon run at 10863 sec
      end if
      sundir = sundir / norm2(sundir)

      if (ldebug .and. myid .eq. 0) &
        & print *, 'Initial sundirection =', sundir, &
                 & ': sza', angle_between_two_vec(sundir, first_normal), 'rad'
      if (myid .eq. 0) &
        & print *, 'Initial sundirection = ', sundir, &
                 & ': sza', rad2deg(angle_between_two_vec(sundir, first_normal)), 'deg'

      rot_angle = zero
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-sundir_rot_phi", rot_angle, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        U = [first_normal(2), -first_normal(1), zero]
        U = U / norm2(U)
        V = cross_3d(first_normal, U)
        Mrot = rotation_matrix_world_to_local_basis(first_normal, U, V)
        rot_sundir = matmul(Mrot, sundir)
        rot_sundir = rotate_angle_x(rot_sundir, rot_angle)
        if (ldebug .and. myid .eq. 0) print *, 'U', U
        if (ldebug .and. myid .eq. 0) print *, 'V', V
        if (ldebug .and. myid .eq. 0) print *, 'rot_sundir', rot_sundir
        Mrot = rotation_matrix_local_basis_to_world(first_normal, U, V)
        rot_sundir = matmul(Mrot, rot_sundir)
        if (myid .eq. 0) print *, 'rotated sundirection =', rot_sundir, &
          ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)), 'deg'
        sundir = rot_sundir
      end if

      rot_angle = zero
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-sundir_rot_theta", rot_angle, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        if (dot_product(first_normal, sundir) .gt. one - epsilon(U)) then
          first_normal(1) = first_normal(1) + epsilon(first_normal)
          first_normal(2) = first_normal(2) - epsilon(first_normal)
          first_normal = first_normal / norm2(first_normal)
        end if
        U = cross_3d(first_normal, sundir)
        Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), U)
        rot_sundir = matmul(Mrot, sundir)
        if (ldebug .and. myid .eq. 0) print *, 'S', sundir, norm2(sundir)
        if (ldebug .and. myid .eq. 0) print *, 'U', U, norm2(U)
        if (ldebug .and. myid .eq. 0) print *, 'rot_sundir', rot_sundir
        if (myid .eq. 0) print *, 'rotated sundirection = ', rot_sundir, &
          ': sza', rad2deg(angle_between_two_vec(rot_sundir, first_normal)), 'deg'
        sundir = rot_sundir
      end if

      call mpi_bcast(sundir, 3_mpiint, imp_ireals, 0_mpiint, comm, ierr); call CHKERR(ierr)
      if (ldebug .and. myid .eq. 0) print *, 'determine initial sundirection ... done'
    end subroutine
  end subroutine
end module
