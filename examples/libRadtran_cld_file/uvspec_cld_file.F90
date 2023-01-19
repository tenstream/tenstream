module m_example_uvspec_cld_file
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero, pprts_get_result_toZero

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, &
                               i0, i1, i2, zero, one, default_str_len

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, &
                                hydrostat_dp, load_atmfile, t_bg_atm, print_tenstr_atm

  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: &
    & CHKERR, &
    & domain_decompose_2d_petsc, &
    & get_petsc_opt, &
    & imp_bcast, &
    & resize_arr, &
    & reverse, &
    & spherical_2_cartesian
  use m_netcdfio, only: ncload, ncwrite, get_global_attribute, set_global_attribute, set_attribute

  use m_petsc_helpers, only: getvecpointer, restorevecpointer

  use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D, &
                               rank0_f90vec_to_plex, dmplex_gvec_from_f90_array, plex_gvec_tozero

  use m_plex_grid, only: t_plexgrid, &
                         setup_plexgrid, create_plex_section

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: init_plex_rt_solver

  use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

  use m_tenstream_interpolation, only: interp_1d
  use m_search, only: find_real_location
  implicit none

contains
  subroutine example_uvspec_cld_file_pprts(specint, comm, &
                                           cldfile, atm_filename, outfile, &
                                           albedo_th, albedo_sol, &
                                           lsolar, lthermal, &
                                           phi0, theta0, &
                                           scene_shift_x, scene_shift_y, scene_shift_it)

    character(len=*), intent(in) :: specint                  ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    integer(iintegers), intent(in) :: scene_shift_x, scene_shift_y, scene_shift_it

    real(ireals), dimension(:, :, :), allocatable, target :: lwc, reliq ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:, :, :), allocatable, target :: plev, tlev ! will have local shape nzp+1, nxp, nyp
    real(ireals), dimension(:), allocatable :: hhl ! dim Nz+1
    real(ireals), pointer :: z(:, :, :, :) => null(), z1d(:) => null() ! dim Nz+1

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    type(t_bg_atm), allocatable :: bg_atm

    class(t_solver), allocatable :: pprts_solver

    real(ireals) :: dx, dy, z_location
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: is, ie, js, je
    integer(iintegers) :: nxp, nyp, nzp ! local sizes of domain, nzp being number of layers
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    integer(iintegers) :: k

    call mpi_comm_rank(comm, myid, ierr)
    call load_input(comm, cldfile, dx, dy, is, ie, js, je, hhl, lwc, reliq, ierr); call CHKERR(ierr)

    ! Determine Domain Decomposition
    call domain_decompose_2d_petsc(comm, &
      & Nx_global=size(lwc, dim=2, kind=iintegers), &
      & Ny_global=size(lwc, dim=3, kind=iintegers), &
      & Nx_local=nxp, &
      & Ny_local=nyp, &
      & xStart=is, &
      & yStart=js, &
      & nxproc=nxproc, &
      & nyproc=nyproc, &
      & ierr=ierr); call CHKERR(ierr)
    is = is + 1 ! fortran based indices
    js = js + 1 ! fortran based indices
    ie = is + nxp - 1
    je = js + nyp - 1

    print *, 'myid', myid, 'is,ie', is, ie, 'js,je', js, je

    ! Load the background atmosphere file and interpolate pressure and temperature from that
    nzp = size(hhl) - 1
    allocate (plev(nzp + 1, nxp, nyp), tlev(nzp + 1, nxp, nyp))
    call load_atmfile(comm, atm_filename, bg_atm)
    do k = 1, nzp + 1
      z_location = find_real_location(bg_atm%zt, hhl(k))
      plev(k, :, :) = interp_1d(z_location, bg_atm%plev)
      tlev(k, :, :) = interp_1d(z_location, bg_atm%tlev)
    end do
    deallocate (bg_atm)

    if (myid .eq. 0) then
      do k = 1, nzp + 1
        print *, k, 'plev', plev(k, is, js), 'Tlev', tlev(k, is, js)
      end do
    end if

    call allocate_pprts_solver_from_commandline(pprts_solver, default_solver='3_10', ierr=ierr); call CHKERR(ierr)

    do k = 1, scene_shift_it
      lwc = cshift(lwc, scene_shift_x, dim=2)
      reliq = cshift(reliq, scene_shift_x, dim=2)
      lwc = cshift(lwc, scene_shift_y, dim=3)
      reliq = cshift(reliq, scene_shift_y, dim=3)
      call run_rrtmg_lw_sw(&
        & specint, &
        & pprts_solver, &
        & nxproc, nyproc, &
        & atm_filename, &
        & dx, dy, phi0, theta0, &
        & plev, tlev, &
        & lwc(:, is:ie, js:je), &
        & reliq(:, is:ie, js:je), &
        & albedo_th, albedo_sol, &
        & lsolar, lthermal, &
        & edir, edn, eup, abso)
    end do

    call do_output()

    call specint_pprts_destroy(specint, pprts_solver, lfinalizepetsc=.true., ierr=ierr)

  contains

    subroutine do_output()
      character(len=default_str_len) :: groups(2), dimnames(4)
      logical :: lspectral_output, lflg
      real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf
      !real(ireals), allocatable, dimension(:, :, :, :) :: tmp_in, tmp_out ! tmp var to gather spectral output
      real(ireals), allocatable, dimension(:, :, :, :) :: gedir_s, gedn_s, geup_s, gabso_s ! global arrays for spectral output
      integer(iintegers) :: uid, k, Nspectral

      groups(1) = trim(outfile)

      associate (&
          & Cdir => pprts_solver%C_dir,    &
          & Cdiff => pprts_solver%C_diff,  &
          & C => pprts_solver%C_one,      &
          & C1 => pprts_solver%C_one1,    &
          & Ca => pprts_solver%C_one_atm, &
          & Ca1 => pprts_solver%C_one_atm1_box)

        dimnames(1) = 'zlev'
        dimnames(2) = 'nx'
        dimnames(3) = 'ny'
        if (allocated(edir)) then
          call gather_all_toZero(C1, edir, gedir)
          if (myid .eq. 0) then
            print *, 'dumping direct radiation with local and global shape', shape(edir), ':', shape(gedir)
            groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames); call CHKERR(ierr)
            call set_attribute(groups(1), 'edir', 'units', 'W/m2', ierr); call CHKERR(ierr)
          end if
        end if
        call gather_all_toZero(C1, edn, gedn)
        call gather_all_toZero(C1, eup, geup)
        call gather_all_toZero(C, abso, gabso)
        if (myid .eq. 0) then

          call set_global_attribute(groups(1), 'theta', pprts_solver%sun%theta, ierr); call CHKERR(ierr)
          call set_global_attribute(groups(1), 'phi', pprts_solver%sun%phi, ierr); call CHKERR(ierr)
          call set_global_attribute(groups(1), 'dx', pprts_solver%atm%dx, ierr); call CHKERR(ierr)
          call set_global_attribute(groups(1), 'dy', pprts_solver%atm%dy, ierr); call CHKERR(ierr)

          print *, 'dumping edn radiation with local and global shape', shape(edn), ':', shape(gedn)
          groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames); call CHKERR(ierr)
          call set_attribute(groups(1), 'edn', 'units', 'W/m2', ierr); call CHKERR(ierr)

          print *, 'dumping eup radiation with local and global shape', shape(eup), ':', shape(geup)
          groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames); call CHKERR(ierr)
          call set_attribute(groups(1), 'eup', 'units', 'W/m2', ierr); call CHKERR(ierr)

          print *, 'dumping abso radiation with local and global shape', shape(abso), ':', shape(gabso)
          dimnames(1) = 'zlay'
          groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames); call CHKERR(ierr)
          call set_attribute(groups(1), 'abso', 'units', 'W/m3', ierr); call CHKERR(ierr)

          print *, 'dumping z coords'
          call getVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
          dimnames(1) = 'nlev'
          groups(2) = 'zlev'; call ncwrite(groups, z(0, Ca1%zs:Ca1%ze, Ca1%xs, Ca1%ys), ierr, dimnames=dimnames(1:1))
          call CHKERR(ierr)
          call set_attribute(groups(1), 'zlev', 'units', 'm', ierr); call CHKERR(ierr)
          dimnames(1) = 'nlay'
          groups(2) = 'zlay'; call ncwrite(groups, &
 & (z(0, Ca1%zs:Ca1%ze - 1, Ca1%xs, Ca1%ys) + z(0, Ca1%zs + 1:Ca1%ze, Ca1%xs, Ca1%ys))*.5_ireals, &
 & ierr, dimnames=dimnames(1:1))
          call CHKERR(ierr)
          call set_attribute(groups(1), 'zlay', 'units', 'm', ierr); call CHKERR(ierr)
          call restoreVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
        end if

        lspectral_output = .false.
        call get_petsc_opt(PETSC_NULL_CHARACTER, '-spectral_output', lspectral_output, lflg, ierr); call CHKERR(ierr)

        if (lspectral_output) then
          Nspectral = count(pprts_solver%solutions(:)%lset)
          if (myid .eq. 0) then
            print *, 'found ', Nspectral, '3d fields'

            dimnames(1) = 'zlev'
            dimnames(2) = 'nx'
            dimnames(3) = 'ny'
            dimnames(4) = 'spectral'

            if (allocated(edir)) then
              allocate (gedir_s(Cdir%glob_zm, Cdir%glob_xm, Cdir%glob_ym, Nspectral))
            end if
            allocate (gedn_s(Cdiff%glob_zm, Cdiff%glob_xm, Cdiff%glob_ym, Nspectral))
            allocate (geup_s(Cdiff%glob_zm, Cdiff%glob_xm, Cdiff%glob_ym, Nspectral))
            allocate (gabso_s(C%glob_zm, C%glob_xm, C%glob_ym, Nspectral))
          end if

          k = 1
          do uid = lbound(pprts_solver%solutions, 1), ubound(pprts_solver%solutions, 1)
            if (pprts_solver%solutions(uid)%lset) then
              if (allocated(edir)) then
                call pprts_get_result_toZero(pprts_solver, gedn, geup, gabso, gedir=gedir, opt_solution_uid=uid)
                if (myid .eq. 0) gedir_s(:, :, :, k) = gedir
              else
                call pprts_get_result_toZero(pprts_solver, gedn, geup, gabso, opt_solution_uid=uid)
              end if
              if (myid .eq. 0) then
                gedn_s(:, :, :, k) = gedn
                geup_s(:, :, :, k) = geup
                gabso_s(:, :, :, k) = gabso
              end if
              k = k + 1
            end if
          end do

          if (myid .eq. 0) then
            if (allocated(gedir_s)) then
              groups(2) = 'edir_spectral'; call ncwrite(groups, gedir_s, ierr, dimnames=dimnames); call CHKERR(ierr)
              call set_attribute(groups(1), 'edir_spectral', 'units', 'W/m2/band', ierr); call CHKERR(ierr)
            end if
            groups(2) = 'edn_spectral'; call ncwrite(groups, gedn_s, ierr, dimnames=dimnames); call CHKERR(ierr)
            call set_attribute(groups(1), 'edn_spectral', 'units', 'W/m2/band', ierr); call CHKERR(ierr)
            groups(2) = 'eup_spectral'; call ncwrite(groups, geup_s, ierr, dimnames=dimnames); call CHKERR(ierr)
            call set_attribute(groups(1), 'eup_spectral', 'units', 'W/m2/band', ierr); call CHKERR(ierr)
            dimnames(1) = 'zlay'
            groups(2) = 'abso_spectral'; call ncwrite(groups, gabso_s, ierr, dimnames=dimnames); call CHKERR(ierr)
            call set_attribute(groups(1), 'abso_spectral', 'units', 'W/m3/band', ierr); call CHKERR(ierr)
          end if

        end if

      end associate

    end subroutine
  end subroutine

  subroutine load_input(comm, cldfile, dx, dy, is, ie, js, je, hhl, lwc, reliq, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    real(ireals), intent(out) :: dx, dy
    integer(iintegers), intent(out) :: is, ie, js, je
    real(ireals), allocatable, intent(out) :: hhl(:)
    real(ireals), allocatable, intent(out) :: lwc(:, :, :)
    real(ireals), allocatable, intent(out) :: reliq(:, :, :)
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    real(ireals), dimension(:, :, :), allocatable :: tmp ! used to resize input data
    character(len=default_str_len) :: groups(2)
    logical :: lflg

    call mpi_comm_rank(comm, myid, ierr)
    ! Load LibRadtran Cloud File
    if (myid .eq. 0) then
      call get_global_attribute(cldfile, 'dx', dx, ierr); call CHKERR(ierr)
      call get_global_attribute(cldfile, 'dy', dy, ierr); call CHKERR(ierr)
      groups(1) = trim(cldfile)
      groups(2) = trim('lwc'); call ncload(groups, lwc, ierr); call CHKERR(ierr)
      groups(2) = trim('reff'); call ncload(groups, reliq, ierr); call CHKERR(ierr)
      groups(2) = trim('z'); call ncload(groups, hhl, ierr); call CHKERR(ierr)

      is = lbound(lwc, dim=2); ie = ubound(lwc, dim=2)
      js = lbound(lwc, dim=3); je = ubound(lwc, dim=3)
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-xs', is, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-xe', ie, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-ys', js, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(PETSC_NULL_CHARACTER, '-ye', je, lflg, ierr); call CHKERR(ierr)
      allocate (tmp(size(lwc, dim=1), size(lwc, dim=2), size(lwc, dim=3)))
      tmp = lwc
      deallocate (lwc)
      allocate (lwc(size(tmp, dim=1), ie - is + 1, je - js + 1), source=tmp(:, is:ie, js:je))
      deallocate (tmp)
      allocate (tmp(size(reliq, dim=1), size(reliq, dim=2), size(reliq, dim=3)))
      tmp = reliq
      deallocate (reliq)
      allocate (reliq(size(tmp, dim=1), ie - is + 1, je - js + 1), source=tmp(:, is:ie, js:je))
      deallocate (tmp)
    end if
    call imp_bcast(comm, dx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dy, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, lwc, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, reliq, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, hhl, 0_mpiint, ierr); call CHKERR(ierr)

    if (size(lwc, dim=2) .eq. 1) call resize_arr(3_iintegers, lwc, dim=2, lrepeat=.true.)
    if (size(reliq, dim=2) .eq. 1) call resize_arr(3_iintegers, reliq, dim=2, lrepeat=.true.)

    if (size(lwc, dim=3) .eq. 1) call resize_arr(3_iintegers, lwc, dim=3, lrepeat=.true.)
    if (size(reliq, dim=3) .eq. 1) call resize_arr(3_iintegers, reliq, dim=3, lrepeat=.true.)

    dx = dx * 1e+3_ireals
    dy = dy * 1e+3_ireals
    hhl = hhl * 1e+3_ireals

    if (myid .eq. 0) then
      print *, 'Loaded LibRadtran Cloud File with:'
      print *, 'dx, dy:', dx, dy
      print *, 'hhl', hhl
      print *, 'shape lwc ', shape(lwc)
      print *, 'shape reliq', shape(reliq)
    end if
  end subroutine

  subroutine run_rrtmg_lw_sw(specint, pprts_solver, nxproc, nyproc, atm_filename, dx, dy, phi0, theta0, &
                             plev, tlev, lwc, reliq, albedo_th, albedo_sol, lsolar, lthermal, &
                             edir, edn, eup, abso)
    character(len=*), intent(in) :: specint                  ! name of module to use for spectral integration
    class(t_solver) :: pprts_solver
    integer(iintegers), intent(in) :: nxproc(:), nyproc(:)
    real(ireals), intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in), dimension(:, :, :), contiguous, target :: lwc, reliq ! dim(Nz,Nx,Ny)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lsolar, lthermal ! switches if solar or thermal computations should be done
    real(ireals), dimension(:, :, :), contiguous, target, intent(in) :: plev ! pressure on layer interfaces [hPa]   dim=nzp+1,nxp,nyp
    real(ireals), dimension(:, :, :), contiguous, target, intent(in) :: tlev ! Temperature on layer interfaces [K]  dim=nzp+1,nxp,nyp

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: comm, myid, ierr

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    ! real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=*), intent(in) :: atm_filename

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, plwc, preliq
    real(ireals) :: sundir(3)

    logical, parameter :: ldebug = .true.

    type(t_tenstr_atm) :: atm

    comm = mpi_comm_world
    call mpi_comm_rank(comm, myid, ierr)

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    ! h2ovmr = zero
    ! o3vmr  = zero
    ! co2vmr = zero
    ! ch4vmr = zero
    ! n2ovmr = zero
    ! o2vmr  = zero

    if (myid .eq. 0 .and. ldebug) print *, 'Setup Atmosphere...'

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    plwc(1:size(lwc, 1), 1:size(lwc, 2) * size(lwc, 3)) => lwc
    preliq(1:size(reliq, 1), 1:size(reliq, 2) * size(reliq, 3)) => reliq

    sundir = spherical_2_cartesian(phi0, theta0)

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_lwc=plwc, d_reliq=preliq)

    !if(myid.eq.0) call print_tenstr_atm(atm)

    call specint_pprts(specint, comm, pprts_solver, atm, &
                       size(plev, 2, kind=iintegers), size(plev, 3, kind=iintegers), &
                       dx, dy, sundir, &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          if (allocated(edir)) then
            print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
          else
            print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
          end if
        end do
      end if

      if (allocated(edir)) &
        print *, 'surface :: direct flux', edir(nlev, 1, 1)
      print *, 'surface :: downw flux ', edn(nlev, 1, 1)
      print *, 'surface :: upward fl  ', eup(nlev, 1, 1)
      print *, 'surface :: absorption ', abso(nlev - 1, 1, 1)

      if (allocated(edir)) &
        print *, 'TOA :: direct flux', edir(1, 1, 1)
      print *, 'TOA :: downw flux ', edn(1, 1, 1)
      print *, 'TOA :: upward fl  ', eup(1, 1, 1)
      print *, 'TOA :: absorption ', abso(1, 1, 1)

    end if

    ! Tidy up
    call destroy_tenstr_atm(atm)
  end subroutine

  subroutine example_uvspec_cld_file_with_plexrt(comm, &
                                                 cldfile, atm_filename, outfile, &
                                                 albedo_th, albedo_sol, &
                                                 lsolar, lthermal, &
                                                 phi0, theta0, &
                                                 Tsrfc, dTdz)

    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: Tsrfc, dTdz

    real(ireals), dimension(:, :, :), allocatable, target :: glob_lwc, glob_reliq ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:, :), allocatable, target :: plev, tlev ! will have local shape nzp+1, Ncol
    real(ireals), dimension(:), allocatable :: hhl ! dim Nz+1
    character(len=default_str_len) :: groups(2)

    real(ireals), allocatable, dimension(:, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    !real(ireals),allocatable, dimension(:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    real(ireals) :: dx, dy
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k

    integer(iintegers) :: Nx_global, Ny_global, Nz
    integer(iintegers) :: fStart, fEnd, Ncol, Nlev

    type(tDM) :: dm2d, dm2d_dist, dm3d
    type(tPetscSF) :: migration_sf
    type(tPetscSection) :: parCellSection

    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    class(t_plex_solver), allocatable :: solver
    type(t_tenstr_atm) :: atm
    real(ireals), pointer :: xlwc(:), xreff(:), col_lwc(:, :), col_reff(:, :)
    type(tVec) :: vlwc, vreff

    real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

    call mpi_comm_rank(comm, myid, ierr)

    if (myid .eq. 0) then
      ! Load LibRadtran Cloud File
      call get_global_attribute(cldfile, 'dx', dx, ierr); call CHKERR(ierr)
      call get_global_attribute(cldfile, 'dy', dy, ierr); call CHKERR(ierr)
      groups(1) = trim(cldfile)
      groups(2) = trim('lwc'); call ncload(groups, glob_lwc, ierr); call CHKERR(ierr)
      groups(2) = trim('reff'); call ncload(groups, glob_reliq, ierr); call CHKERR(ierr)
      groups(2) = trim('z'); call ncload(groups, hhl, ierr); call CHKERR(ierr)

      if (size(glob_lwc, dim=2) .eq. 1) call resize_arr(3_iintegers, glob_lwc, dim=2, lrepeat=.true.)
      print *, 'shape lwc', shape(glob_lwc)
      if (size(glob_reliq, dim=2) .eq. 1) call resize_arr(3_iintegers, glob_reliq, dim=2, lrepeat=.true.)
      print *, 'shape lwc', shape(glob_lwc)

      if (size(glob_lwc, dim=3) .eq. 1) call resize_arr(3_iintegers, glob_lwc, dim=3, lrepeat=.true.)
      print *, 'shape lwc', shape(glob_lwc)
      if (size(glob_reliq, dim=3) .eq. 1) call resize_arr(3_iintegers, glob_reliq, dim=3, lrepeat=.true.)
      print *, 'shape lwc', shape(glob_lwc)

      dx = dx * 1e+3_ireals
      dy = dy * 1e+3_ireals
      hhl = hhl * 1e+3_ireals

      Nz = size(glob_lwc, dim=1)
      Nx_global = size(glob_lwc, dim=2)
      Ny_global = size(glob_lwc, dim=3)

      print *, 'Loaded LibRadtran Cloud File with:'
      print *, 'dx, dy:', dx, dy
      print *, 'global Nx/Ny', Nx_global, Ny_global
      print *, 'hhl', hhl
      print *, 'shape lwc ', shape(glob_lwc)
      print *, 'shape reliq', shape(glob_reliq)
    end if

    call imp_bcast(comm, Nx_global, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Ny_global, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Nz, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dy, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, hhl, 0_mpiint, ierr); call CHKERR(ierr)

    call create_2d_regular_plex(comm, Nx_global + 1, Ny_global + 1, dm2d, dm2d_dist, &
                                opt_migration_sf=migration_sf, opt_dx=dx)

    call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart

    !if(myid.eq.0) then
    print *, myid, 'Local Domain sizes are:', Ncol, Nz
    !endif

    ! Start with a dynamics grid starting at 1000 hPa with a specified lapse rate and surface temperature
    allocate (plev(Nz + 1, fStart:fEnd - 1), tlev(Nz + 1, fStart:fEnd - 1))
    plev(1, :) = 1000_ireals
    tlev(1, :) = Tsrfc
    do k = 2, Nz + 1
      tlev(k, :) = tlev(k - 1, :) + dTdz * (hhl(k) - hhl(k - 1))
      plev(k, :) = plev(k - 1, :) - hydrostat_dp(hhl(k) - hhl(k - 1), plev(k - 1, :), (tlev(k, :) + tlev(k - 1, :)) / 2)
    end do

    if (myid .eq. 0) then
      do k = 1, Nz + 1
        print *, k, 'plev', plev(k, fStart), 'Tlev', tlev(k, fStart)
      end do
    end if
    hhl = reverse(hhl)

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          plev, tlev, atm)

    Nlev = size(atm%plev, 1, kind=iintegers)
    call dmplex_2D_to_3D(dm2d_dist, Nlev, reverse(atm%zt(:, i1)), [zero, zero, -huge(one)], dm3d, zindex)

    call setup_plexgrid(dm2d_dist, dm3d, Nlev - 1, zindex, plex, reverse(atm%zt(:, i1)))
    deallocate (zindex)

    if (myid .eq. 0) then
      ! increase the size by two in the x direction,
      ! i.e. from boxes to wedges with 2 elements per box
      call resize_arr(size(glob_lwc, 2, kind=iintegers) * 2, glob_lwc, dim=2, fillVal=-1._ireals)
      call resize_arr(size(glob_reliq, 2, kind=iintegers) * 2, glob_reliq, dim=2, fillVal=-1._ireals)
      glob_lwc(:, 1:size(glob_lwc, dim=2):2, :) = glob_lwc(:, 1:size(glob_lwc, dim=2) / 2, :)
      glob_lwc(:, 2:size(glob_lwc, dim=2):2, :) = glob_lwc(:, 1:size(glob_lwc, dim=2):2, :)
      glob_reliq(:, 1:size(glob_reliq, dim=2):2, :) = glob_reliq(:, 1:size(glob_reliq, dim=2) / 2, :)
      glob_reliq(:, 2:size(glob_reliq, dim=2):2, :) = glob_reliq(:, 1:size(glob_reliq, dim=2):2, :)

      col_lwc(1:size(glob_lwc, 1), 1:size(glob_lwc, 2) * size(glob_lwc, 3)) => glob_lwc
      col_reff(1:size(glob_lwc, 1), 1:size(glob_lwc, 2) * size(glob_lwc, 3)) => glob_reliq
      print *, myid, 'global shape col lwc', shape(col_lwc)
    end if

    call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, col_lwc, parCellSection, vlwc)
    call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, col_reff, parCellSection, vreff)
    nullify (col_lwc)
    nullify (col_reff)

    call VecGetArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(vreff, xreff, ierr); call CHKERR(ierr)

    col_lwc(1:Nz, fstart:fEnd - 1) => xlwc
    col_reff(1:Nz, fstart:fEnd - 1) => xreff

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          plev, tlev, atm, &
                          d_lwc=col_lwc, d_reliq=col_reff)

    nullify (col_lwc)
    nullify (col_reff)
    call VecRestoreArrayReadF90(vlwc, xlwc, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(vreff, xreff, ierr); call CHKERR(ierr)

    ! Finished preparing input, lets do the computations

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)

    sundir = spherical_2_cartesian(phi0, theta0)
    print *, 'sundir', sundir

    call plexrt_rrtmg(solver, atm, sundir, &
                      albedo_thermal=albedo_th, albedo_solar=albedo_sol, &
                      lthermal=lthermal, lsolar=lsolar, &
                      edir=edir, edn=edn, eup=eup, abso=abso)

    call dump_results()

    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

    call destroy_plexrt_rrtmg(solver, lfinalizepetsc=.false.)
  contains
    subroutine transfer_arr(col_inp, arr3d, idof_start, idof_end)
      real(ireals), intent(in) :: col_inp(:, :)
      real(ireals), intent(out) :: arr3d(:, :, :)
      integer(iintegers), intent(in) :: idof_start, idof_end
      integer(iintegers) :: i, j, icol, idof
      do j = 1, Ny_global
        do i = 1, Nx_global
          icol = Nx_global * (j - 1) + i
          arr3d(:, i, j) = 0
          do idof = idof_start, idof_end
            arr3d(:, i, j) = arr3d(:, i, j) + col_inp(:, 2 * (icol - 1) + idof)
          end do
          arr3d(:, i, j) = arr3d(:, i, j) / real(idof_end - idof_start + 1, ireals)
        end do
      end do
    end subroutine
    subroutine dump_var(var, varname)
      real(ireals), allocatable, intent(in) :: var(:, :)
      character(len=*), intent(in) :: varname
      type(tPetscSection) :: flxSection, r0flxSection
      type(tVec) :: v_var
      type(tVec) :: r0var

      real(ireals), pointer :: xarr(:), xxarr(:, :)
      real(ireals), dimension(:, :, :), allocatable :: oarr

      groups(1) = trim(outfile)
      if (myid .eq. 0) allocate (oarr(size(var, dim=1), Nx_global, Ny_global))

      if (.not. allocated(var)) return

      call create_plex_section(dm2d_dist, 'face_section', i1, &
                               [i0], [size(var, dim=1, kind=iintegers)], [i0], [i0], flxSection)

      call dmplex_gVec_from_f90_array(comm, var, v_var)
      call plex_gVec_toZero(dm2d_dist, migration_sf, flxSection, v_var, &
                            r0flxSection, r0var)

      call VecGetArrayF90(r0var, xarr, ierr); call CHKERR(ierr)
      if (myid .eq. 0) then
        xxarr(1:size(var, dim=1), 1:Nx_global * Ny_global * 2) => xarr

        call transfer_arr(xxarr, oarr, i1, i1)
        groups(2) = trim(varname)//'_leftvertices'; call ncwrite(groups, oarr, ierr); call CHKERR(ierr)
        call transfer_arr(xxarr, oarr, i2, i2)
        groups(2) = trim(varname)//'_rightvertices'; call ncwrite(groups, oarr, ierr); call chkerr(ierr)
        call transfer_arr(xxarr, oarr, i1, i2)
        groups(2) = trim(varname); call ncwrite(groups, oarr, ierr); call chkerr(ierr)

        nullify (xxarr)
      end if
      call VecRestoreArrayF90(r0var, xarr, ierr); call CHKERR(ierr)

      call VecDestroy(v_var, ierr); call CHKERR(ierr)
      call VecDestroy(r0var, ierr); call CHKERR(ierr)
    end subroutine
    subroutine dump_results()

      call dump_var(edir, 'edir')
      call dump_var(edn, 'edn')
      call dump_var(eup, 'eup')
      call dump_var(abso, 'abso')

      if (myid .eq. 0) then
        groups(1) = trim(outfile)
        groups(2) = 'lwc'; call ncwrite(groups, glob_lwc, ierr); call CHKERR(ierr)
        groups(2) = 'reff'; call ncwrite(groups, glob_reliq, ierr); call CHKERR(ierr)
      end if
    end subroutine
  end subroutine

end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only: mpiint, share_dir
  use m_helper_functions, only: get_petsc_opt
  use m_example_uvspec_cld_file

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg
  character(len=10*default_str_len) :: cldfile, outfile
  real(ireals) :: Ag, phi0, theta0, Tsrfc, dTdz
  character(len=default_str_len) :: atm_filename, specint
  integer(iintegers) :: scene_shift_x, scene_shift_y, scene_shift_it

  logical :: luse_plexrt

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  specint = 'no default set'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-specint', specint, lflg, ierr); call CHKERR(ierr)

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need to supply a cloud filename... please call with -cld <libRadtran_cloud_file.nc>')

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  atm_filename = share_dir//'tenstream_default.atm'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  Ag = .1
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-solar", lsolar, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)

  scene_shift_x = 0
  scene_shift_y = 0
  scene_shift_it = 1
  call get_petsc_opt(PETSC_NULL_CHARACTER, &
    & '-scene_shift_x', scene_shift_x, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, &
    & '-scene_shift_y', scene_shift_y, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, &
    & '-scene_shift_it', scene_shift_it, lflg, ierr); call CHKERR(ierr)

  phi0 = 270
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 60
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  Tsrfc = 288
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Tsrfc", Tsrfc, lflg, ierr); call CHKERR(ierr)
  dTdz = -6.5_ireals * 1e-3_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dTdz", dTdz, lflg, ierr); call CHKERR(ierr)

  luse_plexrt = .false.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-use_plexrt", luse_plexrt, lflg, ierr); call CHKERR(ierr)

  if (luse_plexrt) then
    call CHKERR(int(scene_shift_x, mpiint), 'not supported option')
    call CHKERR(int(scene_shift_y, mpiint), 'not supported option')
    call CHKERR(int(scene_shift_it - 1, mpiint), 'not supported option')
    call example_uvspec_cld_file_with_plexrt(&
      & mpi_comm_world, cldfile, atm_filename, outfile, &
      & zero, Ag, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)
  else
    call example_uvspec_cld_file_pprts(&
      & specint, &
      & mpi_comm_world, cldfile, atm_filename, outfile, &
      & zero, Ag, lsolar, lthermal, phi0, theta0, &
      & scene_shift_x, scene_shift_y, scene_shift_it)
  end if
  call mpi_finalize(ierr)
end program
