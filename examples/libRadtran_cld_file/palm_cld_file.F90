module m_example_palm_cld_file
  !
  ! Read PALM (LES) 3d output and compute 3d radiative transfer (solar + thermal)
  ! with the TenStream / specint machinery for selected timesteps.
  !
  ! This is the PALM sibling of `uclales_cld_file.F90`. The main differences to PALM:
  !
  !   * PALM stores potential temperature (`theta` [K], ref 1000 hPa), not T, and the
  !     pressure field `p` is only the *perturbation* pressure of the pressure solver
  !     -- it is NOT the absolute pressure. We therefore reconstruct a hydrostatic
  !     pressure profile per column from a prescribed surface pressure (`-psrfc`,
  !     in hPa) and the temperature profile (fixed point iteration).
  !
  !   * PALM output uses a Cartesian grid with terrain/buildings masking. Grid cells
  !     below the topography are flagged with the netcdf `_FillValue`. We turn every
  !     solid cell (terrain + buildings, taken from `zusi` = zu(nzb_s_inner)) into a
  !     TenStream "building" box, exactly like `ex_pprts_specint_buildings_from_file`.
  !
  !   * BEWARE OF UNITS. The LES netcdf attributes are frequently wrong, e.g. a field
  !     labelled `kg/kg` is in fact `g/kg`, or `degree_` is actually degree Celsius.
  !     All unit conversions are therefore exposed as command line options
  !     (`-qv_scale`, `-lwc_scale`, ...) so you can override them without recompiling.
  !
  ! Example call:
  !
  !   mpirun -np 8 ./ex_palm_cld_file \
  !     -specint repwvl \
  !     -cld /path/to/..._3d.000.nc \
  !     -out palm_rad.nc \
  !     -atm afglus_100m.dat \
  !     -tstart 1 -tend 7 -tinc 1 \
  !     -phi 180 -theta 47 \
  !     -Ag 0.15 \
  !     -ix0 400 -ix1 623 -iy0 300 -iy1 523 -ztop 4000
  !
#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#endif
  use mpi
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, &
                               i0, i1, zero, one, default_str_len, share_dir, &
                               R_DRY_AIR, CP_DRY_AIR

  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  use m_tenstr_atm, only: &
    & destroy_tenstr_atm, &
    & print_tenstr_atm, &
    & reff_from_lwc_and_N, &
    & setup_tenstr_atm, &
    & hydrostat_plev, &
    & t_tenstr_atm

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & domain_decompose_2d_petsc, &
    & deg2rad, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_sum, &
    & imp_bcast, &
    & imp_scan_sum, &
    & ind_1d_to_nd, &
    & is_inrange, &
    & meanval, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_netcdfio, only: ncload, ncwrite, set_attribute, get_global_attribute, list_global_attributes

  use m_buildings, only: &
    & t_pprts_buildings, &
    & init_buildings, &
    & clone_buildings, &
    & check_buildings_consistency, &
    & faceidx_by_cell_plus_offset

  use m_boxmc_geometry, only: PPRTS_TOP_FACE

  implicit none

  private
  public :: example_palm_cld_file_pprts

  ! molar mass ratio dry air / water vapour, to convert specific humidity to volume mixing ratio
  real(ireals), parameter :: MDRY_O_MH2O = 1.60771_ireals

  ! everything below this value in the LES netcdf is considered a fill value / masked (solid) point
  real(ireals), parameter :: FILL_THRESHOLD = -1.e5_ireals

  ! _FillValue written into the heating-rate field inside solid (terrain/building) cells
  real(ireals), parameter :: HR_FILL = -9999._ireals

contains

  ! (liquid water) potential temperature to temperature, exner from pressure p (same unit as p0)
  elemental function Tpot2T(Tpot, p, p0) result(T)
    real(ireals), intent(in) :: Tpot, p, p0
    real(ireals) :: T
    real(ireals), parameter :: kappa = R_DRY_AIR / CP_DRY_AIR
    T = Tpot * (p / p0)**kappa
  end function

  !> Read grid/coordinate meta data (rank 0 reads, then broadcast)
  subroutine load_meta_data(comm, cldfile, time, zw_full, zu_full, Nx_full, Ny_full, dx, dy, origin_z, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    real(ireals), allocatable, intent(out) :: time(:), zw_full(:), zu_full(:)
    integer(iintegers), intent(out) :: Nx_full, Ny_full
    real(ireals), intent(out) :: dx, dy, origin_z
    integer(mpiint), intent(out) :: ierr

    logical :: lflg
    integer(mpiint) :: myid, ierr2
    real(ireals), allocatable :: dimx(:), dimy(:)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      call list_global_attributes(cldfile, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, 'time'], time, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, 'zw_3d'], zw_full, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, 'zu_3d'], zu_full, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, 'x'], dimx, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, 'y'], dimy, ierr); call CHKERR(ierr)
      Nx_full = size(dimx)
      Ny_full = size(dimy)
      dx = dimx(2) - dimx(1)
      dy = dimy(2) - dimy(1)
      call get_petsc_opt('', "-dx", dx, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt('', "-dy", dy, lflg, ierr); call CHKERR(ierr)

      origin_z = zero
      call get_global_attribute(cldfile, 'origin_z', origin_z, ierr2) ! optional, ignore errors
      call get_petsc_opt('', "-surface_height", origin_z, lflg, ierr); call CHKERR(ierr)

      print *, 'PALM file: ', trim(cldfile)
      print *, '  timesteps        :', size(time)
      print *, '  Nx, Ny           :', size(dimx), size(dimy)
      print *, '  dx, dy           :', dx, dy
      print *, '  Nz (zw_3d/zu_3d) :', size(zw_full), size(zu_full)
      print *, '  z range          :', zw_full(1), '..', zw_full(size(zw_full))
      print *, '  surface_height   :', origin_z, ' [m ASL] (used as flat base for the RT grid)'
    end if
    call imp_bcast(comm, time, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, zw_full, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, zu_full, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Nx_full, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Ny_full, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dy, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, origin_z, 0_mpiint, ierr); call CHKERR(ierr)
    ierr = 0
  end subroutine

  !> Read the topography height field zusi(x,y) [m] for the cropped region and derive
  !> the number of solid (terrain/building) layers per column, kterr(x,y)
  subroutine load_topography(comm, cldfile, ix0, iy0, Nx, Ny, zc, zusi, kterr, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    integer(iintegers), intent(in) :: ix0, iy0, Nx, Ny
    real(ireals), intent(in) :: zc(:) ! layer center heights [m], dim(nlay)
    real(ireals), allocatable, intent(out) :: zusi(:, :) ! dim(Nx, Ny)
    integer(iintegers), allocatable, intent(out) :: kterr(:, :) ! dim(Nx, Ny)
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer(iintegers) :: i, j, k, nlay, kmax
    real(ireals), allocatable :: tmp2d(:, :)

    nlay = size(zc)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    allocate (kterr(Nx, Ny))
    if (myid .eq. 0) then
      call ncload([character(len=default_str_len) :: cldfile, 'zusi'], tmp2d, ierr, &
        & ostart=[integer :: int(ix0), int(iy0)], ocount=[integer :: int(Nx), int(Ny)]); call CHKERR(ierr)
      allocate (zusi(Nx, Ny), source=tmp2d)

      kmax = 0
      do j = 1, Ny
        do i = 1, Nx
          k = 0
          do while (k .lt. nlay)
            if (zc(k + 1) .gt. zusi(i, j)) exit
            k = k + 1
          end do
          kterr(i, j) = k
          kmax = max(kmax, k)
        end do
      end do
      print *, 'Topography: zusi min/max      :', minval(zusi), maxval(zusi)
      print *, 'Topography: solid layers min/max:', minval(kterr), maxval(kterr), ' (of nlay=', nlay, ')'
      if (kmax .gt. nlay - 2) call CHKERR(1_mpiint, &
        & 'Topography reaches (almost) the model top -- increase -ztop so that some clear layers remain above the terrain')
    end if
    call imp_bcast(comm, zusi, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, kterr, 0_mpiint, ierr); call CHKERR(ierr)
    ierr = 0
  end subroutine

  !> Read one timestep of the LES cloud/thermodynamic fields for THIS RANK's subdomain
  !> and convert them to what setup_tenstr_atm expects:
  !>   plev [hPa] (levels), tlev [K] (levels), tlay [K] (layers),
  !>   h2ovmr [vol mix ratio] (layers), lwc [g/kg] (layers), reliq [micron] (layers)
  !> Each rank opens the file itself and reads only its own hyperslab -- nothing is
  !> broadcast, so memory stays O(local subdomain) even for large PALM domains.
  !> gx0/gy0 are the 1-based global netcdf x/y start of this subdomain; is0/js0 the
  !> 1-based crop-relative start (to index into the global kterr field).
  subroutine load_timestep_data(comm, cldfile, it, gx0, gy0, is0, js0, nxp, nyp, &
    & zlev, zc, kterr, psrfc_hPa, qv_scale, lwc_scale, cldN, &
    & plev, tlev, tlay, h2ovmr, lwc, reliq, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    integer(iintegers), intent(in) :: it, gx0, gy0, is0, js0, nxp, nyp
    real(ireals), intent(in) :: zlev(:), zc(:) ! [m], dim(nlev), dim(nlay)
    integer(iintegers), intent(in) :: kterr(:, :) ! global (crop) field, dim(Nx, Ny)
    real(ireals), intent(in) :: psrfc_hPa, qv_scale, lwc_scale, cldN
    real(ireals), allocatable, intent(inout) :: plev(:, :, :), tlev(:, :, :) ! dim(nlev, nxp, nyp)
    real(ireals), allocatable, intent(inout) :: tlay(:, :, :), h2ovmr(:, :, :), lwc(:, :, :), reliq(:, :, :) ! dim(nlay, nxp, nyp)
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer(iintegers) :: i, j, k, nlay, nlev, kt, iter
    integer :: ostart(4), ocount(4)

    real(ireals), allocatable :: theta(:, :, :, :), qv(:, :, :, :), ql(:, :, :, :) ! (nxp,nyp,nlay,1)
    real(ireals) :: theta_col(size(zc)), tlay_col(size(zc)), tlev_col(size(zlev))
    real(ireals) :: plev_pa(size(zlev)), dp_pa(size(zc))
    real(ireals) :: q_col(size(zc)), ql_col(size(zc))
    real(ireals) :: play_pa, rho_air, lwc_gm3, w, qq
    real(ireals), parameter :: p0_pa = 1.e5_ireals
    integer(iintegers), parameter :: niter = 15

    nlay = size(zc)
    nlev = size(zlev)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ! PALM 3d fields are dimensioned (time, z, y, x) in CDL, i.e. (x, y, z, time) in Fortran.
    ! Scalar level 1 (zu_3d==0) is the surface boundary point, so our layer kk maps to level kk+1.
    ostart = [integer :: int(gx0), int(gy0), 2, int(it)]
    ocount = [integer :: int(nxp), int(nyp), int(nlay), 1]
    call ncload([character(len=default_str_len) :: cldfile, 'theta'], theta, ierr, ostart=ostart, ocount=ocount)
    call CHKERR(ierr)
    call ncload([character(len=default_str_len) :: cldfile, 'qv'], qv, ierr, ostart=ostart, ocount=ocount); call CHKERR(ierr)
    call ncload([character(len=default_str_len) :: cldfile, 'ql'], ql, ierr, ostart=ostart, ocount=ocount); call CHKERR(ierr)

    do j = 1, nyp
      do i = 1, nxp
        kt = kterr(is0 + i - 1, js0 + j - 1)

        ! --- gather the raw column, replacing masked / solid points by the nearest valid value above ---
        do k = 1, nlay
          theta_col(k) = theta(i, j, k, 1)
          q_col(k) = qv(i, j, k, 1) * qv_scale     ! -> kg/kg (qv_scale fixes a g/kg vs kg/kg mislabel)
          ql_col(k) = ql(i, j, k, 1) * lwc_scale   ! -> g/kg
        end do
        do k = nlay - 1, 1, -1 ! fill masked points from above
          if (theta_col(k) .lt. FILL_THRESHOLD) theta_col(k) = theta_col(k + 1)
          if (q_col(k) .lt. FILL_THRESHOLD .or. q_col(k) .lt. zero) q_col(k) = q_col(k + 1)
          if (ql_col(k) .lt. FILL_THRESHOLD) ql_col(k) = ql_col(k + 1)
        end do
        do k = 1, kt ! solid cells: no cloud, take thermodynamics from first clear layer above
          theta_col(k) = theta_col(kt + 1)
          q_col(k) = q_col(kt + 1)
          ql_col(k) = zero
        end do
        ql_col = max(ql_col, zero)

        ! --- hydrostatic pressure from surface pressure + fixed-point on T(theta, p) ---
        tlay_col = theta_col ! first guess
        do iter = 1, niter
          call hydrostat_plev(psrfc_hPa * 1.e2_ireals, tlay_col, zlev, plev_pa, dp_pa)
          do k = 1, nlay
            play_pa = 0.5_ireals * (plev_pa(k) + plev_pa(k + 1))
            tlay_col(k) = Tpot2T(theta_col(k), play_pa, p0_pa)
          end do
        end do

        ! --- temperature on levels (linear in height, extrapolate the two ends) ---
        do k = 2, nlay
          w = (zlev(k) - zc(k - 1)) / (zc(k) - zc(k - 1))
          tlev_col(k) = tlay_col(k - 1) * (one - w) + tlay_col(k) * w
        end do
        tlev_col(1) = tlay_col(1) + (tlay_col(1) - tlev_col(2))
        tlev_col(nlev) = tlay_col(nlay) + (tlay_col(nlay) - tlev_col(nlay - 1))

        plev(:, i, j) = plev_pa * 1.e-2_ireals ! Pa -> hPa
        tlev(:, i, j) = tlev_col
        tlay(:, i, j) = tlay_col

        do k = 1, nlay
          qq = min(max(q_col(k), zero), 0.999_ireals)
          h2ovmr(k, i, j) = max(MDRY_O_MH2O * qq / (one - qq), 1.e-8_ireals)

          lwc(k, i, j) = ql_col(k) ! [g/kg]

          play_pa = 0.5_ireals * (plev_pa(k) + plev_pa(k + 1))
          rho_air = play_pa / (R_DRY_AIR * tlay_col(k)) ! [kg/m3]
          lwc_gm3 = ql_col(k) * rho_air                  ! [g/kg]*[kg/m3] = [g/m3]
          if (lwc_gm3 .gt. 1.e-9_ireals) then
            reliq(k, i, j) = min(max(reff_from_lwc_and_N(lwc_gm3, cldN), 2.5_ireals), 60._ireals)
          else
            reliq(k, i, j) = 10._ireals
          end if
        end do
      end do
    end do

    if (myid .eq. 0) then
      print *, 'timestep', it, ' (rank 0 subdomain) min/max plev  ', minval(plev), maxval(plev)
      print *, 'timestep', it, ' (rank 0 subdomain) min/max tlev  ', minval(tlev), maxval(tlev)
      print *, 'timestep', it, ' (rank 0 subdomain) min/max h2ovmr', minval(h2ovmr), maxval(h2ovmr)
      print *, 'timestep', it, ' (rank 0 subdomain) min/max lwc   ', minval(lwc), maxval(lwc)
      print *, 'timestep', it, ' (rank 0 subdomain) min/max reliq ', minval(reliq), maxval(reliq)
    end if
    ierr = 0
  end subroutine

  !> Effective broadband shortwave surface albedo per column, taken from PALM`s own
  !> radiation output: alb = rad_sw_out / rad_sw_in at the lowest illuminated (non-fill)
  !> level of each column. PALM has no explicit albedo variable -- the _static driver
  !> only classifies surface types (vegetation/pavement/soil/water) from which PALM
  !> derives albedo internally. Columns without sunlight fall back to ag_default.
  subroutine load_surface_albedo(cldfile, it, gx0, gy0, nxp, nyp, nz_full, ag_default, alb2d, ierr)
    character(len=*), intent(in) :: cldfile
    integer(iintegers), intent(in) :: it, gx0, gy0, nxp, nyp, nz_full
    real(ireals), intent(in) :: ag_default
    real(ireals), pointer, intent(inout) :: alb2d(:, :) ! (nxp,nyp), must be associated
    integer(mpiint), intent(out) :: ierr

    real(ireals), allocatable :: swin(:, :, :, :), swout(:, :, :, :) ! (nxp,nyp,nz_full,1)
    integer :: ost(4), oc(4)
    integer(iintegers) :: i, j, k
    real(ireals) :: si, so

    ! rad_sw_* are on zw_3d (level 1 == z=0); dims (time,z,y,x) -> Fortran (x,y,z,time)
    ost = [integer :: int(gx0), int(gy0), 1, int(it)]
    oc = [integer :: int(nxp), int(nyp), int(nz_full), 1]
    call ncload([character(len=default_str_len) :: cldfile, 'rad_sw_in'], swin, ierr, ostart=ost, ocount=oc); call CHKERR(ierr)
    call ncload([character(len=default_str_len) :: cldfile, 'rad_sw_out'], swout, ierr, ostart=ost, ocount=oc); call CHKERR(ierr)

    do j = 1, nyp
      do i = 1, nxp
        alb2d(i, j) = ag_default
        do k = 1, nz_full ! bottom-up: first non-fill, sunlit level == the surface
          si = swin(i, j, k, 1)
          so = swout(i, j, k, 1)
          if (si .gt. 1._ireals .and. si .lt. 1.e6_ireals) then
            alb2d(i, j) = min(max(so / si, 0.02_ireals), 0.9_ireals)
            exit
          end if
        end do
      end do
    end do
    ierr = 0
  end subroutine

  subroutine run_lw_sw(&
      & specint, comm, pprts_solver, atm, atm_filename, &
      & nxproc, nyproc, dx, dy, phi0, theta0, &
      & albedo_th, albedo_sol, lsolar, lthermal, time, &
      & Nx_glob, Ny_glob, kterr, &
      & building_albedo_sol, building_albedo_th, building_temp, &
      & solar_albedo_2d, &
      & plev, tlev, tlay, h2ovmr, lwc, reliq, &
      & buildings_solar, buildings_thermal, &
      & edir, edn, eup, abso, hr)
    character(len=*), intent(in) :: specint
    integer(mpiint), intent(in) :: comm
    class(t_solver) :: pprts_solver
    type(t_tenstr_atm), intent(inout) :: atm
    character(len=*), intent(in) :: atm_filename
    integer(iintegers), intent(in) :: nxproc(:), nyproc(:)
    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: time
    integer(iintegers), intent(in) :: Nx_glob, Ny_glob
    integer(iintegers), intent(in) :: kterr(:, :) ! solid layers per column, dim(Nx_glob, Ny_glob)
    real(ireals), intent(in) :: building_albedo_sol, building_albedo_th, building_temp
    real(ireals), pointer, intent(in) :: solar_albedo_2d(:, :) ! (nxp,nyp) per-column SW surface albedo, or null() to use the scalar albedo_sol

    real(ireals), intent(in), dimension(:, :, :), contiguous, target :: plev, tlev ! (nlev, nxp, nyp)
    real(ireals), intent(in), dimension(:, :, :), contiguous, target :: tlay, h2ovmr, lwc, reliq ! (nlay, nxp, nyp)

    type(t_pprts_buildings), allocatable, intent(inout) :: buildings_solar, buildings_thermal

    real(ireals), allocatable, dimension(:, :, :), intent(out) :: edir, edn, eup, abso ! [W/m2], [W/m3]
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: hr ! heating rate [K/day]

    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k, nlev, nxp, nyp
    integer(iintegers) :: i, j, kk, ka, kt, icol, nlay_m
    real(ireals) :: play_pa
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, ptlay, ph2ovmr, plwc, preliq
    real(ireals) :: sundir(3)
    real(ireals) :: timeofday, tod_offset, tod_theta, tod_phi
    logical :: lflg
    logical, parameter :: ldebug = .true.

    call mpi_comm_rank(comm, myid, ierr)

    nxp = size(plev, 2)
    nyp = size(plev, 3)
    nlev = size(plev, 1)

    pplev(1:size(plev, 1), 1:nxp * nyp) => plev
    ptlev(1:size(tlev, 1), 1:nxp * nyp) => tlev
    ptlay(1:size(tlay, 1), 1:nxp * nyp) => tlay
    ph2ovmr(1:size(h2ovmr, 1), 1:nxp * nyp) => h2ovmr
    plwc(1:size(lwc, 1), 1:nxp * nyp) => lwc
    preliq(1:size(reliq, 1), 1:nxp * nyp) => reliq

    sundir = spherical_2_cartesian(phi0, theta0)

    ! optional simple diurnal cycle, identical to the uclales example: -tod_offset gives
    ! the time of day (in fractions of a day) at model time == 0.
    tod_offset = 0
    call get_petsc_opt('', "-tod_offset", tod_offset, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      timeofday = modulo(time / 86400._ireals + tod_offset, 1._ireals)
      tod_phi = modulo(timeofday * 360._ireals, 360._ireals)
      tod_theta = 1._ireals - max(0._ireals, sin(deg2rad(timeofday * 360._ireals - 90)))
      sundir = spherical_2_cartesian(tod_phi, tod_theta * 90)
      if (myid .eq. 0) print *, 'diurnal cycle: phi0', tod_phi, 'theta0', tod_theta * 90, 'sundir', sundir
    end if

    if (myid .eq. 0 .and. ldebug) print *, 'Setup Atmosphere...'
    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_tlay=ptlay, &
                          d_h2ovmr=ph2ovmr, &
                          d_lwc=plwc, d_reliq=preliq)

    ! first call in the run: only build the grid structures, then create the building boxes
    if (.not. allocated(buildings_solar)) then
      call specint_pprts(specint, comm, pprts_solver, atm, nxp, nyp, dx, dy, sundir, &
                         albedo_th, albedo_sol, lthermal, lsolar, edir, edn, eup, abso, &
                         nxproc=nxproc, nyproc=nyproc, lonly_initialize=.true.)
      call setup_buildings()
    end if

    ! (re)set each building face temperature to the air temperature of the nearest air
    ! voxel (the first cell above the terrain top of that column); updated every timestep.
    if (allocated(buildings_thermal)) call set_building_temps_from_air()

    call specint_pprts(specint, comm, pprts_solver, atm, nxp, nyp, dx, dy, sundir, &
                       albedo_th, albedo_sol, lthermal, lsolar, edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=time, &
                       solar_albedo_2d=solar_albedo_2d, &
                       opt_buildings_solar=buildings_solar, &
                       opt_buildings_thermal=buildings_thermal)

    if (myid .eq. 0 .and. ldebug) then
      do k = 1, ubound(edn, 1)
        if (allocated(edir)) then
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1)
        else
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1)
        end if
      end do
      k = ubound(edn, 1)
      if (allocated(edir)) print *, 'surface :: direct flux', edir(k, 1, 1)
      print *, 'surface :: downw flux ', edn(k, 1, 1)
      print *, 'surface :: upward flux', eup(k, 1, 1)
    end if

    ! heating rate:  dT/dt [K/day] = abso [W/m3] / (rho * cp) * 86400,  rho = play/(R_dry*Tlay)
    ! atm%plev/%tlay are surface-first while abso is top-down, hence the ka = nlay_m-kk+1 flip.
    !
    ! abso is left as the solver returns it: the topmost solid (terrain/building) cell of
    ! every column carries the whole surface-absorbed flux (~500 W/m2, abso ~ 20 W/m3) --
    ! correct energy accounting for the ground/roof, but not an air heating rate, so hr is
    ! set to HR_FILL (_FillValue) inside all solid cells (bottom `kt` layers of the column).
    allocate (hr, mold=abso)
    hr = HR_FILL
    nlay_m = size(abso, 1)
    associate (C1 => pprts_solver%C_one)
      do j = 1, nyp
        do i = 1, nxp
          icol = i + (j - 1) * nxp
          kt = kterr(C1%xs + i, C1%ys + j) ! solid layers from the surface up, this column
          do kk = 1, nlay_m - kt ! air cells only; solid cells at the bottom keep HR_FILL
            ka = nlay_m - kk + 1
            play_pa = 0.5_ireals * (atm%plev(ka, icol) + atm%plev(ka + 1, icol)) * 1.e2_ireals
            hr(kk, i, j) = abso(kk, i, j) / (play_pa / (R_DRY_AIR * atm%tlay(ka, icol)) * CP_DRY_AIR) * 86400._ireals
          end do
        end do
      end do
    end associate

  contains

    !> Turn every solid (terrain/building) cell of the cropped domain into a TenStream
    !> building box. Mirrors `ex_pprts_specint_buildings_from_file`, but the geometry
    !> comes from kterr(i,j) (number of solid layers per column) instead of a netcdf field.
    subroutine setup_buildings()
      integer(iintegers) :: i, j, kk, gk
      integer(iintegers) :: Nbuildings, Nfaces
      integer(mpiint) :: ierr
      logical :: lflg
      real(ireals) :: v

      associate (C1 => pprts_solver%C_one)

        Nbuildings = 0
        do j = 1, Ny_glob
          do i = 1, Nx_glob
            do kk = 1, kterr(i, j)
              if (have_box(C1%glob_zm - kk + 1, i, j)) Nbuildings = Nbuildings + 1
            end do
          end do
        end do
        Nfaces = Nbuildings * 6
        print *, 'rank '//toStr(myid)//' has '//toStr(Nbuildings)//' building cells ('//toStr(Nfaces)//' faces)'

        call init_buildings(buildings_solar, &
          & [integer(iintegers) :: 6, C1%zm, C1%xm, C1%ym], Nfaces, ierr); call CHKERR(ierr)

        Nbuildings = 0
        do j = 1, Ny_glob
          do i = 1, Nx_glob
            do kk = 1, kterr(i, j)
              gk = C1%glob_zm - kk + 1
              if (have_box(gk, i, j)) then
                Nbuildings = Nbuildings + 1
                call fill_cell_with_building(buildings_solar, Nbuildings, gk, i, j)
                if (associated(solar_albedo_2d)) then ! PALM per-column SW surface albedo on all 6 faces
                  buildings_solar%albedo((Nbuildings - 1) * 6 + 1:Nbuildings * 6) = &
                    & solar_albedo_2d(i - C1%xs, j - C1%ys)
                else
                  buildings_solar%albedo((Nbuildings - 1) * 6 + 1:Nbuildings * 6) = building_albedo_sol
                end if
              end if
            end do
          end do
        end do
        call get_petsc_opt('', "-override_buildings_albedo", v, lflg, ierr); call CHKERR(ierr)
        if (lflg) buildings_solar%albedo(:) = v

        call clone_buildings(buildings_solar, buildings_thermal, l_copy_data=.true., ierr=ierr); call CHKERR(ierr)
        buildings_thermal%albedo(:) = building_albedo_th
        if (.not. allocated(buildings_thermal%temp)) allocate (buildings_thermal%temp(Nfaces))
        buildings_thermal%temp(:) = building_temp
        call get_petsc_opt('', "-override_buildings_temperature", v, lflg, ierr); call CHKERR(ierr)
        if (lflg) buildings_thermal%temp(:) = v

        call check_buildings_consistency(buildings_solar, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)
        call check_buildings_consistency(buildings_thermal, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    !> Set every building face temperature to the air temperature of the nearest air
    !> voxel: the first merged-grid layer above the terrain top of that column
    !> (atm%tlay is surface-first). -override_buildings_temperature still forces a constant.
    subroutine set_building_temps_from_air()
      integer(iintegers) :: m, idx(4), il, jl, ka, natm
      integer(mpiint) :: ie
      real(ireals) :: v
      logical :: lflg

      call get_petsc_opt('', "-override_buildings_temperature", v, lflg, ie); call CHKERR(ie)
      if (lflg) then
        buildings_thermal%temp(:) = v
        return
      end if

      natm = size(atm%tlay, 1)
      associate (C1 => pprts_solver%C_one)
        do m = 1, size(buildings_thermal%iface)
          call ind_1d_to_nd(buildings_thermal%da_offsets, buildings_thermal%iface(m), idx) ! [face,k,i,j] local
          il = idx(3)
          jl = idx(4)
          ka = min(kterr(C1%xs + il, C1%ys + jl) + 1, natm) ! atm layer just above the terrain top
          buildings_thermal%temp(m) = atm%tlay(ka, il + (jl - 1) * nxp)
        end do
      end associate
    end subroutine

    subroutine fill_cell_with_building(B, m, gk, gi, gj)
      type(t_pprts_buildings), intent(inout) :: B
      integer(iintegers), intent(in) :: m, gk, gi, gj
      integer(iintegers) :: k, i, j, iface
      k = gk - pprts_solver%C_one%zs
      i = gi - pprts_solver%C_one%xs
      j = gj - pprts_solver%C_one%ys
      do iface = 1, 6
        B%iface((m - 1) * 6 + iface) = faceidx_by_cell_plus_offset(B%da_offsets, k, i, j, iface)
      end do
    end subroutine

    logical function have_box(gk, gi, gj)
      integer(iintegers), intent(in) :: gk, gi, gj
      have_box = all([ &
        & is_inrange(gk, pprts_solver%C_one%zs + 1, pprts_solver%C_one%ze + 1), &
        & is_inrange(gi, pprts_solver%C_one%xs + 1, pprts_solver%C_one%xe + 1), &
        & is_inrange(gj, pprts_solver%C_one%ys + 1, pprts_solver%C_one%ye + 1)])
    end function
  end subroutine

  subroutine example_palm_cld_file_pprts(specint, comm, &
                                         cldfile, atm_filename, outfile, &
                                         albedo_th, albedo_sol, lsolar, lthermal, &
                                         phi0, theta0, tstart, tend, tinc, &
                                         ix0_in, ix1_in, iy0_in, iy1_in, ztop, &
                                         psrfc_hPa, qv_scale, lwc_scale, cldN, &
                                         building_albedo_sol, building_albedo_th, building_temp, deflate, &
                                         lAg_from_palm)

    character(len=*), intent(in) :: specint
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    logical, intent(in) :: lAg_from_palm ! use PALM`s rad_sw_out/rad_sw_in as a per-column SW surface albedo
    real(ireals), intent(in) :: phi0, theta0
    integer(iintegers), intent(in) :: tstart, tend, tinc
    integer(iintegers), intent(in) :: ix0_in, ix1_in, iy0_in, iy1_in
    real(ireals), intent(in) :: ztop, psrfc_hPa, qv_scale, lwc_scale, cldN
    real(ireals), intent(in) :: building_albedo_sol, building_albedo_th, building_temp
    integer(iintegers), intent(in) :: deflate ! netcdf deflate level for the output (0 = none, fastest to (re)read)

    character(len=default_str_len) :: groups(3), dimnames(3)
    character(len=10*default_str_len) :: outfile_bldg ! big building-face arrays go into a sibling file
    integer(iintegers) :: nout

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso, hr
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso, ghr
    real(ireals), pointer :: p_alb2d(:, :) ! per-column SW surface albedo from PALM (null unless -Ag_from_palm)

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm
    type(t_pprts_buildings), allocatable :: buildings_solar, buildings_thermal

    real(ireals) :: dx, dy, origin_z
    real(ireals), allocatable :: time(:), zw_full(:), zu_full(:), zlev(:), zc(:), zusi(:, :)
    integer(iintegers), allocatable :: kterr(:, :)

    integer(mpiint) :: myid, ierr
    integer(iintegers) :: Nx_full, Ny_full, Nx, Ny, nlay, nlev
    integer(iintegers) :: ix0, ix1, iy0, iy1
    integer(iintegers) :: is, js, nxp, nyp, it, k
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals), allocatable, dimension(:, :, :) :: plev, tlev, tlay, h2ovmr, lwc, reliq

    call mpi_comm_rank(comm, myid, ierr)

    call load_meta_data(comm, cldfile, time, zw_full, zu_full, Nx_full, Ny_full, dx, dy, origin_z, ierr); call CHKERR(ierr)

    ! resolve crop window (1-based, inclusive, global indices)
    ix0 = merge(ix0_in, 1_iintegers, ix0_in .gt. 0)
    iy0 = merge(iy0_in, 1_iintegers, iy0_in .gt. 0)
    ix1 = merge(ix1_in, Nx_full, ix1_in .gt. 0)
    iy1 = merge(iy1_in, Ny_full, iy1_in .gt. 0)
    ix1 = min(ix1, Nx_full); iy1 = min(iy1, Ny_full)
    if (ix1 .lt. ix0 .or. iy1 .lt. iy0) call CHKERR(1_mpiint, 'invalid crop window')
    Nx = ix1 - ix0 + 1
    Ny = iy1 - iy0 + 1

    ! vertical extent: use all PALM levels up to -ztop
    nlay = 0
    do k = 1, size(zw_full) - 1
      if (zw_full(k + 1) .le. ztop) nlay = k
    end do
    if (nlay .lt. 2) call CHKERR(1_mpiint, '-ztop too low, need at least 2 layers')
    nlev = nlay + 1
    allocate (zlev(nlev), source=zw_full(1:nlev))
    allocate (zc(nlay))
    zc = 0.5_ireals * (zlev(1:nlay) + zlev(2:nlev))

    if (myid .eq. 0) then
      print *, 'crop window x:', ix0, ix1, ' (', Nx, ') y:', iy0, iy1, ' (', Ny, ')'
      print *, 'using', nlay, 'dynamics layers up to z =', zlev(nlev), ' [m]'
    end if

    if (tstart .lt. 1 .or. tend .gt. size(time) .or. tstart .gt. tend) &
      & call CHKERR(1_mpiint, 'invalid -tstart/-tend, file has '//toStr(size(time))//' timesteps')

    call load_topography(comm, cldfile, ix0, iy0, Nx, Ny, zc, zusi, kterr, ierr); call CHKERR(ierr)

    call domain_decompose_2d_petsc(comm, Nx_global=Nx, Ny_global=Ny, &
      & Nx_local=nxp, Ny_local=nyp, xStart=is, yStart=js, &
      & nxproc=nxproc, nyproc=nyproc, ierr=ierr); call CHKERR(ierr)
    is = is + 1; js = js + 1 ! fortran 1-based, crop-relative

    call allocate_pprts_solver_from_commandline(pprts_solver, default_solver='3_10', ierr=ierr); call CHKERR(ierr)

    ! per-rank subdomain arrays only -- nothing global is kept, so this scales to large PALM domains
    allocate (plev(nlev, nxp, nyp))
    allocate (tlev(nlev, nxp, nyp))
    allocate (tlay(nlay, nxp, nyp))
    allocate (h2ovmr(nlay, nxp, nyp))
    allocate (lwc(nlay, nxp, nyp))
    allocate (reliq(nlay, nxp, nyp))

    p_alb2d => null()
    if (lAg_from_palm) allocate (p_alb2d(nxp, nyp))

    ! Put the (potentially huge, ~1e8 element) building-face arrays in a sibling file so
    ! that opening the main file in ncview/paraview does not have to scan them.
    nout = len_trim(outfile)
    if (nout .ge. 4 .and. outfile(max(1, nout - 2):nout) .eq. '.nc') then
      outfile_bldg = outfile(1:nout - 3)//'_buildings.nc'
    else
      outfile_bldg = trim(outfile)//'_buildings.nc'
    end if

    do it = tstart, tend, tinc
      if (myid .eq. 0) print *, ''
      if (myid .eq. 0) print *, '==== timestep', it, ' (t =', time(it), 's) ===='

      ! each rank reads its own hyperslab: global netcdf start = crop origin + local start
      call load_timestep_data(comm, cldfile, it, ix0 + is - 1, iy0 + js - 1, is, js, nxp, nyp, &
        & zlev, zc, kterr, psrfc_hPa, qv_scale, lwc_scale, cldN, &
        & plev, tlev, tlay, h2ovmr, lwc, reliq, ierr); call CHKERR(ierr)

      if (lAg_from_palm) then
        call load_surface_albedo(cldfile, it, ix0 + is - 1, iy0 + js - 1, nxp, nyp, &
          & size(zw_full, kind=iintegers), albedo_sol, p_alb2d, ierr); call CHKERR(ierr)
        if (myid .eq. 0) print *, 'PALM surface SW albedo (rank 0 subdomain): min/max', minval(p_alb2d), maxval(p_alb2d)
      end if

      call run_lw_sw(specint, comm, pprts_solver, atm, atm_filename, &
        & nxproc, nyproc, dx, dy, phi0, theta0, &
        & albedo_th, albedo_sol, lsolar, lthermal, time(it), &
        & Nx, Ny, kterr, &
        & building_albedo_sol, building_albedo_th, building_temp, &
        & p_alb2d, &
        & plev, tlev, tlay, &
        & h2ovmr, lwc, reliq, &
        & buildings_solar, buildings_thermal, &
        & edir, edn, eup, abso, hr)

      if (len_trim(outfile) .gt. 0) call dump_timestep(it)
    end do

    if (associated(p_alb2d)) deallocate (p_alb2d)
    call specint_pprts_destroy(specint, pprts_solver, lfinalizepetsc=.true., ierr=ierr)
    call destroy_tenstr_atm(atm)

  contains

    subroutine dump_timestep(it)
      integer(iintegers), intent(in) :: it
      groups(1) = trim(outfile)
      groups(3) = trim(toStr(it))
      associate (C1 => pprts_solver%C_one1, C => pprts_solver%C_one, Ca1 => pprts_solver%C_one_atm1_box)
        dimnames(1) = 'zlev'; dimnames(2) = 'nx'; dimnames(3) = 'ny'
        if (allocated(edir)) then
          call gather_all_toZero(C1, edir, gedir)
          if (myid .eq. 0) then
            groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames, deflate_lvl=int(deflate)); call CHKERR(ierr)
          end if
        end if
        call gather_all_toZero(C1, edn, gedn)
        call gather_all_toZero(C1, eup, geup)
        call gather_all_toZero(C, abso, gabso)
        call gather_all_toZero(C, hr, ghr)
        if (myid .eq. 0) then
          groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames, deflate_lvl=int(deflate)); call CHKERR(ierr)
          groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames, deflate_lvl=int(deflate)); call CHKERR(ierr)
          dimnames(1) = 'zlay'
          groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames, deflate_lvl=int(deflate)); call CHKERR(ierr)
          groups(2) = 'hr'
          call ncwrite(groups, ghr, ierr, dimnames=dimnames, deflate_lvl=int(deflate), fill_value=HR_FILL)
          call CHKERR(ierr)
          dimnames(1) = 'nlev'
          groups(2) = 'zlev'; call ncwrite(groups, &
 & pprts_solver%atm%hhl(0, Ca1%zs:Ca1%ze, Ca1%xs, Ca1%ys), ierr, dimnames=dimnames(1:1), deflate_lvl=int(deflate))
          call CHKERR(ierr)
          dimnames(1) = 'nlay'
          groups(2) = 'zlay'; call ncwrite(groups, &
 & (pprts_solver%atm%hhl(0, Ca1%zs:Ca1%ze - 1, Ca1%xs, Ca1%ys) + &
 &  pprts_solver%atm%hhl(0, Ca1%zs + 1:Ca1%ze, Ca1%xs, Ca1%ys)) * 0.5_ireals, &
 & ierr, dimnames=dimnames(1:1), deflate_lvl=int(deflate)); call CHKERR(ierr)

          if (allocated(edir)) call put_attrs(outfile, 'edir', it, 'W m-2', 'downward direct solar irradiance on horizontal plane')
          call put_attrs(outfile, 'edn', it, 'W m-2', 'downward diffuse irradiance (solar diffuse + thermal)')
          call put_attrs(outfile, 'eup', it, 'W m-2', 'upward diffuse irradiance (solar diffuse + thermal)')
          call put_attrs(outfile, 'abso', it, 'W m-3', &
            & 'radiative flux-divergence heating rate per volume (solid cells hold the surface absorption)')
          call put_attrs(outfile, 'hr', it, 'K day-1', &
            & 'radiative heating rate of the air (_FillValue inside solid terrain/building cells)')
          call put_attrs(outfile, 'zlev', it, 'm', 'height of merged-grid level interfaces above the RT-grid base')
          call put_attrs(outfile, 'zlay', it, 'm', 'height of merged-grid layer centers above the RT-grid base')
        end if
      end associate

      ! Building face fluxes: every rank owns a disjoint set of faces. Instead of one
      ! variable per rank, write a single concatenated variable per field plus a
      ! `buildings_idx(4, Nglobal)` companion holding the global [face, k, i, j] index
      ! of each entry (k counted from the surface upward). This mirrors
      ! `dump_input_buildings` in specint/specint_pprts.F90 and lets a reader restore
      ! the global ordering / map the flat list back onto the 3d grid.
      call mpi_barrier(comm, ierr) ! let rank-0 finish its serial writes before the collective ones
      if (allocated(buildings_solar)) call dump_buildings_faces(buildings_solar, 'buildings', lsolar, it)
      if (allocated(buildings_thermal)) call dump_buildings_faces(buildings_thermal, 'buildings_thermal', .false., it)
    end subroutine

    subroutine dump_buildings_faces(B, tag, lhave_edir, it)
      type(t_pprts_buildings), intent(in) :: B
      character(len=*), intent(in) :: tag
      logical, intent(in) :: lhave_edir
      integer(iintegers), intent(in) :: it

      integer(iintegers) :: Nloc, Nglob, bStart, m, idx(4)
      integer(iintegers), allocatable :: bidx(:, :)
      character(len=default_str_len) :: g(3), dn(2), facedim

      associate (C1 => pprts_solver%C_one)
        Nloc = size(B%iface)
        call imp_scan_sum(comm, Nloc, bStart, ierr); call CHKERR(ierr)
        bStart = 1 + bStart - Nloc ! this rank's 1-based start in the global concatenation
        call imp_allreduce_sum(comm, Nloc, Nglob)

        allocate (bidx(4, Nloc))
        do m = 1, Nloc
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx) ! [face, k, i, j] in local buildings-DA indices
          idx(3:4) = idx(3:4) + [C1%xs, C1%ys] ! -> global (crop-relative) i, j
          idx(2) = C1%zm - idx(2) + 1          ! -> k counted from the surface upward
          bidx(:, m) = idx
        end do

        g(1) = trim(outfile_bldg)
        g(3) = trim(toStr(it))
        facedim = trim(tag)//'_nfaces'

        g(2) = trim(tag)//'_idx'
        dn = [character(len=default_str_len) :: 'fkij', facedim]
        call ncwrite(comm=comm, groups=g, arr=bidx, ierr=ierr, &
          & arr_shape=[integer :: 4, int(Nglob)], dimnames=dn, &
          & startp=[integer :: 1, int(bStart)], countp=shape(bidx), deflate_lvl=int(deflate), verbose=.false.)
        call CHKERR(ierr)

        if (lhave_edir .and. allocated(B%edir)) &
          & call write_face_var(trim(tag)//'_edir', facedim, B%edir, Nglob, bStart, it)
        if (allocated(B%incoming)) call write_face_var(trim(tag)//'_incoming', facedim, B%incoming, Nglob, bStart, it)
        if (allocated(B%outgoing)) call write_face_var(trim(tag)//'_outgoing', facedim, B%outgoing, Nglob, bStart, it)
      end associate

      call mpi_barrier(comm, ierr) ! all collective writes done -> rank 0 may open serially for attributes
      if (myid .eq. 0) then
        call put_attrs(outfile_bldg, trim(tag)//'_idx', it, '1', &
          & 'per-face global index: [face(1=top,2=bottom,3=left,4=right,5=rear,6=front), k_from_surface, i, j]')
        if (lhave_edir .and. allocated(B%edir)) &
          & call put_attrs(outfile_bldg, trim(tag)//'_edir', it, 'W m-2', 'direct solar irradiance on the building face')
        if (allocated(B%incoming)) &
          & call put_attrs(outfile_bldg, trim(tag)//'_incoming', it, 'W m-2', 'radiation incident on the building face')
        if (allocated(B%outgoing)) &
          & call put_attrs(outfile_bldg, trim(tag)//'_outgoing', it, 'W m-2', 'radiation leaving the building face')
      end if
    end subroutine

    subroutine put_attrs(fname, varbase, it, units, long_name)
      character(len=*), intent(in) :: fname, varbase, units, long_name
      integer(iintegers), intent(in) :: it
      character(len=default_str_len) :: vname
      integer(mpiint) :: ie
      vname = trim(varbase)//'.'//trim(toStr(it))
      call set_attribute(fname, vname, 'units', units, ie); call CHKERR(ie)
      call set_attribute(fname, vname, 'long_name', long_name, ie); call CHKERR(ie)
    end subroutine

    subroutine write_face_var(varname, facedim, arr, Nglob, bStart, it)
      character(len=*), intent(in) :: varname, facedim
      real(ireals), intent(in) :: arr(:)
      integer(iintegers), intent(in) :: Nglob, bStart, it
      character(len=default_str_len) :: g(3), dn(1)
      g(1) = trim(outfile_bldg)
      g(2) = trim(varname)
      g(3) = trim(toStr(it))
      dn(1) = trim(facedim)
      call ncwrite(comm=comm, groups=g, arr=arr, ierr=ierr, &
        & arr_shape=[integer :: int(Nglob)], dimnames=dn, &
        & startp=[integer :: int(bStart)], countp=shape(arr), deflate_lvl=int(deflate), verbose=.false.)
      call CHKERR(ierr)
    end subroutine
  end subroutine
end module

program main
  use mpi
  use m_data_parameters, only: init_mpi_data_parameters, mpiint, iintegers, ireals, &
                               default_str_len, share_dir
  use m_helper_functions, only: CHKERR, get_petsc_opt
  use m_tenstream_options, only: read_commandline_options
  use m_example_palm_cld_file, only: example_palm_cld_file_pprts

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg, lAg_from_palm
  character(len=10*default_str_len) :: cldfile, outfile
  character(len=default_str_len) :: atm_filename, specint
  real(ireals) :: Ag, Ag_solar, Ag_thermal, phi0, theta0
  real(ireals) :: ztop, psrfc_hPa, qv_scale, lwc_scale, cldN
  real(ireals) :: building_albedo_sol, building_albedo_th, building_temp
  integer(iintegers) :: tstart, tend, tinc
  integer(iintegers) :: ix0, ix1, iy0, iy1
  integer(iintegers) :: deflate

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  specint = 'no default set'
  call get_petsc_opt('', '-specint', specint, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need -specint <repwvl|ecckd|rrtmg|...>')

  cldfile = 'unset'
  call get_petsc_opt('', '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need -cld <PALM_..._3d.000.nc>')

  outfile = ''
  call get_petsc_opt('', '-out', outfile, lflg, ierr); call CHKERR(ierr)

  atm_filename = share_dir//'tenstream_default.atm'
  call get_petsc_opt('', '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  Ag = 0.1_ireals
  call get_petsc_opt('', "-Ag", Ag, lflg, ierr); call CHKERR(ierr)
  Ag_solar = Ag; Ag_thermal = Ag
  call get_petsc_opt('', "-Ag_solar", Ag_solar, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-Ag_thermal", Ag_thermal, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  lthermal = .true.
  call get_petsc_opt('', "-solar", lsolar, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)

  phi0 = 180
  call get_petsc_opt('', "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 40
  call get_petsc_opt('', "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  tstart = 1; tend = 1; tinc = 1
  call get_petsc_opt('', "-tstart", tstart, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-tend", tend, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-tinc", tinc, lflg, ierr); call CHKERR(ierr)

  ix0 = -1; ix1 = -1; iy0 = -1; iy1 = -1
  call get_petsc_opt('', "-ix0", ix0, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-ix1", ix1, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-iy0", iy0, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-iy1", iy1, lflg, ierr); call CHKERR(ierr)

  ztop = 1.e9_ireals
  call get_petsc_opt('', "-ztop", ztop, lflg, ierr); call CHKERR(ierr)

  psrfc_hPa = 1000._ireals
  call get_petsc_opt('', "-psrfc", psrfc_hPa, lflg, ierr); call CHKERR(ierr)

  ! LES unit fudge factors (attributes often lie): raw_field * scale -> expected unit
  qv_scale = 1._ireals   ! qv given as kg/kg -> set 1e-3 if the file is actually g/kg
  lwc_scale = 1.e3_ireals ! ql given as kg/kg -> g/kg ; set 1.0 if the file is already g/kg
  call get_petsc_opt('', "-qv_scale", qv_scale, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-lwc_scale", lwc_scale, lflg, ierr); call CHKERR(ierr)

  cldN = 200._ireals ! cloud droplet number concentration [1/cm3] for the effective radius
  call get_petsc_opt('', "-cldN", cldN, lflg, ierr); call CHKERR(ierr)

  building_albedo_sol = 0.2_ireals
  building_albedo_th = 0.05_ireals
  building_temp = 288._ireals
  call get_petsc_opt('', "-building_albedo", building_albedo_sol, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-building_albedo_solar", building_albedo_sol, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-building_albedo_thermal", building_albedo_th, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt('', "-building_temp", building_temp, lflg, ierr); call CHKERR(ierr)

  deflate = 0 ! netcdf deflate level; 0 keeps the file fast to (re)open in ncview/paraview
  call get_petsc_opt('', "-deflate", deflate, lflg, ierr); call CHKERR(ierr)

  lAg_from_palm = .true. ! derive the SW surface albedo per column from PALM rad_sw_out/rad_sw_in
  call get_petsc_opt('', "-Ag_from_palm", lAg_from_palm, lflg, ierr); call CHKERR(ierr) ! -Ag_from_palm false -> use scalar -Ag

  call example_palm_cld_file_pprts(specint, mpi_comm_world, &
    & cldfile, atm_filename, outfile, &
    & Ag_thermal, Ag_solar, lsolar, lthermal, &
    & phi0, theta0, tstart, tend, tinc, &
    & ix0, ix1, iy0, iy1, ztop, &
    & psrfc_hPa, qv_scale, lwc_scale, cldN, &
    & building_albedo_sol, building_albedo_th, building_temp, deflate, &
    & lAg_from_palm)

  call mpi_finalize(ierr)
end program
