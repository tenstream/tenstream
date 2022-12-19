module m_example_pprts_specint_lw_sw

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, default_str_len

  use m_helper_functions, only: linspace, CHKERR, spherical_2_cartesian, meanval, get_petsc_opt

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, abso2hr

  use m_petsc_helpers, only: getvecpointer, restorevecpointer
  use m_netcdfio, only: ncwrite

  implicit none

contains
  subroutine ex_pprts_specint_lw_sw(specint, comm, nxp, nyp, nzp, dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & gedir, gedn, geup, gabso, &
      & vlwc, viwc, vreff, vreice, outfile)
    character(len=*), intent(in) :: specint           ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(in) :: nxp, nyp, nzp   ! local domain size for each rank
    real(ireals), intent(in) :: dx, dy                ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0          ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lthermal, lsolar           ! switches to enable/disable spectral integration
    character(len=*), intent(in) :: atm_filename ! ='afglus_100m.dat'
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: gedir, gedn, geup, gabso
    real(ireals), intent(in), optional :: vlwc, viwc            ! liquid/ice water content to be set in a layer
    real(ireals), intent(in), optional :: vreff, vreice         ! liquid/ice effective radius
    character(len=*), intent(in), optional :: outfile

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso, hr ! [nlev_merged(-1), nxp, nyp]

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, myid, N_ranks_x, N_ranks_y, ierr

    real(ireals), dimension(nzp + 1, nxp, nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp + 1, nxp, nyp), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    real(ireals), dimension(nzp, nxp, nyp), target :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp, nxp, nyp), target :: lwc, reliq, iwc, reice

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    logical :: lflg

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev, icld, iter, icollapse
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)
    character(len=default_str_len) :: groups(2), dimnames(3)
    real(ireals), pointer :: z(:, :, :, :) => null(), z1d(:) => null() ! dim Nz+1

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, plwc, preliq, piwc, preice
    real(ireals), pointer, dimension(:, :) :: ph2ovmr, po3vmr, pco2vmr, pch4vmr, pn2ovmr, po2vmr

    real(ireals) :: vmr, sundir(3)

    logical, parameter :: ldebug = .true.

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm

    call MPI_COMM_SIZE(comm, numnodes, ierr)
    call MPI_COMM_RANK(comm, myid, ierr)

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if (N_ranks_y * N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    end if
    if (ldebug .and. myid .eq. 0) print *, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate (nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate (nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 50K
    do k = 1, nzp + 1
      plev(k, :, :) = linspace(k, [1e3_ireals, 100._ireals], nzp + 1)
      tlev(k, :, :) = linspace(k, [288._ireals, 250._ireals], nzp + 1)
    end do

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    h2ovmr = .007
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-h2o", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) h2ovmr = vmr

    o3vmr = 3e-8
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-o3", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) o3vmr = vmr

    co2vmr = 400e-6
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-co2", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) co2vmr = vmr

    ch4vmr = 1.7e-6
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-ch4", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) ch4vmr = vmr

    n2ovmr = 3.2e-7
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-n2o", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) n2ovmr = vmr

    o2vmr = .2
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-o2", vmr, lflg, ierr); call CHKERR(ierr)
    if (lflg) o2vmr = vmr

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    if (present(vlwc)) then
      icld = int(real(nzp + 1) / 2)
      lwc(icld, :, :) = vlwc
      if (present(vreff)) then
        reliq(icld, :, :) = vreff
      else
        reliq(icld, :, :) = 10
      end if
      tlev(icld + 1, :, :) = tlev(icld, :, :)
    end if

    iwc = 0
    reice = 0

    if (present(viwc)) then
      icld = nzp
      iwc(icld, :, :) = viwc
      if (present(vreice)) then
        reice(icld, :, :) = vreice
      else
        reice(icld, :, :) = 60
      end if
      tlev(icld + 1, :, :) = tlev(icld, :, :)
    end if

    if (myid .eq. 0 .and. ldebug) print *, 'Setup Atmosphere...'

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    plwc(1:size(lwc, 1), 1:size(lwc, 2) * size(lwc, 3)) => lwc
    preliq(1:size(reliq, 1), 1:size(reliq, 2) * size(reliq, 3)) => reliq
    piwc(1:size(iwc, 1), 1:size(iwc, 2) * size(iwc, 3)) => iwc
    preice(1:size(reice, 1), 1:size(reice, 2) * size(reice, 3)) => reice
    pco2vmr(1:size(co2vmr, 1), 1:size(co2vmr, 2) * size(co2vmr, 3)) => co2vmr
    ph2ovmr(1:size(h2ovmr, 1), 1:size(h2ovmr, 2) * size(h2ovmr, 3)) => h2ovmr
    po3vmr(1:size(o3vmr, 1), 1:size(o3vmr, 2) * size(o3vmr, 3)) => o3vmr
    pch4vmr(1:size(ch4vmr, 1), 1:size(ch4vmr, 2) * size(ch4vmr, 3)) => ch4vmr
    pn2ovmr(1:size(n2ovmr, 1), 1:size(n2ovmr, 2) * size(n2ovmr, 3)) => n2ovmr
    po2vmr(1:size(o2vmr, 1), 1:size(o2vmr, 2) * size(o2vmr, 3)) => o2vmr

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_lwc=plwc, d_reliq=preliq, &
                          d_iwc=piwc, d_reice=preice, &
                          d_co2vmr=pco2vmr, &
                          d_h2ovmr=ph2ovmr, &
                          d_o3vmr=po3vmr, &
                          d_ch4vmr=pch4vmr, &
                          d_n2ovmr=pn2ovmr, &
                          d_o2vmr=po2vmr)

    sundir = spherical_2_cartesian(phi0, theta0)

    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10', ierr); call CHKERR(ierr)

    icollapse = 1
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-icollapse", icollapse, lflg, ierr)

    iter = 1
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-iter", iter, lflg, ierr)

    do k = 1, iter
      call specint_pprts(specint, comm, pprts_solver, atm, nxp, nyp, &
                         dx, dy, sundir, &
                         albedo_th, albedo_sol + 0.1_ireals * iter, &
                         lthermal, lsolar, &
                         edir, edn, eup, abso, &
                         nxproc=nxproc, nyproc=nyproc, &
                         icollapse=icollapse, &
                         opt_time=real(k, ireals))
    end do

    allocate (hr(size(abso, 1), size(abso, 2), size(abso, 3)))
    call abso2hr(atm, abso, hr, ierr); call CHKERR(ierr)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          if (allocated(edir)) then
            print *, k, 'edir', meanval(edir(k, :, :)), 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), &
              & 'abso', meanval(abso(min(nlev - 1, k), :, :)), 'hr', meanval(hr(min(nlev - 1, k), :, :)) * 3600 * 24
          else
            print *, k, 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), &
              & 'abso', meanval(abso(min(nlev - 1, k), :, :)), 'hr', meanval(hr(min(nlev - 1, k), :, :)) * 3600 * 24
          end if
        end do
      end if

      if (allocated(edir)) &
        print *, 'surface :: direct flux', meanval(edir(nlev, :, :))
      print *, 'surface :: downw flux ', meanval(edn(nlev, :, :))
      print *, 'surface :: upward fl  ', meanval(eup(nlev, :, :))
      print *, 'surface :: absorption ', meanval(abso(nlev - 1, :, :))

      if (allocated(edir)) &
        print *, 'TOA :: direct flux', meanval(edir(1, :, :))
      print *, 'TOA :: downw flux ', meanval(edn(1, :, :))
      print *, 'TOA :: upward fl  ', meanval(eup(1, :, :))
      print *, 'TOA :: absorption ', meanval(abso(1, :, :))
    end if

    if (allocated(edir)) &
      & call gather_all_toZero(pprts_solver%C_one1, edir, gedir)
    call gather_all_toZero(pprts_solver%C_one1, edn, gedn)
    call gather_all_toZero(pprts_solver%C_one1, eup, geup)
    call gather_all_toZero(pprts_solver%C_one, abso, gabso)

    if (myid .eq. 0_mpiint .and. present(outfile)) then
      dimnames(1) = 'zlev'
      dimnames(2) = 'nx'
      dimnames(3) = 'ny'
      groups(1) = trim(outfile)
      if (lsolar) then
        groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames); call CHKERR(ierr)
      end if
      groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames); call CHKERR(ierr)
      groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames); call CHKERR(ierr)
      dimnames(1) = 'zlay'
      groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames); call CHKERR(ierr)

      print *, 'dumping z coords'
      associate (Ca1 => pprts_solver%C_one_atm1_box)
        call getVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
        dimnames(1) = 'nlev'
        groups(2) = 'zlev'
        call ncwrite(groups, z(0, Ca1%zs:Ca1%ze, Ca1%xs, Ca1%ys), ierr, dimnames=dimnames(1:1))
        call CHKERR(ierr)
        dimnames(1) = 'nlay'
        groups(2) = 'zlay'
        call ncwrite(groups, &
                     & (z(0, Ca1%zs:Ca1%ze - 1, Ca1%xs, Ca1%ys) &
                     & + z(0, Ca1%zs + 1:Ca1%ze, Ca1%xs, Ca1%ys) &
                     & )*.5_ireals, &
                     & ierr, dimnames=dimnames(1:1))
        call CHKERR(ierr)
        call restoreVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
      end associate
    end if

    ! Tidy up
    call specint_pprts_destroy(specint, pprts_solver, lfinalizepetsc=.true., ierr=ierr)
    call destroy_tenstr_atm(atm)
  end subroutine

end module
