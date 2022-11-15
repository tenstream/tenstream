module test_pprts_specint
  use iso_fortran_env, only: real32, real64

  use m_data_parameters, only: &
    & default_str_len, &
    & iintegers, &
    & init_mpi_data_parameters, &
    & ireals, &
    & mpiint, &
    & one, &
    & share_dir, &
    & zero

  use m_helper_functions, only: CHKERR, spherical_2_cartesian

  use m_tenstream_options, only: read_commandline_options

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  use m_pprts_base, only: t_solver_3_10

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver
  type(t_tenstr_atm) :: atm

  integer(iintegers), parameter :: nxp = 3, nyp = 3, nzp = 10 ! local domain size for each rank

  ! MPI variables and domain decomposition sizes
  integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y, ierr

  real(ireals), parameter :: dx = 100, dy = dx              ! horizontal grid spacing in [m]
  real(ireals), parameter :: phi0 = 180, theta0 = 60        ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
  real(ireals), parameter :: albedo_th = .1, albedo_sol = .3 ! broadband ground albedo for solar and thermal spectrum
  real(ireals), parameter :: atolerance = 1                  ! absolute tolerance when regression testing fluxes

  integer(iintegers), allocatable :: nxproc(:), nyproc(:)

  real(ireals), dimension(nzp + 1, nxp, nyp), target :: plev ! pressure on layer interfaces [hPa]
  real(ireals), dimension(nzp + 1, nxp, nyp), target :: tlev ! Temperature on layer interfaces [K]

  ! Layer values for the atmospheric constituents -- those are actually all
  ! optional and if not provided, will be taken from the background profile file (atm_filename)
  ! see interface of `pprts_rrtmg()` for units
  ! real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

  ! Liquid water cloud content [g/kg] and effective radius in micron
  real(ireals), dimension(nzp, nxp, nyp), target :: lwc, reliq

  integer(iintegers) :: k, nlev, icld

  real(ireals), pointer, dimension(:, :) :: pplev, ptlev, plwc, preliq

  logical, parameter :: ldebug = .true.
  logical :: lthermal, lsolar

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if (N_ranks_y * N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    end if
    if (myid .eq. 0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate (nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate (nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 50K
    do k = 1, nzp + 1
      plev(k, :, :) = 1000_ireals - real(k - 1, ireals) * 500._ireals / real(nzp, ireals)
      tlev(k, :, :) = 288._ireals - real(k - 1, ireals) * 50._ireals / real(nzp, ireals)
    end do

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    icld = int(real(nzp + 1) / 2)
    lwc(icld, :, :) = 1e-2_ireals
    reliq(icld, :, :) = 10._ireals

    tlev(icld, :, :) = 288._ireals
    tlev(icld + 1, :, :) = tlev(icld, :, :)

    ! Setup Atmosphere
    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(plev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    plwc(1:size(lwc, 1), 1:size(lwc, 2) * size(lwc, 3)) => lwc
    preliq(1:size(reliq, 1), 1:size(reliq, 2) * size(reliq, 3)) => reliq

    call setup_tenstr_atm(comm, .false., './afglus_100m.dat', &
                          pplev, ptlev, atm, &
                          d_lwc=plwc, d_reliq=preliq)
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    deallocate (nxproc)
    deallocate (nyproc)
    call destroy_tenstr_atm(atm)
  end subroutine teardown

  @test(npes=[1, 2])
  subroutine specint_rrtm_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this
    character(len=*), parameter :: specint = 'rrtmg'

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation:'
    lthermal = .false.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, &
                       opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual(298.65, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(165.81, edn(nlev, 1, 1), atolerance, 'solar at surface :: downw flux not correct')
    @mpiassertEqual(139.34, eup(nlev, 1, 1), atolerance, 'solar at surface :: upward fl  not correct')
    @mpiassertEqual(1.425e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'solar at surface :: absorption not correct')

    @mpiassertEqual(684.1109, edir(1, 1, 1), atolerance, 'solar at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar at TOA :: downw flux not correct')
    @mpiassertEqual(213.85, eup(1, 1, 1), atolerance, 'solar at TOA :: upward fl  not correct')
    @mpiassertEqual(7.552e-08, abso(1, 1, 1), atolerance * 1e-2, 'solar at TOA :: absorption not correct')

    @mpiassertEqual(511.60, edir(nlev - icld, 1, 1), atolerance, 'solar at icloud :: direct flux not correct')
    @mpiassertEqual(323.27, edir(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: direct flux not correct')
    @mpiassertEqual(167.71, edn(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: downw flux not correct')
    @mpiassertEqual(194.19, eup(nlev - icld, 1, 1), atolerance, 'solar at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0224, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation:'
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done'

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual(334.38, edn(nlev, 1, 1), atolerance, 'thermal at surface :: downw flux not correct')
    @mpiassertEqual(384.50, eup(nlev, 1, 1), atolerance, 'thermal at surface :: upward fl  not correct')
    @mpiassertEqual(-1.487e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'thermal at surface :: absorption not correct')

    @mpiassertEqual(0.0, edn(1, 1, 1), atolerance, 'thermal at TOA :: downw flux not correct')
    @mpiassertEqual(250.67, eup(1, 1, 1), atolerance, 'thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(-8.303e-07, abso(1, 1, 1), atolerance * 1e-2, 'thermal at TOA :: absorption not correct')

    @mpiassertEqual(323.49, edn(nlev - icld + 1, 1, 1), atolerance, 'thermal at icloud :: downw flux not correct')
    @mpiassertEqual(384.82, eup(nlev - icld, 1, 1), atolerance, 'thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(-0.199, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'thermal at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation:'
    lthermal = .true.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual(298.65, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(165.81 + 334.38, edn(nlev, 1, 1), atolerance, 'solar+thermal at surface :: downw flux not correct')
    @mpiassertEqual(139.34 + 384.50, eup(nlev, 1, 1), atolerance, 'solar+thermal at surface :: upward fl  not correct')
    @mpiassertEqual(1.425E-02-1.487E-02, abso(nlev-1,1,1), atolerance*1e-2, 'solar+thermal at surface :: absorption not correct')

    @mpiassertEqual(684.1109, edir(1, 1, 1), atolerance, 'solar+thermal at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar+thermal at TOA :: downw flux not correct')
    @mpiassertEqual(213.85 + 250.67, eup(1, 1, 1), atolerance, 'solar+thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(7.552e-08 - 8.303e-07, abso(1, 1, 1), atolerance * 1e-2, 'solar+thermal at TOA :: absorption not correct')

    @mpiassertEqual(511.60, edir(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: direct flux not correct')
    @mpiassertEqual(323.27, edir(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud+1 :: direct flux not correct')
    @mpiassertEqual(168.64 + 323.49, edn(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud :: downw flux not correct')
    @mpiassertEqual(194.18 + 384.82, eup(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0243 - 0.199, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar+thermal at icloud :: absorption not correct')

    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[1, 2])
  subroutine specint_repwvl_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this
    character(len=*), parameter :: specint = 'repwvl'

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation:'
    lthermal = .false.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, &
                       opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual( 63.05, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(373.33, edn(nlev, 1, 1), atolerance, 'solar at surface :: downw flux not correct')
    @mpiassertEqual(130.91, eup(nlev, 1, 1), atolerance, 'solar at surface :: upward fl  not correct')
    @mpiassertEqual(1.332e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'solar at surface :: absorption not correct')

    @mpiassertEqual(684.27, edir(1, 1, 1), atolerance, 'solar at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar at TOA :: downw flux not correct')
    @mpiassertEqual(241.15, eup(1, 1, 1), atolerance, 'solar at TOA :: upward fl  not correct')
    @mpiassertEqual(7.552e-08, abso(1, 1, 1), atolerance * 1e-2, 'solar at TOA :: absorption not correct')

    @mpiassertEqual(516.43, edir(nlev - icld, 1, 1), atolerance, 'solar at icloud :: direct flux not correct')
    @mpiassertEqual( 67.97, edir(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: direct flux not correct')
    @mpiassertEqual(391.26, edn(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: downw flux not correct')
    @mpiassertEqual(224.00, eup(nlev - icld, 1, 1), atolerance, 'solar at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0226, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation:'
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done'

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual(338.72, edn(nlev, 1, 1), atolerance, 'thermal at surface :: downw flux not correct')
    @mpiassertEqual(384.44, eup(nlev, 1, 1), atolerance, 'thermal at surface :: upward fl  not correct')
    @mpiassertEqual(-1.333e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'thermal at surface :: absorption not correct')

    @mpiassertEqual(0.0, edn(1, 1, 1), atolerance, 'thermal at TOA :: downw flux not correct')
    @mpiassertEqual(248.75, eup(1, 1, 1), atolerance, 'thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(-4.924e-16, abso(1, 1, 1), atolerance * 1e-2, 'thermal at TOA :: absorption not correct')

    @mpiassertEqual(329.87, edn(nlev - icld + 1, 1, 1), atolerance, 'thermal at icloud :: downw flux not correct')
    @mpiassertEqual(381.32, eup(nlev - icld, 1, 1), atolerance, 'thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(-0.203, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'thermal at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation:'
    lthermal = .true.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual( 63.05, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(373.33 + 338.72, edn(nlev, 1, 1), atolerance, 'solar+thermal at surface :: downw flux not correct')
    @mpiassertEqual(130.91 + 384.44, eup(nlev, 1, 1), atolerance, 'solar+thermal at surface :: upward fl  not correct')
    @mpiassertEqual(1.332e-02-1.333e-02, abso(nlev-1,1,1), atolerance*1e-2, 'solar+thermal at surface :: absorption not correct')

    @mpiassertEqual(684.2709, edir(1, 1, 1), atolerance, 'solar+thermal at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar+thermal at TOA :: downw flux not correct')
    @mpiassertEqual(241.15 + 248.75, eup(1, 1, 1), atolerance, 'solar+thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(7.552e-08 -4.924e-16, abso(1, 1, 1), atolerance * 1e-2, 'solar+thermal at TOA :: absorption not correct')

    @mpiassertEqual(516.43, edir(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: direct flux not correct')
    @mpiassertEqual( 67.97, edir(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud+1 :: direct flux not correct')
    @mpiassertEqual(391.26 + 329.87, edn(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud :: downw flux not correct')
    @mpiassertEqual(224.00 + 381.32, eup(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0226 - 0.203, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar+thermal at icloud :: absorption not correct')

    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[1, 2])
  subroutine specint_ecckd_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this
    character(len=*), parameter :: specint = 'ecckd'

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation:'
    lthermal = .false.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, &
                       opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual( 93.13, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(344.54, edn(nlev, 1, 1), atolerance, 'solar at surface :: downw flux not correct')
    @mpiassertEqual(131.30, eup(nlev, 1, 1), atolerance, 'solar at surface :: upward fl  not correct')
    @mpiassertEqual(1.227e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'solar at surface :: absorption not correct')

    @mpiassertEqual(680.50, edir(1, 1, 1), atolerance, 'solar at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar at TOA :: downw flux not correct')
    @mpiassertEqual(227.62, eup(1, 1, 1), atolerance, 'solar at TOA :: upward fl  not correct')
    @mpiassertEqual(3.174e-10, abso(1, 1, 1), atolerance * 1e-2, 'solar at TOA :: absorption not correct')

    @mpiassertEqual(508.90, edir(nlev - icld, 1, 1), atolerance, 'solar at icloud :: direct flux not correct')
    @mpiassertEqual(100.95, edir(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: direct flux not correct')
    @mpiassertEqual(359.13, edn(nlev - icld + 1, 1, 1), atolerance, 'solar at icloud+1 :: downw flux not correct')
    @mpiassertEqual(210.88, eup(nlev - icld, 1, 1), atolerance, 'solar at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0310, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation:'
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done'

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual(329.91, edn(nlev, 1, 1), atolerance, 'thermal at surface :: downw flux not correct')
    @mpiassertEqual(384.06, eup(nlev, 1, 1), atolerance, 'thermal at surface :: upward fl  not correct')
    @mpiassertEqual(-1.589e-02, abso(nlev - 1, 1, 1), atolerance * 1e-2, 'thermal at surface :: absorption not correct')

    @mpiassertEqual(0.0, edn(1, 1, 1), atolerance, 'thermal at TOA :: downw flux not correct')
    @mpiassertEqual(243.80, eup(1, 1, 1), atolerance, 'thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(-1.402e-08, abso(1, 1, 1), atolerance * 1e-2, 'thermal at TOA :: absorption not correct')

    @mpiassertEqual(315.03, edn(nlev - icld + 1, 1, 1), atolerance, 'thermal at icloud :: downw flux not correct')
    @mpiassertEqual(378.14, eup(nlev - icld, 1, 1), atolerance, 'thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(-0.173, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'thermal at icloud :: absorption not correct')

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation:'
    lthermal = .true.; lsolar = .true.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    @mpiassertEqual( 93.13, edir(nlev, 1, 1), atolerance, 'solar at surface :: direct flux not correct')
    @mpiassertEqual(344.54 + 329.91, edn(nlev, 1, 1), atolerance, 'solar+thermal at surface :: downw flux not correct')
    @mpiassertEqual(131.30 + 384.06, eup(nlev, 1, 1), atolerance, 'solar+thermal at surface :: upward fl  not correct')
    @mpiassertEqual(1.227e-02-1.589e-02, abso(nlev-1,1,1), atolerance*1e-2, 'solar+thermal at surface :: absorption not correct')

    @mpiassertEqual(680.50, edir(1, 1, 1), atolerance, 'solar+thermal at TOA :: direct flux not correct')
    @mpiassertEqual(0, edn(1, 1, 1), atolerance, 'solar+thermal at TOA :: downw flux not correct')
    @mpiassertEqual(227.62 + 243.80, eup(1, 1, 1), atolerance, 'solar+thermal at TOA :: upward fl  not correct')
    @mpiassertEqual(3.174e-10 -1.402e-08, abso(1, 1, 1), atolerance * 1e-2, 'solar+thermal at TOA :: absorption not correct')

    @mpiassertEqual(508.90, edir(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: direct flux not correct')
    @mpiassertEqual(100.95, edir(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud+1 :: direct flux not correct')
    @mpiassertEqual(359.13 + 315.03, edn(nlev - icld + 1, 1, 1), atolerance, 'solar+thermal at icloud :: downw flux not correct')
    @mpiassertEqual(210.88 + 378.14, eup(nlev - icld, 1, 1), atolerance, 'solar+thermal at icloud :: upward fl  not correct')
    @mpiassertEqual(0.0310 - 0.173, abso(nlev - icld, 1, 1), atolerance * 1e-2, 'solar+thermal at icloud :: absorption not correct')

    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr); call CHKERR(ierr)
  end subroutine
end module
