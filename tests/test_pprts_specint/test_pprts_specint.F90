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

  use m_helper_functions, only: CHKERR, get_petsc_opt, spherical_2_cartesian, toStr

  use m_tenstream_options, only: read_commandline_options

  use m_netcdfio, only: ncwrite, ncload

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
  real(ireals), parameter :: atolerance = .1                 ! absolute tolerance when regression testing fluxes

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

  character(len=default_str_len) :: specint
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
    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr); call CHKERR(ierr)
    call destroy_tenstr_atm(atm)
  end subroutine teardown

  !@ generated result files with:
  !  PETSC_OPTIONS="-gen_test_results" mpirun -wd generated/test_pprts_specint/ $(pwd)/bin/pfunit_test_pprts_specint
  subroutine assert_results(this, fname, edir, edn, eup, abso, lthermal, lsolar, ierr)
    class(MpiTestMethod), intent(inout) :: this
    character(len=*), intent(in) :: fname
    real(ireals), allocatable, dimension(:, :, :), intent(in) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    logical, intent(in) :: lthermal, lsolar
    integer(mpiint), intent(out) :: ierr
    logical :: lgen_test_results, lflg
    character(len=default_str_len) :: groups(2)
    real(ireals), allocatable, dimension(:, :, :) :: tedir, tedn, teup, tabso ! [nlev_merged(-1), nxp, nyp]
    groups(1) = "result."//trim(fname)//".lw"//toStr(lthermal)//".sw"//toStr(lsolar)//".nc" ! test output data filename

    lgen_test_results = .false.
    call get_petsc_opt('', "-gen_test_results", lgen_test_results, lflg, ierr); call CHKERR(ierr)
    if (lgen_test_results .and. (numnodes .eq. 1)) then
      groups(2) = 'edir'; call ncwrite(groups, edir, ierr); call CHKERR(ierr)
      groups(2) = 'edn'; call ncwrite(groups, edn, ierr); call CHKERR(ierr)
      groups(2) = 'eup'; call ncwrite(groups, eup, ierr); call CHKERR(ierr)
      groups(2) = 'abso'; call ncwrite(groups, abso, ierr); call CHKERR(ierr)
    end if

    groups(2) = 'edir'; call ncload(groups, tedir, ierr); call CHKERR(ierr)
    groups(2) = 'edn'; call ncload(groups, tedn, ierr); call CHKERR(ierr)
    groups(2) = 'eup'; call ncload(groups, teup, ierr); call CHKERR(ierr)
    groups(2) = 'abso'; call ncload(groups, tabso, ierr); call CHKERR(ierr)

    @mpiassertEqual(tedir, edir, atolerance, 'edir not as expected')
    @mpiassertEqual(tedn, edn, atolerance, 'edn not as expected')
    @mpiassertEqual(teup, eup, atolerance, 'eup not as expected')
    @mpiassertEqual(tabso, abso, atolerance*1e-2, 'abso not as expected')
  end subroutine

  @test(npes=[1,2])
  subroutine specint_rrtm_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    specint = 'rrtmg'

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation: '//trim(specint)
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

    call assert_results(this, 'rrtm', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation: '//trim(specint)
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done '//trim(specint)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    call assert_results(this, 'rrtm', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation: '//trim(specint)
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

    call assert_results(this, 'rrtm', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[1,2])
  subroutine specint_repwvl_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    specint = 'repwvl'

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation: '//trim(specint)
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

    call assert_results(this, 'repwvl', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation: '//trim(specint)
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done '//trim(specint)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    call assert_results(this, 'repwvl', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation: '//trim(specint)
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

    call assert_results(this, 'repwvl', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)
  end subroutine

  @test(npes=[1,2])
  subroutine specint_ecckd_lw_sw(this)
    class(MpiTestMethod), intent(inout) :: this

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    specint = 'ecckd'

    ! For comparison, compute lw and sw separately
    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar Radiation: '//trim(specint)
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

    call assert_results(this, 'ecckd', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation: '//trim(specint)
    lthermal = .true.; lsolar = .false.

    call specint_pprts(specint, comm, solver, atm, nxp, nyp, &
                       dx, dy, spherical_2_cartesian(phi0, theta0), &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Thermal Radiation done '//trim(specint)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
        end do
      end if
    end if

    call assert_results(this, 'ecckd', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Computing Solar AND Thermal Radiation: '//trim(specint)
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

    call assert_results(this, 'ecckd', edir, edn, eup, abso, lthermal, lsolar, ierr); call CHKERR(ierr)
  end subroutine
end module
