module m_py_rrtm_lw_sw
  use m_data_parameters, only : ireals, mpiint, iintegers, &
    init_mpi_data_parameters, myid
  use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

  implicit none

  ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
  ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
  ! the size of the merged grid. If you only want to use heating rates on the
  ! dynamics grid, use the lower layers, i.e.,
  !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
  ! or:
  !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
  real(ireals), allocatable, dimension(:,:,:) :: edir,edn,eup,abso   ! [nlyr(+1), local_nx, local_ny ]

contains

  subroutine rrtmg(nlay, nxp, nyp, nprocx, nprocy, &
      comm, dx, dy, phi0, theta0,                     &
      albedo_thermal, albedo_solar, atm_filename,     &
      lthermal, lsolar,                               &
      d_plev, d_tlev, nxproc, nyproc,                 &
      d_tlay, d_h2ovmr, d_o3vmr,                      &
      d_co2vmr, d_ch4vmr, d_n2ovmr,  d_o2vmr,         &
      d_lwc, d_reliq, d_iwc, d_reice,                 &
      icollapse, opt_time, solar_albedo_2d)

    integer(iintegers) :: nlay, nxp, nyp, nprocx, nprocy
    integer(mpiint), intent(in) :: comm ! MPI Communicator

    real(ireals), intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=250), intent(in) :: atm_filename

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! dim(nlay_dynamics+1, nxp, nyp)
    real(ireals),intent(in) :: d_plev(nlay+1, nxp, nyp) ! pressure on layer interfaces [hPa]
    real(ireals),intent(in) :: d_tlev(nlay+1, nxp, nyp) ! Temperature on layer interfaces [K]

    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
    integer(iintegers),intent(in) :: nxproc(nprocx), nyproc(nprocy)

    ! all have dim(nlay_dynamics, nxp, nyp)
    real(ireals),intent(in) :: d_tlay   (nlay, nxp, nyp) ! layer mean temperature [K]
    real(ireals),intent(in) :: d_h2ovmr (nlay, nxp, nyp) ! watervapor volume mixing ratio [e.g. 1e-3]
    real(ireals),intent(in) :: d_o3vmr  (nlay, nxp, nyp) ! ozone volume mixing ratio      [e.g. .1e-6]
    real(ireals),intent(in) :: d_co2vmr (nlay, nxp, nyp) ! CO2 volume mixing ratio        [e.g. 407e-6]
    real(ireals),intent(in) :: d_ch4vmr (nlay, nxp, nyp) ! methane volume mixing ratio    [e.g. 2e-6]
    real(ireals),intent(in) :: d_n2ovmr (nlay, nxp, nyp) ! no2 volume mixing ratio        [e.g. .32]
    real(ireals),intent(in) :: d_o2vmr  (nlay, nxp, nyp) ! oxygen volume mixing ratio     [e.g. .2]
    real(ireals),intent(in) :: d_lwc    (nlay, nxp, nyp) ! liq water content              [g/kg]
    real(ireals),intent(in) :: d_reliq  (nlay, nxp, nyp) ! effective radius               [micron]
    real(ireals),intent(in) :: d_iwc    (nlay, nxp, nyp) ! ice water content              [g/kg]
    real(ireals),intent(in) :: d_reice  (nlay, nxp, nyp) ! ice effective radius           [micron]

    integer(iintegers),intent(in) :: icollapse ! experimental, dont use it if you dont know what you are doing otherwise set to 1

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo
    real(ireals), intent(in) :: opt_time, solar_albedo_2d(nxp, nyp)

    call init_mpi_data_parameters(comm)

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0,  &
      albedo_thermal, albedo_solar, atm_filename,     &
      lthermal, lsolar,                               &
      edir,edn,eup,abso,                              &
      d_plev, d_tlev, d_tlay, d_h2ovmr, d_o3vmr,      &
      d_co2vmr, d_ch4vmr, d_n2ovmr,  d_o2vmr,         &
      d_lwc, d_reliq, d_iwc, d_reice,                 &
      nxproc, nyproc, icollapse,                      &
      opt_time, solar_albedo_2d)

  end subroutine

  subroutine rrtmg_minimal(nlay, nxp, nyp, nprocx, nprocy, &
      comm, dx, dy, phi0, theta0,                          &
      albedo_thermal, albedo_solar, atm_filename,          &
      lthermal, lsolar,                                    &
      d_plev, d_tlev, d_lwc, d_reliq, d_iwc, d_reice,      &
      nxproc, nyproc)

    integer(iintegers) :: nlay, nxp, nyp, nprocx, nprocy
    integer(mpiint), intent(in) :: comm ! MPI Communicator

    real(ireals), intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=250), intent(in) :: atm_filename

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! dim(nlay_dynamics+1, nxp, nyp)
    real(ireals),intent(in) :: d_plev(nlay+1, nxp, nyp) ! pressure on layer interfaces [hPa]
    real(ireals),intent(in) :: d_tlev(nlay+1, nxp, nyp) ! Temperature on layer interfaces [K]

    ! all have dim(nlay_dynamics, nxp, nyp)
    real(ireals),intent(in) :: d_lwc    (nlay, nxp, nyp) ! liq water content              [g/kg]
    real(ireals),intent(in) :: d_reliq  (nlay, nxp, nyp) ! effective radius               [micron]
    real(ireals),intent(in) :: d_iwc    (nlay, nxp, nyp) ! ice water content              [g/kg]
    real(ireals),intent(in) :: d_reice  (nlay, nxp, nyp) ! ice effective radius           [micron]

    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
    integer(iintegers),intent(in) :: nxproc(nprocx), nyproc(nprocy)


    call init_mpi_data_parameters(comm)

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0,  &
      albedo_thermal, albedo_solar, atm_filename,     &
      lthermal, lsolar,                               &
      edir,edn,eup,abso,                              &
      d_plev=d_plev, d_tlev=d_tlev,                   &
      d_lwc=d_lwc, d_reliq=d_reliq,                   &
      d_iwc=d_iwc, d_reice=d_reice,                   &
      nxproc=nxproc, nyproc=nyproc)
  end subroutine

  subroutine destroy_rrtmg
    call destroy_tenstream_rrtmg()
    deallocate(edir, edn, eup, abso)
  end subroutine


end module
