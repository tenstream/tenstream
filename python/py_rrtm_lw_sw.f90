module m_py_rrtm_lw_sw
  use m_data_parameters, only : ireals, mpiint, iintegers, init_mpi_data_parameters, myid, default_str_len
  use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

  implicit none

  ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
  ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
  ! the size of the merged grid. If you only want to use heating rates on the
  ! dynamics grid, use the lower layers, i.e.,
  !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
  ! or:
  !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
  double precision, allocatable, dimension(:,:,:) :: edir,edn,eup,abso   ! [nlyr(+1), local_nx, local_ny ]

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

    double precision, intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    double precision, intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    double precision, intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=*), intent(in) :: atm_filename

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! dim(nlay_dynamics+1, nxp, nyp)
    double precision,intent(in) :: d_plev(nlay+1, nxp, nyp) ! pressure on layer interfaces [hPa]
    double precision,intent(in) :: d_tlev(nlay+1, nxp, nyp) ! Temperature on layer interfaces [K]

    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
    integer(iintegers),intent(in) :: nxproc(nprocx), nyproc(nprocy)

    ! all have dim(nlay_dynamics, nxp, nyp)
    double precision,intent(in) :: d_tlay   (nlay, nxp, nyp) ! layer mean temperature [K]
    double precision,intent(in) :: d_h2ovmr (nlay, nxp, nyp) ! watervapor volume mixing ratio [e.g. 1e-3]
    double precision,intent(in) :: d_o3vmr  (nlay, nxp, nyp) ! ozone volume mixing ratio      [e.g. .1e-6]
    double precision,intent(in) :: d_co2vmr (nlay, nxp, nyp) ! CO2 volume mixing ratio        [e.g. 407e-6]
    double precision,intent(in) :: d_ch4vmr (nlay, nxp, nyp) ! methane volume mixing ratio    [e.g. 2e-6]
    double precision,intent(in) :: d_n2ovmr (nlay, nxp, nyp) ! no2 volume mixing ratio        [e.g. .32]
    double precision,intent(in) :: d_o2vmr  (nlay, nxp, nyp) ! oxygen volume mixing ratio     [e.g. .2]
    double precision,intent(in) :: d_lwc    (nlay, nxp, nyp) ! liq water content              [g/kg]
    double precision,intent(in) :: d_reliq  (nlay, nxp, nyp) ! effective radius               [micron]
    double precision,intent(in) :: d_iwc    (nlay, nxp, nyp) ! ice water content              [g/kg]
    double precision,intent(in) :: d_reice  (nlay, nxp, nyp) ! ice effective radius           [micron]

    integer(iintegers),intent(in) :: icollapse ! experimental, dont use it if you dont know what you are doing otherwise set to 1

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo
    double precision, intent(in) :: opt_time, solar_albedo_2d(nxp, nyp)

    real(ireals), allocatable, dimension(:,:,:) :: edir_ir,edn_ir,eup_ir,abso_ir ! [nlyr(+1), local_nx, local_ny ]

    character(default_str_len) :: atmfname
    atmfname = atm_filename

    call init_mpi_data_parameters(comm)

    call tenstream_rrtmg(comm, real(dx,kind=ireals), real(dy,kind=ireals), &
      real(phi0,kind=ireals), real(theta0,kind=ireals),                    &
      real(albedo_thermal,kind=ireals), real(albedo_solar,kind=ireals),    &
      atmfname,  lthermal, lsolar,                                         &
      edir_ir,edn_ir,eup_ir,abso_ir,                                       &
      real(d_plev,kind=ireals), real(d_tlev,kind=ireals),                  &
      real(d_tlay,kind=ireals), real(d_h2ovmr,kind=ireals),                &
      real(d_o3vmr,kind=ireals),                                           &
      real(d_co2vmr,kind=ireals), real(d_ch4vmr,kind=ireals),              &
      real(d_n2ovmr,kind=ireals),  real(d_o2vmr,kind=ireals),              &
      real(d_lwc,kind=ireals), real(d_reliq,kind=ireals),                  &
      real(d_iwc,kind=ireals), real(d_reice,kind=ireals),                  &
      nxproc, nyproc, icollapse,                                           &
      real(opt_time,kind=ireals), real(solar_albedo_2d,kind=ireals))

    if(allocated(edir_ir)) then
      if(allocated(edir)) deallocate(edir)
      allocate(edir(size(edir_ir,dim=1),size(edir_ir,dim=2),size(edir_ir,dim=3)))
      edir = edir_ir
    endif
    if(allocated(edn_ir)) then
      if(allocated(edn)) deallocate(edn)
      allocate(edn(size(edn_ir,dim=1),size(edn_ir,dim=2),size(edn_ir,dim=3)))
      edn = edn_ir
    endif
    if(allocated(eup_ir)) then
      if(allocated(eup)) deallocate(eup)
      allocate(eup(size(eup_ir,dim=1),size(eup_ir,dim=2),size(eup_ir,dim=3)))
      eup = eup_ir
    endif
    if(allocated(abso_ir)) then
      if(allocated(abso)) deallocate(abso)
      allocate(abso(size(abso_ir,dim=1),size(abso_ir,dim=2),size(abso_ir,dim=3)))
      abso = abso_ir
    endif
  end subroutine

  subroutine rrtmg_minimal(nlay, nxp, nyp, nprocx, nprocy, &
      comm, dx, dy, phi0, theta0,                          &
      albedo_thermal, albedo_solar, atm_filename,          &
      lthermal, lsolar,                                    &
      d_plev, d_tlev, d_lwc, d_reliq, d_iwc, d_reice,      &
      nxproc, nyproc)

    integer(iintegers) :: nlay, nxp, nyp, nprocx, nprocy
    integer(mpiint), intent(in) :: comm ! MPI Communicator

    double precision, intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    double precision, intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    double precision, intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=*), intent(in) :: atm_filename

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! dim(nlay_dynamics+1, nxp, nyp)
    double precision,intent(in) :: d_plev(nlay+1, nxp, nyp) ! pressure on layer interfaces [hPa]
    double precision,intent(in) :: d_tlev(nlay+1, nxp, nyp) ! Temperature on layer interfaces [K]

    ! all have dim(nlay_dynamics, nxp, nyp)
    double precision,intent(in) :: d_lwc    (nlay, nxp, nyp) ! liq water content              [g/kg]
    double precision,intent(in) :: d_reliq  (nlay, nxp, nyp) ! effective radius               [micron]
    double precision,intent(in) :: d_iwc    (nlay, nxp, nyp) ! ice water content              [g/kg]
    double precision,intent(in) :: d_reice  (nlay, nxp, nyp) ! ice effective radius           [micron]

    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
    integer(iintegers),intent(in) :: nxproc(nprocx), nyproc(nprocy)

    real(ireals), allocatable, dimension(:,:,:) :: edir_ir,edn_ir,eup_ir,abso_ir ! [nlyr(+1), local_nx, local_ny ]

    character(default_str_len) :: atmfname
    atmfname = atm_filename

    call init_mpi_data_parameters(comm)

    call tenstream_rrtmg(comm, real(dx,kind=ireals), real(dy,kind=ireals), &
      real(phi0,kind=ireals), real(theta0,kind=ireals),                    &
      real(albedo_thermal,kind=ireals), real(albedo_solar,kind=ireals),    &
      atmfname,  lthermal, lsolar,                                         &
      edir_ir,edn_ir,eup_ir,abso_ir,                                       &
      d_plev=real(d_plev,kind=ireals), d_tlev=real(d_tlev,kind=ireals),    &
      d_lwc=real(d_lwc,kind=ireals), d_reliq=real(d_reliq,kind=ireals),    &
      d_iwc=real(d_iwc,kind=ireals), d_reice=real(d_reice,kind=ireals),    &
      nxproc=nxproc, nyproc=nyproc)

    if(allocated(edir_ir)) then
      if(allocated(edir)) deallocate(edir)
      allocate(edir(size(edir_ir,dim=1),size(edir_ir,dim=2),size(edir_ir,dim=3)))
      edir = edir_ir
    endif
    if(allocated(edn_ir)) then
      if(allocated(edn)) deallocate(edn)
      allocate(edn(size(edn_ir,dim=1),size(edn_ir,dim=2),size(edn_ir,dim=3)))
      edn = edn_ir
    endif
    if(allocated(eup_ir)) then
      if(allocated(eup)) deallocate(eup)
      allocate(eup(size(eup_ir,dim=1),size(eup_ir,dim=2),size(eup_ir,dim=3)))
      eup = eup_ir
    endif
    if(allocated(abso_ir)) then
      if(allocated(abso)) deallocate(abso)
      allocate(abso(size(abso_ir,dim=1),size(abso_ir,dim=2),size(abso_ir,dim=3)))
      abso = abso_ir
    endif
  end subroutine

  subroutine destroy_rrtmg
    call destroy_tenstream_rrtmg(lfinalizepetsc=.True.)
    deallocate(edir, edn, eup, abso)
  end subroutine


end module
