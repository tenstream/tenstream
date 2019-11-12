module m_py_rrtm_lw_sw
  use m_data_parameters, only : ireals, mpiint, iintegers, init_mpi_data_parameters, default_str_len
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg
  use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_3_16
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  implicit none

  private
  public :: edir,edn,eup,abso, rrtmg, destroy_rrtmg

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
      d_lwc, d_reliq )

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
    !double precision,intent(in) :: d_tlay   (nlay, nxp, nyp) ! layer mean temperature [K]
    !double precision,intent(in) :: d_h2ovmr (nlay, nxp, nyp) ! watervapor volume mixing ratio [e.g. 1e-3]
    !double precision,intent(in) :: d_o3vmr  (nlay, nxp, nyp) ! ozone volume mixing ratio      [e.g. .1e-6]
    !double precision,intent(in) :: d_co2vmr (nlay, nxp, nyp) ! CO2 volume mixing ratio        [e.g. 407e-6]
    !double precision,intent(in) :: d_ch4vmr (nlay, nxp, nyp) ! methane volume mixing ratio    [e.g. 2e-6]
    !double precision,intent(in) :: d_n2ovmr (nlay, nxp, nyp) ! no2 volume mixing ratio        [e.g. .32]
    !double precision,intent(in) :: d_o2vmr  (nlay, nxp, nyp) ! oxygen volume mixing ratio     [e.g. .2]
    double precision,intent(in) :: d_lwc    (nlay, nxp, nyp) ! liq water content              [g/kg]
    double precision,intent(in) :: d_reliq  (nlay, nxp, nyp) ! effective radius               [micron]
    !double precision,intent(in) :: d_iwc    (nlay, nxp, nyp) ! ice water content              [g/kg]
    !double precision,intent(in) :: d_reice  (nlay, nxp, nyp) ! ice effective radius           [micron]

    real(ireals), allocatable, dimension(:,:,:) :: edir_ir,edn_ir,eup_ir,abso_ir ! [nlyr(+1), local_nx, local_ny ]

    integer(iintegers) :: icol, Ncol, i,j
    real(ireals), allocatable, dimension(:,:) :: plev, tlev, lwc, reliq
    character(default_str_len) :: atmfname

    type(t_solver_3_16) :: pprts_solver
    type(t_tenstr_atm) :: atm

    atmfname = atm_filename

    call init_mpi_data_parameters(comm)

    Ncol = nxp*nyp

    allocate(plev(nlay+1,Ncol))
    allocate(tlev(nlay+1,Ncol))
    allocate(lwc(nlay   ,Ncol))
    allocate(reliq(nlay ,Ncol))

    do j=1,nyp
      do i=1,nxp
        icol = i + (j-1)*nxp
        plev (:, icol) = real(d_plev (:,i,j), ireals)
        tlev (:, icol) = real(d_tlev (:,i,j), ireals)
        lwc  (:, icol) = real(d_lwc  (:,i,j), ireals)
        reliq(:, icol) = real(d_reliq(:,i,j), ireals)
      enddo
    enddo

    call setup_tenstr_atm(comm, .False., atmfname, &
      plev, tlev, atm, &
      d_lwc=lwc, d_reliq=reliq)

    call pprts_rrtmg(comm, pprts_solver, atm, nxp, nyp,                 &
      real(dx,kind=ireals), real(dy,kind=ireals),                       &
      real(phi0,kind=ireals), real(theta0,kind=ireals),                 &
      real(albedo_thermal,kind=ireals), real(albedo_solar,kind=ireals), &
      lthermal, lsolar,                                                 &
      edir_ir,edn_ir,eup_ir,abso_ir,                                    &
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

    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.True.)
  end subroutine

  subroutine destroy_rrtmg
    deallocate(edir, edn, eup, abso)
  end subroutine
end module
