module m_tenstr_disort
  use m_data_parameters, only: mpiint
  use m_search, only: find_real_location
  use m_helper_functions, only: CHKERR

  implicit none

  private
  public :: default_flx_computation

contains
  subroutine default_flx_computation(&
      mu0, S0, Ag, &
      lthermal, wvnm, Bfracs, &
      dtau, ssalb, gasym, temper, &
      RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
      nstreams, lverbose)

    real, intent(in)    :: mu0   ! cos(solar zenith angle)
    real, intent(in)    :: S0    ! solar constant
    real, intent(in)    :: Ag    ! lambertian surface albedo
    logical, intent(in) :: lthermal ! do thermal computations ?
    real, intent(in)    :: wvnm(2) ! Wavenumbers low and high [inv cm], ignored if not lthermal
    real, dimension(:), intent(in)  :: Bfracs      ! Planck fractions for vertical weighting (needed for RRTMG), otherwise use 1 (nlay)
    real, dimension(:), intent(in)  :: dtau        ! vertical optical thicknesses (nlay)
    real, dimension(:), intent(in)  :: ssalb       ! single scatter albedo        (nlay)
    real, dimension(:), intent(in)  :: gasym       ! asymmetry parameter          (nlay)
    real, dimension(:), intent(in)  :: temper      ! level temperatures           (nlay+1)
    real, dimension(:), intent(out) :: RFLDIR      ! Direct-beam flux (without delta-M scaling)
    real, dimension(:), intent(out) :: RFLDN       ! Diffuse down-flux (total minus direct-beam) (without delta-M scaling)
    real, dimension(:), intent(out) :: FLUP        ! Diffuse up-flux
    real, dimension(:), intent(out) :: DFDT        ! Flux divergence  d(net flux)/d(optical depth)
    real, dimension(:), intent(out) :: UAVG        ! Mean intensity (including the direct beam)
    integer, optional, intent(in)   :: nstreams    ! number of streams to use (default: 16)
    logical, optional, intent(in)   :: lverbose    ! verbose output of disort (default: False)

    integer :: nlyr, ntau, nstr, nmom, numu

    real, allocatable :: utau(:)
    real, allocatable :: pmom(:,:)
    real, allocatable :: bemst(:)            ! Surface directional emissivity (computational angles)
    real, allocatable :: emust(:)            ! Surface directional emissivity (user angles)
    real, allocatable :: h_lyr(:)            ! spherical layer height
    real, allocatable :: UMU(:)              ! cosines of output polar angles
    real, allocatable :: UU(:,:,:)           ! Intensity ( if ONLYFL = FALSE;  zero otherwise )
    real, allocatable :: ALBMED(:)           ! Albedo of the medium as a function of incident beam angle cosine
    ! UMU(IU)  (IBCND = 1 case only)
    real, allocatable :: TRNMED(:)           ! Transmissivity of the medium as a function of incident
    ! beam angle cosine UMU(IU)  (IBCND = 1 case only)
    real, allocatable :: rho_accurate(:,:)   ! analytic brdf results
    real, allocatable, dimension(:,:,:) :: & ! brdf fourier components matrix
      rhoq, rhou


    ! Default vars will be overridden in beneath code:
    logical :: usrtau            ! output only at levels, not at user defined taus
    logical :: usrang            ! output only at computational polar angles
    integer, parameter :: nphi=1 ! number of azimuthal angles of return intensities, zero only valid if ONLYFL
    real    :: phi(nphi)
    integer :: ibcnd             ! general boundary condition case
    real    :: fbeam             ! intensity of incident beam
    real    :: umu0              ! Polar angle cosine of incident beam
    real    :: phi0              ! Azimuth angle of incident beam (0 to 360 degrees)
    real    :: fisot             ! Intensity of top-boundary isotropic illumination, zero for outer space
    logical :: lamber            ! isotropically reflecting bottom boundary
    real    :: albedo            ! isotropic albedo value
    real    :: btemp             ! Temperature of bottom boundary [K]
    real    :: ttemp             ! Temperature of top boundary [K]
    real    :: temis             ! Emissivity of top boundary

    logical :: deltaM            ! use delta-M method
    logical :: pseudo_sphere     ! employ pseudo spherical computations
    real    :: earth_radius

    logical :: plank             ! include thermal emission
    real    :: wvnmlo, wvnmhi    ! Wavenumbers [inv cm] of spectral interval of interest
    logical :: onlyfl            ! return fluxes, flux divergences, and mean intensities
    real    :: accur             ! Convergence criterion for azimuthal series
    logical :: prnt(7)           ! debug output level
    integer :: maxcly            ! Max. number of computational layers
    integer :: maxulv            ! Max. number of user levels

    integer :: maxumu            ! Max. number of user polar angles
    integer :: maxcmu            ! Max. number of computaional polar angles
    integer :: maxphi            ! Max. number of user azimuth angles
    character(len=127) HEADER


    nstr = 16
    if(present(nstreams)) nstr = nstreams


    prnt = .False.
    if(present(lverbose)) then
      prnt(2) = lverbose
    endif

    nmom = nstr+2
    numu = nstr

    nlyr = size(dtau)
    ntau = nlyr+1

    fbeam  = S0
    umu0   = mu0
    albedo = Ag

    call set_default_values()

    call DISORT_rrtmg( nlyr, nmom, numu, &
      MAXUMU, MAXPHI, MAXULV, &
      USRANG, USRTAU, IBCND, ONLYFL, PRNT, &
      PLANK, Bfracs, LAMBER, DELTAM, PSEUDO_SPHERE, &
      DTAU, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI, &
      UTAU, UMU0, PHI0, UMU, PHI, FBEAM, &
      FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, &
      EARTH_RADIUS, H_LYR, &
      RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST, &
      ACCUR,  trim(HEADER), &
      RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, &
      ALBMED, TRNMED )

    call check_bad_mu0()

    contains
      subroutine check_bad_mu0()
        integer :: loc
        loc = nint(find_real_location(umu, umu0))
        if(abs(umu(loc)-umu0).le.epsilon(umu)*100) then
          print *,'our choice of umu might be bad. '// &
            'You could try to use a different number of disort streams. '// &
            'Try to set a value where nstr/2 is an even number.'
          print *,'sampling points: umu', umu
          print *,'closest umu sampling point:', umu(loc)
          print *,'solar mu0', umu0
          call CHKERR(1_mpiint, 'bad disort umu sampling points')
        endif

      end subroutine
      subroutine set_default_values()
        allocate( utau(ntau) )
        allocate( pmom(0:nmom,nlyr) )
        allocate( bemst(nstr/2) )
        allocate( emust(numu) )
        allocate( h_lyr(0:nlyr) )
        allocate( UMU(numu) )
        allocate( UU(numu,ntau,nphi) )
        allocate( ALBMED(numu) )
        allocate( TRNMED(numu) )
        allocate( rho_accurate(numu,nphi) )
        allocate( rhoq(nstr/2,0:nstr/2,0:(nstr-1)) )
        allocate( rhou(nstr/2,0:nstr/2,0:(nstr-1)) )

        pmom(0,:) = 1
        pmom(1,:) = gasym
        pmom(2:nmom,:) = 0

        usrtau=.False.
        utau = 0
        usrang=.False.

        phi = 0
        ibcnd=0

        phi0   = 0
        fisot=0
        lamber=.True.

        btemp  = temper(size(temper))
        ttemp  = 0
        temis  = 0
        bemst  = 0
        emust  = 0
        deltaM = .True.

        pseudo_sphere = .False.
        earth_radius = 6371.0
        h_lyr = 0

        plank  = lthermal
        wvnmlo = wvnm(1)
        wvnmhi = wvnm(2)

        onlyfl = .True.
        accur  = 0

        header = ''

        maxcly = nlyr
        maxulv = nlyr+1
        maxumu = nstr
        maxcmu = nstr
        maxphi = nphi
      end subroutine
  end subroutine
end module
