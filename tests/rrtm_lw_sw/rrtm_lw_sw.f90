@test(npes =[2,1])
subroutine test_rrtm_lw(this)

    ! Import datatype from the TenStream lib. Depending on how PETSC is
    ! compiled(single or double floats, or long ints), this will determine what
    ! the Tenstream uses.
    use m_data_parameters, only : &
      init_mpi_data_parameters,   &
      iintegers, ireals, mpiint,  &
      zero, one, default_str_len

    ! main entry point for solver, and desctructor
    use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y

    integer(iintegers),parameter :: nxp=3, nyp=3, nzp=10 ! local domain size for each rank
    real(ireals),parameter :: dx=100, dy=dx              ! horizontal grid spacing in [m]
    real(ireals),parameter :: phi0=180, theta0=60        ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals),parameter :: albedo_th=0, albedo_sol=.3 ! broadband ground albedo for solar and thermal spectrum
    real(ireals),parameter :: atolerance = 1             ! absolute tolerance when regression testing fluxes

    real(ireals), dimension(nzp+1,nxp,nyp) :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,nxp,nyp) :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp,nxp,nyp) :: lwc, reliq

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: i,j,k, nlev, icld
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    logical,parameter :: ldebug=.True.
    logical :: lthermal, lsolar

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    N_ranks_y = sqrt(1.*numnodes)
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate(nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 50K
    do k=1,nzp+1
      plev(k,:,:) = 1000_ireals - (k-one)*500._ireals/(nzp)
      tlev(k,:,:) = 288._ireals - (k-one)*50._ireals/(nzp)
    enddo

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    h2ovmr = zero
    o3vmr  = zero
    co2vmr = zero
    ch4vmr = zero
    n2ovmr = zero
    o2vmr  = zero

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    icld = (nzp+1)/2
    lwc  (icld, :,:) = 1e-2
    reliq(icld, :,:) = 10

    tlev (icld  , :,:) = 288
    tlev (icld+1, :,:) = tlev (icld  , :,:)

    ! mean layer temperature is approx. arithmetic mean between layer interfaces
    tlay = (tlev(1:nzp,:,:) + tlev(2:nzp+1,:,:))/2

    ! For comparison, compute lw and sw separately
    if(myid.eq.0 .and. ldebug) print *,'Computing Solar Radiation:'
    lthermal=.False.; lsolar=.True.

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0, albedo_th, albedo_sol, &
      atm_filename, lthermal, lsolar,                                       &
      edir,edn,eup,abso,                                                    &
      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,      &
      nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(313.47, edir(nlev,1,1), atolerance, 'solar at surface :: direct flux not correct')
        @assertEqual(143.32, edn (nlev,1,1), atolerance, 'solar at surface :: downw flux not correct')
        @assertEqual(137.04, eup (nlev,1,1), atolerance, 'solar at surface :: upward fl  not correct')
        @assertEqual(-1.395E-02, abso(nlev-1,1,1), atolerance, 'solar at surface :: absorption not correct')

        @assertEqual(684.1109, edir(1,1,1), atolerance, 'solar at TOA :: direct flux not correct')
        @assertEqual(0       , edn (1,1,1), atolerance, 'solar at TOA :: downw flux not correct')
        @assertEqual(207.18  , eup (1,1,1), atolerance, 'solar at TOA :: upward fl  not correct')
        @assertEqual(-2.063E-04, abso(1,1,1), atolerance, 'solar at TOA :: absorption not correct')

        @assertEqual(502.23 , edir(nlev-icld  ,1,1), atolerance, 'solar at icloud :: direct flux not correct')
        @assertEqual(339.22 , edir(nlev-icld+1,1,1), atolerance, 'solar at icloud+1 :: direct flux not correct')
        @assertEqual(143.68 , edn (nlev-icld+1,1,1), atolerance, 'solar at icloud :: downw flux not correct')
        @assertEqual(190.29 , eup (nlev-icld  ,1,1), atolerance, 'solar at icloud :: upward fl  not correct')
        @assertEqual(-0.0242, abso(nlev-icld  ,1,1), atolerance, 'solar at icloud :: absorption not correct')
    endif


    if(myid.eq.0 .and. ldebug) print *,'Computing Thermal Radiation:'
    lthermal=.True.; lsolar=.False.

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0, albedo_th, albedo_sol, &
      atm_filename, lthermal, lsolar,                                       &
      edir,edn,eup,abso,                                                    &
      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,      &
      nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    if(myid.eq.0 .and. ldebug) print *,'Computing Thermal Radiation done'

    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(325.29, edn (nlev,1,1), atolerance, 'thermal at surface :: downw flux not correct')
        @assertEqual(390.07, eup (nlev,1,1), atolerance, 'thermal at surface :: upward fl  not correct')
        @assertEqual(-2.321E-02, abso(nlev-1,1,1), atolerance, 'thermal at surface :: absorption not correct')

        @assertEqual(  0.0     , edn (1,1,1), atolerance, 'thermal at TOA :: downw flux not correct')
        @assertEqual(254.90    , eup (1,1,1), atolerance, 'thermal at TOA :: upward fl  not correct')
        @assertEqual(-2.193E-04, abso(1,1,1), atolerance, 'thermal at TOA :: absorption not correct')

        @assertEqual(318.18, edn (nlev-icld+1,1,1), atolerance, 'thermal at icloud :: downw flux not correct')
        @assertEqual(386.69, eup (nlev-icld  ,1,1), atolerance, 'thermal at icloud :: upward fl  not correct')
        @assertEqual(-0.255, abso(nlev-icld  ,1,1), atolerance, 'thermal at icloud :: absorption not correct')
    endif

    if(myid.eq.0 .and. ldebug) print *,'Computing Solar AND Thermal Radiation:'
    lthermal=.True.; lsolar=.True.

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0, albedo_th, albedo_sol, &
      atm_filename, lthermal, lsolar,                                       &
      edir,edn,eup,abso,                                                    &
      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,      &
      nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(313.47, edir(nlev,1,1), atolerance, 'solar+thermal at surface :: direct flux not correct')
        @assertEqual(468.61, edn (nlev,1,1), atolerance, 'solar+thermal at surface :: downw flux not correct')
        @assertEqual(527.11, eup (nlev,1,1), atolerance, 'solar+thermal at surface :: upward fl  not correct')
        @assertEqual(-0.03716, abso(nlev-1,1,1), atolerance, 'solar+thermal at surface :: absorption not correct')

        @assertEqual(684.1109, edir(1,1,1), atolerance, 'solar+thermal at TOA :: direct flux not correct')
        @assertEqual(0       , edn (1,1,1), atolerance, 'solar+thermal at TOA :: downw flux not correct')
        @assertEqual(462.08  , eup (1,1,1), atolerance, 'solar+thermal at TOA :: upward fl  not correct')
        @assertEqual(-4.4256e-4, abso(1,1,1), atolerance, 'solar+thermal at TOA :: absorption not correct')

        @assertEqual(502.23 , edir(nlev-icld  ,1,1), atolerance, 'solar+thermal at icloud :: direct flux not correct')
        @assertEqual(339.22 , edir(nlev-icld+1,1,1), atolerance, 'solar+thermal at icloud+1 :: direct flux not correct')
        @assertEqual(461.86 , edn (nlev-icld+1,1,1), atolerance, 'solar+thermal at icloud :: downw flux not correct')
        @assertEqual(576.98 , eup (nlev-icld  ,1,1), atolerance, 'solar+thermal at icloud :: upward fl  not correct')
        @assertEqual(-0.2792, abso(nlev-icld  ,1,1), atolerance, 'solar+thermal at icloud :: absorption not correct')
    endif

    ! Tidy up
    call destroy_tenstream_rrtmg(lfinalizepetsc=.True.)

end subroutine
