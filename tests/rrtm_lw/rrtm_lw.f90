@test(npes =[3,2,1])
subroutine test_rrtm_lw(this)

    use m_data_parameters, only : init_mpi_data_parameters, &
      iintegers, ireals, mpiint, zero, one, default_str_len

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast

    use m_tenstr_rrtm_lw, only : tenstream_rrtm_lw, destroy_tenstream_rrtm_lw

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: nxp=6, nyp=3, nzp=10 ! local domain size for each rank
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: albedo=0, phi0=180, theta0=0
    real(ireals),parameter :: atolerance = 1

    real(ireals), dimension(nzp+1,nxp,nyp) :: plev, tlev
    real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
    real(ireals), dimension(nzp,nxp,nyp) :: lwc, reliq
    real(ireals),allocatable, dimension(:,:,:) :: edn,eup,abso ! [nlev_merged(-1), nxp, nyp]

    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    integer(iintegers) :: k, nlev, icld
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    logical,parameter :: ldebug=.True.

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    allocate(nxproc(numnodes), source=nxp)   ! domain decomp in x direction with all ranks
    allocate(nyproc(1), source=nyp)          ! no domain decomposition in y direction

    call init_mpi_data_parameters(comm)

    do k=1,nzp+1
      plev(k,:,:) = 1000_ireals - (k-one)*500._ireals/(nzp)
      tlev(k,:,:) = 288._ireals - (k-one)*50._ireals/(nzp)
    enddo

    tlay = (tlev(1:nzp,:,:) + tlev(2:nzp+1,:,:))/2

    h2ovmr = zero
    o3vmr  = zero
    co2vmr = zero
    ch4vmr = zero
    n2ovmr = zero
    o2vmr  = zero

    lwc = 0
    reliq = 0

    icld = (nzp+1)/2
    lwc  (icld, :,:) = .1
    reliq(icld, :,:) = 10
    tlev (icld  , :,:) = 288
    tlev (icld+1, :,:) = tlev (icld  , :,:)

    call tenstream_rrtm_lw(comm, dx, dy, phi0, theta0, albedo, atm_filename, &
                               edn,eup,abso,                                 &
                               plev, tlev, d_lwc=lwc, d_reliq=reliq,         &
                               nxproc=nxproc, nyproc=nyproc)
    call destroy_tenstream_rrtm_lw()

    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(377.6775, edn (nlev,1,1), atolerance, 'thermal at surface :: downw flux not correct')
        @assertEqual(390.0707, eup (nlev,1,1), atolerance, 'thermal at surface :: upward fl  not correct')
        @assertEqual(-1.443E-02, abso(nlev-1,1,1), atolerance, 'thermal at surface :: absorption not correct')

        @assertEqual(0             , edn (1,1,1), atolerance, 'thermal at TOA :: downw flux not correct')
        @assertEqual(256.4602      , eup (1,1,1), atolerance, 'thermal at TOA :: upward fl  not correct')
        @assertEqual(-2.1898020E-04, abso(1,1,1), atolerance, 'thermal at TOA :: absorption not correct')

        @assertEqual(389.5942  , edn (nlev-icld+1,1,1), atolerance, 'thermal at icloud :: downw flux not correct')
        @assertEqual(390.0449  , eup (nlev-icld  ,1,1), atolerance, 'thermal at icloud :: upward fl  not correct')
        @assertEqual(-0.3934462, abso(nlev-icld  ,1,1), atolerance, 'thermal at icloud :: absorption not correct')
    endif

end subroutine
