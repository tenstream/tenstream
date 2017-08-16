@test(npes =[2,1])
subroutine test_rrtm_sw(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast, CHKERR

    use m_tenstr_rrtm_sw, only : tenstream_rrtm_sw

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: nxp=6,nyp=3
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0.2
    real(ireals),parameter :: atolerance = 1

    real(ireals),allocatable :: atm(:,:) ! # z(km)  p(mb)  T(K) air(cm-3) o3(cm-3) o2(cm-3)  h2o(cm-3) co2(cm-3) no2(cm-3)
    real(ireals),allocatable,dimension(:,:,:) :: plev                                               ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq                                         ! nlay  , nxp, nyp
    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso          ! [nlyr(+1), global_nx, global_ny ]

    integer(iintegers) :: i,j,k, nlay

    logical,parameter :: ldebug=.True.

    integer(mpiint) :: ierr

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    if(myid.eq.0) then
        call read_ascii_file_2d('afglus_100m.dat', atm, 9, 2, ierr); call CHKERR(ierr)

        nlay = ubound(atm,1)-1

        allocate(plev   (nlay+1, nxp, nyp))
        allocate(tlay   (nlay  , nxp, nyp))
        allocate(h2ovmr (nlay  , nxp, nyp))
        allocate(o3vmr  (nlay  , nxp, nyp))
        allocate(co2vmr (nlay  , nxp, nyp))
        allocate(ch4vmr (nlay  , nxp, nyp))
        allocate(n2ovmr (nlay  , nxp, nyp))
        allocate(o2vmr  (nlay  , nxp, nyp))

        do j=1,nyp
            do i=1,nxp
                plev(:,i,j) = atm(:,2)
                tlay(:,i,j) = meanvec(atm(:,3))
                h2ovmr(:,i,j) = meanvec(atm(:,7)) / meanvec(atm(:,4))
                o3vmr (:,i,j) = meanvec(atm(:,5)) / meanvec(atm(:,4))
                co2vmr(:,i,j) = meanvec(atm(:,8)) / meanvec(atm(:,4))
                ch4vmr(:,i,j) = co2vmr(:,i,j)/1e2
                n2ovmr(:,i,j) = meanvec(atm(:,9)) / meanvec(atm(:,4))
                o2vmr (:,i,j) = meanvec(atm(:,6)) / meanvec(atm(:,4))
            enddo
        enddo

        if(ldebug) then
            print *,'loaded',myid
            if(myid.eq.0) print *,'plev   shape',shape(plev)
            if(myid.eq.0) print *,'h2ovmr shape',shape(h2ovmr)
            if(myid.eq.0) print *,'o3vmr  shape',shape(o3vmr )
            if(myid.eq.0) print *,'co2vmr shape',shape(co2vmr)
            if(myid.eq.0) print *,'n2ovmr shape',shape(n2ovmr)
            if(myid.eq.0) print *,'o2vmr  shape',shape(o2vmr )
        endif
    endif
    call imp_bcast(comm, nlay  , 0_mpiint)
    call imp_bcast(comm, plev  , 0_mpiint)
    call imp_bcast(comm, tlay  , 0_mpiint)
    call imp_bcast(comm, h2ovmr, 0_mpiint)
    call imp_bcast(comm, o3vmr , 0_mpiint)
    call imp_bcast(comm, co2vmr, 0_mpiint)
    call imp_bcast(comm, ch4vmr, 0_mpiint)
    call imp_bcast(comm, n2ovmr, 0_mpiint)
    call imp_bcast(comm, o2vmr , 0_mpiint)

    allocate(lwc   (nlay  ,nxp,nyp))
    allocate(reliq (nlay  ,nxp,nyp))
    lwc = 0

    !icld = minloc(abs(atm(:,1)-5.5),dim=1)
    !lwc(icld:icld+5, nxp/2, nyp/2) = .1
    reliq = 10

    if(myid.eq.0 .and. ldebug) then
        do k=1,nlay
            print *,'plev',plev(k,1,1), 'T', tlay(k,1,1),'CO2',co2vmr(k,1,1),'H2O',h2ovmr(k,1,1)
        enddo
        print *,'plev',plev(nlay+1,1,1)
    endif
    call tenstream_rrtm_sw(comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir, edn, eup, abso)

    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlay+1
                print *,k,'edir', edir(k,1,1), edn(k,1,1), eup(k,1,1), abso(min(nlay,k),1,1)
            enddo
        endif

        @assertEqual(452.7361, edir(nlay+1,1,1), atolerance, 'solar at surface :: divergence not correct')
        @assertEqual(49.50430, edn (nlay+1,1,1), atolerance, 'solar at surface :: downw flux not correct')
        @assertEqual(100.4471, eup (nlay+1,1,1), atolerance, 'solar at surface :: upward fl  not correct')
        @assertEqual(.0159963, abso(nlay  ,1,1), atolerance, 'solar at surface :: absorption not correct')

        @assertEqual(684.1108, edir(1,1,1), atolerance, 'solar at TOA :: divergence not correct')
        @assertEqual(  0.0000, edn (1,1,1), atolerance, 'solar at TOA :: downw flux not correct')
        @assertEqual(131.6273, eup (1,1,1), atolerance, 'solar at TOA :: upward fl  not correct')
        @assertEqual(1.905E-4, abso(1,1,1), atolerance, 'solar at TOA :: absorption not correct')

    endif

end subroutine
