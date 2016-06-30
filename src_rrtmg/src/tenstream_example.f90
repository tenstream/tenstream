program main

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,imp_comm,myid,mpierr,zero, i0, numnodes
    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast

    use m_netcdfIO, only : ncwrite, ncload

    use m_tenstr_rrtm_sw
    implicit none

    !integer(iintegers),parameter :: nxp=32, nyp=nxp*2
    integer(iintegers) :: nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: plev                                                     ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, air  ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq                                               ! nlay  , nxp, nyp
    real(ireals),parameter :: dx=100, dy=dx, phi0=0, theta0=60, albedo=.2

    integer(iintegers) :: i,j,k, nlay, ncerr, icld
    character(len=80) :: output_filename(2)

    real(ireals),allocatable :: atm(:,:) ! # z(km)  p(mb)  T(K) air(cm-3) o3(cm-3) o2(cm-3)  h2o(cm-3) co2(cm-3) no2(cm-3)

    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso          ! [nlyr(+1), global_nx, global_ny ]
    integer(iintegers),parameter :: izout=30
    real(ireals),allocatable :: tmp(:,:,:)
    !character(len=300) :: fname(2)

    call mpi_init(mpierr)
    call init_mpi_data_parameters(MPI_COMM_WORLD)

    !call read_ascii_file_2d('afglus_100m.dat', atm, 9, 2)

    !nlay = ubound(atm,1)-1

    if(myid.eq.0) then
        call ncload(['./hill_input.nc', 'plev']  , plev   , ncerr)
        call ncload(['./hill_input.nc', 'tlay']  , tlay   , ncerr)
        call ncload(['./hill_input.nc', 'air']   , air    , ncerr)
        call ncload(['./hill_input.nc', 'h2ovmr'], h2ovmr , ncerr)
        call ncload(['./hill_input.nc', 'o3vmr'] , o3vmr  , ncerr)
        call ncload(['./hill_input.nc', 'co2vmr'], co2vmr , ncerr)
        call ncload(['./hill_input.nc', 'n2ovmr'], n2ovmr , ncerr)
        call ncload(['./hill_input.nc', 'o2vmr'] , o2vmr  , ncerr)
        h2ovmr = h2ovmr / air
        o3vmr  = o3vmr  / air
        co2vmr = co2vmr / air
        n2ovmr = n2ovmr / air
        o2vmr  = o2vmr  / air
        deallocate(air)

        print *,'loaded',myid
        if(myid.eq.0) print *,'plev shape',shape(plev)
        if(myid.eq.0) print *,'h2ovmr shape',shape(h2ovmr)
        if(myid.eq.0) print *,'o3vmr  shape',shape(o3vmr )
        if(myid.eq.0) print *,'co2vmr shape',shape(co2vmr)
        if(myid.eq.0) print *,'n2ovmr shape',shape(n2ovmr)
        if(myid.eq.0) print *,'o2vmr  shape',shape(o2vmr )
    endif
    call imp_bcast(plev  , 0_mpiint, myid)
    call imp_bcast(tlay  , 0_mpiint, myid)
    call imp_bcast(h2ovmr, 0_mpiint, myid)
    call imp_bcast(o3vmr , 0_mpiint, myid)
    call imp_bcast(co2vmr, 0_mpiint, myid)
    call imp_bcast(n2ovmr, 0_mpiint, myid)
    call imp_bcast(o2vmr , 0_mpiint, myid)

    nlay= ubound(plev,1)-1
    nxp = ubound(plev,2)
    nyp = ubound(plev,3)


    allocate(ch4vmr   (nlay  ,nxp,nyp))
    ch4vmr = co2vmr/1e2


    allocate(lwc   (nlay  ,nxp,nyp))
    allocate(reliq (nlay  ,nxp,nyp))
    lwc = 0

    icld = minloc(abs(atm(:,1)-5.5),dim=1)
    lwc(icld:icld+5, nxp/2-5:nxp/2+5, nyp/2-5:nyp/2+5) = .1
    reliq = 10

    if(myid.eq.0) then
        do k=1,nlay
            print *,'plev',plev(k,1,1), 'T', tlay(k,1,1),'CO2',co2vmr(k,1,1),'H2O',h2ovmr(k,1,1)
        enddo
        print *,'plev',plev(nlay+1,1,1)
    endif

    call tenstream_rrtm_sw(imp_comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir, edn, eup, abso)

    if(myid.eq.0) then
        do k=1,nlay+1
            print *,k,'edir', edir(k,1,1), edn(k,1,1), eup(k,1,1)
        enddo

    icld = minloc(abs(atm(:,1)-10),dim=1)
    call fill_nzout(edir); call ncwrite(['output.nc', 'edir'], tmp, ncerr)
    call fill_nzout(edn ); call ncwrite(['output.nc', 'edn'] , tmp, ncerr)
    call fill_nzout(eup ); call ncwrite(['output.nc', 'eup'] , tmp, ncerr)
    call fill_nzout(abso); call ncwrite(['output.nc', 'abso'], tmp, ncerr)
    call fill_nzout(plev); call ncwrite(['output.nc', 'plev'], tmp, ncerr)
    call ncwrite(['output.nc', 'Qnet'], edir(nlay+1,:,:)+edn(nlay+1,:,:), ncerr)
    call ncwrite(['output.nc', 'psrfc'], plev(nlay+1,:,:), ncerr)
    print *,'done',shape(edir)
    endif
    call mpi_barrier(MPI_COMM_WORLD, mpierr)

    contains 
        subroutine fill_nzout(arr)
            real(ireals), intent(in) :: arr(:,:,:)
            if(.not.allocated(tmp)) allocate(tmp(nxp,nyp,izout))
            do k=1,izout
                tmp(:,:,k) = arr(ubound(arr,1)-k+1,:,:)
            enddo
        end subroutine

end program
