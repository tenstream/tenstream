program main

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,imp_comm,myid,mpierr,zero, i0
    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec

    use m_tenstr_rrtm_sw
    implicit none

    integer(iintegers),parameter :: nxp=16, nyp=nxp
    real(ireals),allocatable,dimension(:,:,:) :: plev                                                     ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr       ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq                                               ! nlay  , nxp, nyp

    real(ireals),parameter :: dx=200, dy=dx, phi0=0, theta0=0, albedo=.0

    integer(iintegers) :: i,j,k, nlay

    real(ireals),allocatable :: atm(:,:) ! # z(km)  p(mb)  T(K) air(cm-3) o3(cm-3) o2(cm-3)  h2o(cm-3) co2(cm-3) no2(cm-3)

    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso          ! [nlyr(+1), local_nx, local_ny ]

    call mpi_init(mpierr)
    call init_mpi_data_parameters(MPI_COMM_WORLD)

    call read_ascii_file_2d('afglus.dat', atm, 9, 2)

    nlay = ubound(atm,1)-1

    allocate(plev  (nlay+1,nxp,nyp))
    allocate(tlay  (nlay  ,nxp,nyp))
    allocate(h2ovmr(nlay  ,nxp,nyp))
    allocate(o3vmr (nlay  ,nxp,nyp))
    allocate(co2vmr(nlay  ,nxp,nyp))
    allocate(ch4vmr(nlay  ,nxp,nyp))
    allocate(n2ovmr(nlay  ,nxp,nyp))
    allocate(o2vmr (nlay  ,nxp,nyp))

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


    allocate(lwc   (nlay  ,nxp,nyp))
    allocate(reliq (nlay  ,nxp,nyp))
    lwc = 0
    lwc(40, 1, 1) = 1e-2
    reliq = 10

    if(myid.eq.0) then
        do k=1,nlay
            print *,'plev',plev(k,1,1), 'T', tlay(k,1,1),'CO2',co2vmr(k,1,1),'H2O',h2ovmr(k,1,1)
        enddo
        print *,'plev',plev(nlay+1,1,1)
    endif

    call tenstr_rrtm_sw(imp_comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir, edn, eup, abso)

    if(myid.eq.0) then
        do k=1,nlay+1
            print *,k,'edir', edir(k,1,1), edn(k,1,1), eup(k,1,1)
        enddo
    endif
    print *,'done'

end program
