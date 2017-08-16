@test(npes =[16])
subroutine tenstream_hill1_x(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, default_str_len

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast, CHKERR

    use m_netcdfIO, only : ncwrite, ncload

    use m_tenstr_rrtm_sw, only : tenstream_rrtm_sw

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0.2
    real(ireals),parameter :: atolerance = 1

    real(ireals),allocatable,dimension(:,:,:) :: plev                                               ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq, air                                    ! nlay  , nxp, nyp
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso                ! nlyr(+1), global_nx, global_ny

    character(default_str_len) :: nc_path(2) ! [ filename, varname ]
    real(ireals),allocatable :: tmp(:,:,:)

    integer(iintegers) :: k
    integer(iintegers) :: nxp,nyp,nlay
    integer(iintegers) :: ncerr

    logical,parameter :: ldebug=.True.

    real(ireals),allocatable :: target_edir(:)
    real(ireals),allocatable :: target_edn(:)
    real(ireals),allocatable :: target_abso(:)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    if(myid.eq.0) then
        nc_path(1) = 'input.nc'
        nc_path(2)='plev'  ;call ncload(nc_path, plev   , ncerr); call CHKERR(ncerr)
        nc_path(2)='tlay'  ;call ncload(nc_path, tlay   , ncerr); call CHKERR(ncerr)
        nc_path(2)='air'   ;call ncload(nc_path, air    , ncerr); call CHKERR(ncerr)
        nc_path(2)='h2ovmr';call ncload(nc_path, h2ovmr , ncerr); call CHKERR(ncerr)
        nc_path(2)='o3vmr' ;call ncload(nc_path, o3vmr  , ncerr); call CHKERR(ncerr)
        nc_path(2)='co2vmr';call ncload(nc_path, co2vmr , ncerr); call CHKERR(ncerr)
        nc_path(2)='n2ovmr';call ncload(nc_path, n2ovmr , ncerr); call CHKERR(ncerr)
        nc_path(2)='o2vmr' ;call ncload(nc_path, o2vmr  , ncerr); call CHKERR(ncerr)
        h2ovmr = h2ovmr / air
        o3vmr  = o3vmr  / air
        co2vmr = co2vmr / air
        n2ovmr = n2ovmr / air
        o2vmr  = o2vmr  / air
        deallocate(air)

        if(myid.eq.0) print *,'plev shape',shape(plev)
        if(myid.eq.0) print *,'h2ovmr shape',shape(h2ovmr)
        if(myid.eq.0) print *,'o3vmr  shape',shape(o3vmr )
        if(myid.eq.0) print *,'co2vmr shape',shape(co2vmr)
        if(myid.eq.0) print *,'n2ovmr shape',shape(n2ovmr)
        if(myid.eq.0) print *,'o2vmr  shape',shape(o2vmr )
    endif
    call imp_bcast(comm, plev  , 0_mpiint)
    call imp_bcast(comm, tlay  , 0_mpiint)
    call imp_bcast(comm, h2ovmr, 0_mpiint)
    call imp_bcast(comm, o3vmr , 0_mpiint)
    call imp_bcast(comm, co2vmr, 0_mpiint)
    call imp_bcast(comm, n2ovmr, 0_mpiint)
    call imp_bcast(comm, o2vmr , 0_mpiint)

    nlay= ubound(plev,1)-1
    nxp = ubound(plev,2)
    nyp = ubound(plev,3)

    allocate(ch4vmr   (nlay  ,nxp,nyp))
    ch4vmr = co2vmr/1e2

    allocate(lwc   (nlay  ,nxp,nyp))
    lwc = 0

    allocate(reliq (nlay  ,nxp,nyp))
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

          nc_path(1) = 'output.nc'
          print *,'writing output to file', nc_path(1)
          nc_path(2)='edir' ; call fill_nzout(edir); call ncwrite(nc_path, tmp, ncerr)
          nc_path(2)='edn'  ; call fill_nzout(edn ); call ncwrite(nc_path, tmp, ncerr)
          nc_path(2)='eup'  ; call fill_nzout(eup ); call ncwrite(nc_path, tmp, ncerr)
          nc_path(2)='abso' ; call fill_nzout(abso); call ncwrite(nc_path, tmp, ncerr)
          nc_path(2)='plev' ; call fill_nzout(plev); call ncwrite(nc_path, tmp, ncerr)
          nc_path(2)='Qnet' ; call ncwrite(nc_path, edir(nlay+1,:,:)+edn(nlay+1,:,:), ncerr)
          nc_path(2)='psrfc'; call ncwrite(nc_path, plev(nlay+1,:,:), ncerr)
          print *,'done',shape(edir)
          print *,'abso _srfc',abso(nlay  ,2,:)
        endif

        allocate( target_edir(nyp) )
        target_edir = [ 479.4954, 469.4137, 461.5663, 453.6571, 446.9970, 442.0156,    &
          438.7523, 437.1892, 437.4877, 439.7622, 443.9711, 450.0184, 457.4526, &
          467.4857, 476.4721, 485.6826, 492.7633, 500.9566, 509.2136, 515.4571, &
          520.4636, 523.8958, 525.6567, 525.7076, 524.1005, 520.9698, 516.4503, &
          510.9528, 505.1550, 496.9903, 488.7550]

        allocate( target_edn(nyp) )
        target_edn = [43.51924, 43.49717, 43.51923, 43.49463, 43.49198, &
          43.50251, 43.52306, 43.54593, 43.57961, 43.62609, 43.67543, 43.72336, &
          43.77139, 43.82191, 43.83408, 43.87674, 43.89641, 43.94455, 43.94963, &
          43.94936, 43.93943, 43.91687, 43.88498, 43.85325, 43.81978, 43.77658, &
          43.73243, 43.69201, 43.66598, 43.58746, 43.54970]

        allocate( target_abso(nyp) )
        target_abso = [1.5701262E-02, 1.5620301E-02, 1.5592689E-02,   &
          1.5966615E-02, 1.6579399E-02, 1.7405437E-02, 1.8453605E-02, &
          1.9726738E-02, 2.1241061E-02, 2.2985294E-02, 2.4917131E-02, &
          2.6942026E-02, 2.8901748E-02, 3.0621707E-02, 3.1606011E-02, &
          3.2174263E-02, 3.1969056E-02, 3.1328052E-02, 2.9930128E-02, &
          2.8154518E-02, 2.6231173E-02, 2.4317959E-02, 2.2520296E-02, &
          2.0905165E-02, 1.9501472E-02, 1.8310279E-02, 1.7332951E-02, &
          1.6557695E-02, 1.6030772E-02, 1.5902195E-02, 1.5793735E-02]

        @assertEqual( target_edir, edir(nlay+1,2,:), atolerance, 'solar at surface :: edir       not correct')
        @assertEqual( target_edn , edn (nlay+1,2,:), atolerance, 'solar at surface :: edn        not correct')
        @assertEqual( target_abso, abso(nlay+1,2,:), atolerance, 'solar at surface :: absorption not correct')
    endif

    contains 
        subroutine fill_nzout(arr)
            real(ireals), intent(in) :: arr(:,:,:)
            integer(iintegers) :: izout = 60
            if(.not.allocated(tmp)) allocate(tmp(nxp,nyp,izout))
            do k=1,izout
                tmp(:,:,k) = arr(ubound(arr,1)-k+1,:,:)
            enddo
        end subroutine

end subroutine
