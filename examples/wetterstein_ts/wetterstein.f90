module m_wetterstein
  use mpi

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, mpierr

  use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
    tenstream_get_result, getvecpointer, restorevecpointer, &
    t_coord

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast, CHKERR

  use m_netcdfIO, only : ncwrite, ncload

  use m_tenstr_rrtm_sw, only : tenstream_rrtm_sw

  implicit none

contains
  subroutine test_rrtm_sw()

    implicit none

    integer(mpiint) :: comm, myid

    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=180, theta0=60
    real(ireals),parameter :: albedo=0.2
    integer(iintegers),parameter :: icollapse=40

    real(ireals),allocatable,dimension(:,:,:) :: plev                                               ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq, air                                    ! nlay  , nxp, nyp
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso                              ! nlyr(+1), global_nx, global_ny

    character(len=80) :: nc_path(2) ! [ filename, varname ]
    real(ireals),allocatable :: tmp(:,:,:)

    integer(iintegers) :: k
    integer(iintegers) :: nxp,nyp,nlay
    integer(iintegers) :: ncerr

    logical,parameter :: ldebug=.True.

    integer(mpiint) :: ierr

    comm = MPI_COMM_WORLD

    call mpi_comm_rank(comm,myid,ierr)

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
    call tenstream_rrtm_sw(comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir, edn, eup, abso, icollapse=icollapse)

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
        nc_path(2)='Qdir' ; call ncwrite(nc_path, edir(nlay+1,:,:), ncerr)
        nc_path(2)='psrfc'; call ncwrite(nc_path, plev(nlay+1,:,:), ncerr)
        print *,'done',shape(edir)
      endif
    endif
    call mpi_barrier(comm, mpierr)

  contains 
    subroutine fill_nzout(arr)
      real(ireals), intent(in) :: arr(:,:,:)
      integer(iintegers) :: izout=40
      if(.not.allocated(tmp)) allocate(tmp(nxp,nyp,izout))
      do k=1,izout
        tmp(:,:,k) = arr(ubound(arr,1)-k+1,:,:)
      enddo
    end subroutine

  end subroutine
end module

program main
  use mpi
  use m_wetterstein

  integer ierr
  call mpi_init(ierr)

  call test_rrtm_sw()

  call mpi_finalize(ierr)
end program

