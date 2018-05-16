module m_wetterstein
  use mpi

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, default_str_len

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast, CHKERR

  use m_netcdfIO, only : ncwrite, ncload

  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg
  use m_pprts, only : t_solver_3_10

  implicit none

contains
  subroutine ex_wetterstein()

    implicit none

    type(t_solver_3_10) :: solver
    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y

    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=180, theta0=60
    real(ireals),parameter :: albedo_sol=0.2, albedo_th=0
    integer(iintegers),parameter :: icollapse=40

    real(ireals),allocatable,dimension(:,:,:) :: plev, tlev                                         ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay                                               ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq, air                                    ! nlay  , nxp, nyp
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso                              ! nlyr(+1), global_nx, global_ny

    character(len=default_str_len), parameter :: atm_filename='afglus_100m.dat'
    character(len=80) :: nc_path(2) ! [ filename, varname ]
    real(ireals),allocatable :: tmp(:,:,:)

    integer(iintegers) :: k
    integer(iintegers) :: nxp,nyp,nlay
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)
    integer(mpiint) :: ncerr

    logical, parameter :: lthermal=.False., lsolar=.True.

    logical, parameter :: ldebug=.True.

    integer(mpiint) :: ierr

    comm = MPI_COMM_WORLD
    call mpi_comm_size(comm, numnodes, ierr)
    call mpi_comm_rank(comm, myid, ierr)

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate(nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    call init_mpi_data_parameters(comm)

    if(myid.eq.0) then
      nc_path(1) = 'input.nc'
      nc_path(2)='plev'  ;call ncload(nc_path, plev   , ncerr); call CHKERR(ncerr)
      nc_path(2)='tlay'  ;call ncload(nc_path, tlay   , ncerr); call CHKERR(ncerr)
      nc_path(2)='air'   ;call ncload(nc_path, air    , ncerr); call CHKERR(ncerr)
      deallocate(air)

      if(myid.eq.0) print *,'plev shape',shape(plev)
    endif
    call imp_bcast(comm, plev  , 0_mpiint)
    call imp_bcast(comm, tlay  , 0_mpiint)

    nlay= ubound(plev,1)-1
    nxp = ubound(plev,2)
    nyp = ubound(plev,3)

    call CHKERR(1_mpiint, 'TODO: need to define tlev and overall, this example is not tested at all')

    allocate(lwc(nlay, nxp, nyp))
    lwc = 0

    allocate(reliq(nlay, nxp, nyp))
    reliq = 10

    if(myid.eq.0 .and. ldebug) then
      do k=1,nlay
        print *,'plev',plev(k,1,1), 'T', tlay(k,1,1)
      enddo
      print *,'plev',plev(nlay+1,1,1)
    endif

    call pprts_rrtmg(comm, solver, dx, dy, phi0, theta0, albedo_th, albedo_sol, &
      atm_filename, lthermal, lsolar,                                       &
      edir,edn,eup,abso,                                                    &
      d_plev=plev, d_tlev=tlev, d_tlay=tlay, d_lwc=lwc, d_reliq=reliq,      &
      nxproc=nxproc, nyproc=nyproc)

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
    call mpi_barrier(comm, ierr)

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

  call ex_wetterstein()

  call mpi_finalize(ierr)
end program

