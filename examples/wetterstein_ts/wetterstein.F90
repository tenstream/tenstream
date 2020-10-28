module m_wetterstein
  use mpi

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, default_str_len

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only : &
    & CHKERR, &
    & domain_decompose_2d_petsc, &
    & imp_bcast, &
    & spherical_2_cartesian

  use m_netcdfIO, only : ncwrite, ncload

  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm

  implicit none

contains
  subroutine ex_wetterstein(comm)
    integer(mpiint), intent(in) :: comm

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, myid

    integer(iintegers) :: Nx, Ny, Nlay   ! global domain size
    integer(iintegers) :: nxp, nyp, xs, ys ! local domain size in x and y aswell as start indices
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    class(t_solver), allocatable :: solver
    type(t_tenstr_atm)  :: atm

    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=180, theta0=60
    real(ireals),parameter :: albedo_sol=0.2, albedo_th=0

    real(ireals),allocatable,dimension(:,:,:), target :: plev, tlev    ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:), target :: lwc, reliq    ! nlay  , nxp, nyp
    real(ireals),pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq ! reshape pointers to convert to column vecs
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! nlyr(+1), global_nx, global_ny

    character(len=default_str_len), parameter :: atm_filename='afglus_100m.dat'
    character(len=80) :: nc_path(2) ! [ filename, varname ]
    real(ireals),allocatable :: tmp(:,:,:)

    real(ireals) :: sundir(3)
    integer(iintegers) :: k
    integer(mpiint) :: ncerr, ierr

    logical, parameter :: lthermal=.False., lsolar=.True.
    logical, parameter :: ldebug=.True.

    call init_mpi_data_parameters(comm)
    call mpi_comm_size(comm, numnodes, ierr)
    call mpi_comm_rank(comm, myid, ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      nc_path(1) = 'input.nc'
      nc_path(2)='plev'  ;call ncload(nc_path, plev   , ncerr); call CHKERR(ncerr)
      nc_path(2)='tlev'  ;call ncload(nc_path, tlev   , ncerr); call CHKERR(ncerr)

      if(myid.eq.0) print *,'plev shape',shape(plev)
    endif
    call imp_bcast(comm, plev  , 0_mpiint)
    call imp_bcast(comm, tlev  , 0_mpiint)

    Nlay = ubound(plev,1)-1
    Nx   = ubound(plev,2)
    Ny   = ubound(plev,3)

    call domain_decompose_2d_petsc(comm, Nx, Ny, &
      & nxp, nyp, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)

    if(myid.eq.0) print *, myid, 'Domain Decomposition on', numnodes, 'ranks will be', nxproc, 'and', nyproc

    call CHKERR(1_mpiint, 'TODO: need to define tlev and overall, this example is not tested at all')

    allocate(lwc(nlay, nxp, nyp))
    lwc = 0

    allocate(reliq(nlay, nxp, nyp))
    reliq = 10

    if(myid.eq.0 .and. ldebug) then
      do k=1,nlay+1
        print *,'plev',plev(k,1,1), 'T', tlev(k,1,1)
      enddo
    endif

    sundir = spherical_2_cartesian(phi0, theta0)

    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm, &
      d_lwc=plwc, d_reliq=preliq)

    call pprts_rrtmg(comm, solver, atm, nxp, nyp, &
      dx, dy, sundir,         &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir, edn, eup, abso,   &
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
  use mpi, only: mpi_init, mpi_finalize, MPI_COMM_WORLD
  use m_wetterstein

  integer ierr
  call mpi_init(ierr)

  call ex_wetterstein(MPI_COMM_WORLD)

  call mpi_finalize(ierr)
end program

