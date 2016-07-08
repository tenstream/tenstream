@test(npes =[8])
subroutine tenstream_hill1_x(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,mpierr,zero,pi

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord,C_dir,C_diff,C_one,C_one1

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast

    use m_netcdfIO, only : ncwrite, ncload

    use m_tenstr_rrtm_sw, only : tenstream_rrtm_sw

#include "petsc/finclude/petscdef.h"
    use petsc 

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=0
    real(ireals),parameter :: albedo=0.2, dz=dx
    real(ireals),parameter :: atolerance = 1
    real(ireals),parameter :: rtolerance = .05

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,B
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals),allocatable,dimension(:,:,:) :: plev                                               ! nlay+1, nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr ! nlay  , nxp, nyp
    real(ireals),allocatable,dimension(:,:,:) :: lwc, reliq, air                                    ! nlay  , nxp, nyp
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! nlyr(+1), global_nx, global_ny

    real(ireals),allocatable :: tmp(:,:,:)

    integer(iintegers) :: i,j,k 
    integer(iintegers) :: nxp,nyp,nlay
    integer(iintegers) :: ncerr

    logical,parameter :: ldebug=.True.

    PetscErrorCode :: ierr

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    if(myid.eq.0) then
        call ncload(['./hill1_x_input.nc', 'plev']  , plev   , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'tlay']  , tlay   , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'air']   , air    , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'h2ovmr'], h2ovmr , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'o3vmr'] , o3vmr  , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'co2vmr'], co2vmr , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'n2ovmr'], n2ovmr , ncerr); CHKERRQ(ncerr)
        call ncload(['./hill1_x_input.nc', 'o2vmr'] , o2vmr  , ncerr); CHKERRQ(ncerr)
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
    call imp_bcast(comm, plev  , 0_mpiint, myid)
    call imp_bcast(comm, tlay  , 0_mpiint, myid)
    call imp_bcast(comm, h2ovmr, 0_mpiint, myid)
    call imp_bcast(comm, o3vmr , 0_mpiint, myid)
    call imp_bcast(comm, co2vmr, 0_mpiint, myid)
    call imp_bcast(comm, n2ovmr, 0_mpiint, myid)
    call imp_bcast(comm, o2vmr , 0_mpiint, myid)

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

    if(myid.eq.0.and.ldebug) then
        if(ldebug) then
            do k=1,nlay+1
                print *,k,'edir', edir(k,1,1), edn(k,1,1), eup(k,1,1), abso(min(nlay,k),1,1)
            enddo
        endif

        print *,'writing output to file'
        call fill_nzout(edir); call ncwrite(['output.nc', 'edir'], tmp, ncerr)
        call fill_nzout(edn ); call ncwrite(['output.nc', 'edn'] , tmp, ncerr)
        call fill_nzout(eup ); call ncwrite(['output.nc', 'eup'] , tmp, ncerr)
        call fill_nzout(abso); call ncwrite(['output.nc', 'abso'], tmp, ncerr)
        call fill_nzout(plev); call ncwrite(['output.nc', 'plev'], tmp, ncerr)
        call ncwrite(['output.nc', 'Qnet'], edir(nlay+1,:,:)+edn(nlay+1,:,:), ncerr)
        call ncwrite(['output.nc', 'psrfc'], plev(nlay+1,:,:), ncerr)
        print *,'done',shape(edir)

        @assertEqual(0.0, edir(nlay+1,1,1), atolerance, 'solar at surface :: edir       not correct')
        @assertEqual(0.0, edn (nlay+1,1,1), atolerance, 'solar at surface :: downw flux not correct')
        @assertEqual(0.0, eup (nlay+1,1,1), atolerance, 'solar at surface :: upward fl  not correct')
        @assertEqual(0.0, abso(nlay  ,1,1), atolerance, 'solar at surface :: absorption not correct')

        @assertEqual(0.0, edir(1,1,1), atolerance, 'solar at TOA :: edir       not correct')
        @assertEqual(0.0, edn (1,1,1), atolerance, 'solar at TOA :: downw flux not correct')
        @assertEqual(0.0, eup (1,1,1), atolerance, 'solar at TOA :: upward fl  not correct')
        @assertEqual(0.0, abso(1,1,1), atolerance, 'solar at TOA :: absorption not correct')

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
