@test(npes =[1]) 
subroutine test_rrtm_lw(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord

    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast

    use m_tenstr_rrtm_lw, only : tenstream_rrtm_lw

#include "petsc/finclude/petscdef.h"
    use petsc 

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: nxp=3, nyp=3, nzp=10 ! local domain size for each rank
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: albedo=0, phi0=180, theta0=0
    real(ireals),parameter :: atolerance = 1

    real(ireals), dimension(nzp+1,nxp,nyp) :: plev, tlev
    real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr 
    real(ireals), dimension(nzp,nxp,nyp) :: lwc, reliq                                         
    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso ! [nlev_merged(-1), nxp, nyp]

    character(len=250),parameter :: atm_filename='afglus_100m.dat'

    integer(iintegers) :: i,j,k, nlev, icld

    logical,parameter :: ldebug=.True.

    PetscErrorCode :: ierr

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

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


!    call tenstream_rrtm_lw(comm, dx, dy, phi0, theta0, albedo, atm_filename, &
!                               edir,edn,eup,abso,                            & 
!                               plev, tlev, tlay, h2ovmr, o3vmr,              &
!                               co2vmr, ch4vmr, n2ovmr,  o2vmr,               &
!                               lwc, reliq)

    call tenstream_rrtm_lw(comm, dx, dy, phi0, theta0, albedo, atm_filename, &
                               edir,edn,eup,abso,                            & 
                               plev, tlev, d_lwc=lwc, d_reliq=reliq)

    nlev = ubound(edn,1)
    if(myid.eq.0) then
        if(ldebug) then
            do k=1,nlev
                print *,k,'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
            enddo
        endif

        @assertEqual(300.4933, edn (nlev,1,1), atolerance, 'thermal at surface :: downw flux not correct')
        @assertEqual(391.1564, eup (nlev,1,1), atolerance, 'thermal at surface :: upward fl  not correct')
        @assertEqual(-3.1456E-02, abso(nlev  ,1,1), atolerance, 'thermal at surface :: absorption not correct')

        @assertEqual(  0.0000, edn (1,1,1), atolerance, 'thermal at TOA :: downw flux not correct')
        @assertEqual(244.0147, eup (1,1,1), atolerance, 'thermal at TOA :: upward fl  not correct')
        @assertEqual(-1.4760E-04, abso(1,1,1), atolerance, 'thermal at TOA :: absorption not correct')

        @assertEqual(163.9768, edn (icld+1,1,1), atolerance, 'thermal at icloud :: downw flux not correct')
        @assertEqual(270.4422, eup (icld  ,1,1), atolerance, 'thermal at icloud :: upward fl  not correct')
        @assertEqual(-0.69361, abso(icld  ,1,1), atolerance, 'thermal at icloud :: absorption not correct')
    endif

end subroutine
