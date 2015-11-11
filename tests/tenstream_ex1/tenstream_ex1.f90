!todo teardown of some stuff does not work correctly... need fix to use automatic parallel testing
@test(npes =[8]) 
subroutine test_tenstream_ex1(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,mpierr,zero,pi

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord,C_dir,C_diff,C_one,C_one1

    use m_tenstream_options, only: read_commandline_options

#include "petsc/finclude/petscdef.h"
    use petsc 

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: k, numnodes, myid, comm

    integer(iintegers),parameter :: nxp=20,nyp=20,nv=10
    real(ireals),parameter :: dx=67,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0, dz=dx
    real(ireals),parameter :: incSolar = -1
    real(ireals),parameter :: tolerance = 1e-3
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,B
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: div_target(nv)
    real(ireals) :: dn_target(nv+1)
    real(ireals) :: up_target(nv+1)

    PetscErrorCode :: ierr

    div_target = [-2.7829354792270369, -8.2537454026863771E-003, 1.3956358754361362E-002, 1.4644200952248162E-002, 1.5171323114716945E-002, 1.5706565902073349E-002, 1.6252617408845316E-002, 1.6796285431804895E-002, 1.5835237341535743E-002, -0.17650992850077310]
    dn_target  = [0.0000000000000000,  202.03307328725239,  219.03714311366440,  235.46022830171623,  252.77735226449187,  271.03260893291355,  290.25932657631955,  310.49131813452840,  331.76306602630132,  354.10961761810569,  377.56656927800549]
    up_target  = [204.07584410995435,  219.65261706673965,  236.10460087064266,  253.46395412430059,  271.76355087489213,  291.03669390961892,  311.31725943277547,  332.63978969874182,  355.03860747633070,  378.44788368251432,  390.07937341058152]
    
    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

!    print *,'tenstream_test',myid,'/',numnodes


    PETSC_COMM_WORLD = comm
    call init_tenstream(comm, nv, nxp,nyp, dx,dy,phi0, theta0, albedo, dz1d=dz1d)

    call init_mpi_data_parameters(comm)
    call read_commandline_options()
!    print *,myid,'my local domain size:',C_one%xm, C_one%ym, C_one%zm

    allocate(kabs(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(ksca(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(g   (C_one%zm , C_one%xm,  C_one%ym ))
    allocate(B   (C_one1%zm, C_one1%xm, C_one1%ym))

    kabs = 1._ireals/nv
    ksca = 1e-8
    g    = zero
    do k=1,C_one1%zm
        B(k,:,:) = (288 - 50*(1 - float(k)/float(nv+1)) )**4 * 5.67e-8/ pi
!        print *,'T',k,288 - 50*(1 - float(k)/float(nv+1))
    enddo

    call set_optical_properties( kabs, ksca, g, B )
    call solve_tenstream(incSolar)

    allocate(fdn  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fup  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fdiv (C_one%zm, C_one%xm, C_one%ym))

    call tenstream_get_result(fdir, fdn, fup, fdiv)

!    if(myid.eq.0) then
!        print *,'B:', B(:,1,1)
!        print *,'kabs:', kabs(:,1,1)
!        print *,'fdn:',  fdn (:,1,1)
!        print *,'fup:',  fup (:,1,1)
!        print *,'fdiv:', fdiv(:,1,1)
!    endif

    call destroy_tenstream()

    @assertEqual(div_target, fdiv(:,1,1), tolerance, 'thermal divergence not correct')
    @assertEqual(dn_target,  fdn (:,1,1), tolerance, 'thermal downw flux not correct')
    @assertEqual(up_target,  fup (:,1,1), tolerance, 'thermal upward fl  not correct')

    ! Check that surface emission is the one that we stick in
    @assertEqual(B(C_one1%zm,:,:)*pi, fup (C_one1%zm,:,:), tolerance, 'Surface Emission not correct')
end subroutine
