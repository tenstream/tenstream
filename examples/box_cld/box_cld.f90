module m_box_cld

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi

  use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
    tenstream_get_result, getvecpointer, restorevecpointer, &
    t_coord,C_diff,C_one,C_one1, get_mem_footprint

  use m_tenstream_options, only: read_commandline_options

  use mpi, only : MPI_COMM_WORLD

  implicit none

  integer(mpiint) :: myid, ierr

  contains
subroutine box_cld()
    implicit none

    integer(iintegers) :: k

    integer(iintegers),parameter :: nxp=30,nyp=10,nv=10
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=90, theta0=60
    real(ireals),parameter :: albedo=0, dz=100
    real(ireals),parameter :: incSolar = 1364
    real(ireals),parameter :: atolerance = 1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: div_target(nv)
    real(ireals) :: dn_target(nv+1)
    real(ireals) :: up_target(nv+1)

    dz1d = dz

    call init_tenstream(MPI_COMM_WORLD, nv, nxp,nyp, dx,dy,phi0, theta0, dz1d=dz1d)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

    allocate(kabs(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(ksca(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(g   (C_one%zm , C_one%xm,  C_one%ym ))

    kabs = .1_ireals/(dz*nv)
    ksca = zero !1e-3_ireals/dz
    g    = zero

    kabs(nv/2,nxp/2,1:nyp) = 1/dz
    ksca(nv/2,nxp/2,1:nyp) = 1/dz
    g   (nv/2,nxp/2,1:nyp) = .9

    call set_optical_properties(albedo, kabs, ksca, g)
    call solve_tenstream(incSolar)

    allocate(fdir (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fdn  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fup  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fdiv (C_one%zm, C_one%xm, C_one%ym))

    call tenstream_get_result(fdir, fdn, fup, fdiv)

    if(myid.eq.0) then
        print *,'kabs:', kabs(:,1,1)
        print *,'fdir:', fdir(:,1,1)
        print *,'fdn:',  fdn (:,1,1)
        print *,'fup:',  fup (:,1,1)
        print *,'fdiv:', fdiv(:,1,1)
    endif
    print *,myid,'Memory:',get_mem_footprint()

    call destroy_tenstream(.True.)
end subroutine
end module

program main
  use m_box_cld

  call box_cld()

  if(myid.eq.0) then
    print *,''
    print *,''
    print *,'Call this example e.g. with options: -show_edir hdf5:edir.h5'
    print *,'and plot results with python:'
    print *,'import h5py as H; h=H.File("edir.h5","r"); edir = h["edir0"][:]'
    print *,'imshow(edir[0,:,:,0].T,interpolation="nearest");' ! has dimension nyp,nxp,nzp,8streams
    print *,'colorbar(); savefig("edir_x0.pdf")'
  endif

end program
