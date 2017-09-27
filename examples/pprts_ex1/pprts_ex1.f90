module m_pprts_ex1

  implicit none


  contains
subroutine pprts_ex1()
    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, zero, pi, myid, numnodes

    use m_pprts, only : init_pprts, set_optical_properties, t_solver_3_6, t_solver_8_10
    use m_tenstream_options, only: read_commandline_options

    use mpi, only : MPI_COMM_WORLD

    implicit none

    integer(iintegers) :: k

    integer(iintegers),parameter :: nxp=30,nyp=10,nv=10
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=90, theta0=60
    real(ireals),parameter :: albedo=0, dz=100
    real(ireals) :: dz1d(nv)
    
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv
   
    type(t_solver_3_6)  :: s1
    type(t_solver_8_10) :: s2


    dz1d = dz

    call init_pprts(MPI_Comm_World, nv, nxp, nyp, dx,dy, phi0, theta0, s1, dz1d)

    call init_pprts(MPI_Comm_World, nv, nxp, nyp, dx,dy, phi0, theta0, s2, dz1d)

    allocate(kabs(s1%C_one%zm , s1%C_one%xm,  s1%C_one%ym ))
    allocate(ksca(s1%C_one%zm , s1%C_one%xm,  s1%C_one%ym ))
    allocate(g   (s1%C_one%zm , s1%C_one%xm,  s1%C_one%ym ))

    kabs = .1_ireals/(dz*nv)
    ksca = zero !1e-3_ireals/dz
    g    = zero

    kabs(nv/2,nxp/2,1:nyp) = 1/dz
    ksca(nv/2,nxp/2,1:nyp) = 1/dz
    g   (nv/2,nxp/2,1:nyp) = .9

    call set_optical_properties(s1, albedo, kabs, ksca, g)
    call set_optical_properties(s2, albedo, kabs, ksca, g)

end subroutine
end module

program main
  use m_pprts_ex1

  call pprts_ex1()

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
