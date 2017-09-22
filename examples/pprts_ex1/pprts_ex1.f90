module m_pprts_ex1

  implicit none


  contains
subroutine pprts_ex1()
    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, zero, pi, myid, numnodes

    use m_pprts, only : init_pprts, t_solver
    use m_tenstream_options, only: read_commandline_options

    use mpi, only : MPI_COMM_WORLD

    implicit none

    integer(iintegers) :: k

    integer(iintegers),parameter :: nxp=30,nyp=10,nv=10
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=90, theta0=60
    real(ireals),parameter :: albedo=0, dz=100
    real(ireals) :: dz1d(nv)
    integer(iintegers) :: dofdiff(3), dofdir(3)
    
   
    type(t_solver), allocatable :: s1, s2


    dz1d = dz
    dofdiff = [2,2,2]
    dofdir = [1,1,1]

    call init_pprts(MPI_Comm_World, nv, nxp, nyp, dx,dy, phi0, theta0, dofdiff, dofdir, s1, dz1d)

    dofdiff = [2,4,4]
    dofdir = [4,2,2]

    call init_pprts(MPI_Comm_World, nv, nxp, nyp, dx,dy, phi0, theta0, dofdiff, dofdir, s2, dz1d)


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
