module m_box_cld

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi

  use m_pprts, only : init_pprts, set_optical_properties, solve_pprts, destroy_pprts, &
    pprts_get_result
  use m_pprts_base, only: t_solver_3_10

  use m_helper_functions, only : get_mem_footprint

  use m_tenstream_options, only: read_commandline_options

  use mpi, only : MPI_COMM_WORLD

  implicit none

  integer(mpiint) :: myid, ierr

  contains
subroutine box_cld()
    implicit none

    integer(iintegers),parameter :: nxp=3,nyp=3,nv=20
    real(ireals),parameter :: dx=500,dy=dx
    real(ireals),parameter :: phi0=90, theta0=0
    real(ireals),parameter :: albedo=0, dz=10
    real(ireals),parameter :: incSolar = 1364
    real(ireals),parameter :: atolerance = 1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    type(t_solver_3_10) :: solver

    dz1d = dz

    call init_pprts(MPI_COMM_WORLD, nv, nxp,nyp, dx,dy,phi0, theta0, solver, dz1d=dz1d)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = .0001_ireals/(dz*nv)
    ksca = .005_ireals/(dz*nv)
    g    = zero

    !kabs(nv/2,nxp/2,1:nyp) = 1/dz
    !ksca(nv/2,nxp/2,1:nyp) = 1/dz
    !g   (nv/2,nxp/2,1:nyp) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call solve_pprts(solver, incSolar)

    allocate(fdir (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    call pprts_get_result(solver, fdn, fup, fdiv, fdir)

    if(myid.eq.0) then
        print *,'kabs:', kabs(:,1,1)
        print *,'fdir:', fdir(:,1,1)
        print *,'fdn:',  fdn (:,1,1)
        print *,'fup:',  fup (:,1,1)
        print *,'fdiv:', fdiv(:,1,1)
    endif
    print *,myid,'Memory:',get_mem_footprint(MPI_COMM_WORLD)

    call destroy_pprts(solver, .True.)
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
