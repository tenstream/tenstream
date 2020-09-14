module m_box_cld

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi

  use m_pprts, only : init_pprts, set_optical_properties, solve_pprts, &
    pprts_get_result
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts

  use m_helper_functions, only : CHKERR, get_mem_footprint, spherical_2_cartesian

  use m_tenstream_options, only: read_commandline_options

  use mpi, only : MPI_COMM_WORLD

  implicit none

  integer(mpiint) :: myid, ierr

  contains
subroutine box_cld()
    implicit none

    integer(iintegers),parameter :: nxp=16,nyp=16,nv=20
    real(ireals),parameter :: dx=500,dy=dx
    real(ireals),parameter :: phi0=180, theta0=40
    real(ireals),parameter :: albedo=.1, dz=100
    real(ireals),parameter :: incSolar = 1364
    integer(iintegers), parameter :: cld_width=0
    real(ireals) :: dz1d(nv)

    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(MPI_COMM_WORLD)

    call allocate_pprts_solver_from_commandline(solver, '8_16', ierr); call CHKERR(ierr)
    dz1d = dz

    sundir = spherical_2_cartesian(phi0, theta0)
    print *,'sundir', sundir

    call init_pprts(MPI_COMM_WORLD, nv, nxp, nyp, dx, dy, sundir, solver, dz1d=dz1d)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = .0001_ireals/(dz*nv)
    ksca = .005_ireals/(dz*nv)
    g    = zero

    kabs(nv/2, nxp/2-cld_width:nxp/2+cld_width, nyp/2-cld_width:nyp/2+cld_width) = 1/dz
    ksca(nv/2, nxp/2-cld_width:nxp/2+cld_width, nyp/2-cld_width:nyp/2+cld_width) = 1/dz
    g   (nv/2, nxp/2-cld_width:nxp/2+cld_width, nyp/2-cld_width:nyp/2+cld_width) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call solve_pprts(solver, lthermal=.False., lsolar=.True., edirTOA=incSolar)

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
