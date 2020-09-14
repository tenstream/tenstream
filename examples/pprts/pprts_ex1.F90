module m_pprts_ex1
    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi
    use m_helper_functions, only : spherical_2_cartesian
    use m_pprts, only : init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles
    use m_pprts_base, only: t_solver, t_solver_1_2, t_solver_3_6, t_solver_3_10, &
      t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18, destroy_pprts
    use m_tenstream_options, only: read_commandline_options

    use mpi, only : MPI_COMM_WORLD

  implicit none

  contains
subroutine pprts_ex1()
    integer(iintegers),parameter :: nxp=9,nyp=9,nv=40
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=270, theta0=20
    real(ireals),parameter :: albedo=0., dz=100
    real(ireals),parameter :: incSolar = 1364
    real(ireals) :: dz1d(nv)

    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver

    integer(iintegers) :: mid_idx_x, mid_idx_y, mid_idx_z
    character(len=80) :: arg

    dz1d = dz

    call get_command_argument(1, arg)

    select case(arg)
      case('-solver_1_2')
        allocate(t_solver_1_2::solver)
      case('-solver_3_6')
        allocate(t_solver_3_6::solver)
      case ('-solver_3_10')
        allocate(t_solver_3_10::solver)
      case ('-solver_8_10')
        allocate(t_solver_8_10::solver)
      case ('-solver_3_16')
        allocate(t_solver_3_16::solver)
      case ('-solver_8_16')
        allocate(t_solver_8_16::solver)
      case ('-solver_8_18')
        allocate(t_solver_8_18::solver)
      case default
        print *,'error, have to provide solver type as argument, e.g.'
        print *,'-solver_1_2'
        print *,'-solver_3_6'
        print *,'-solver_3_10'
        print *,'-solver_8_10'
        print *,'-solver_3_16'
        print *,'-solver_8_16'
        print *,'-solver_8_18'
        stop
    end select

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(mpi_comm_world, nv, nxp, nyp, dx,dy, sundir, solver, dz1d)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = .1_ireals/(dz*nv)
    ksca = 1e-3_ireals/(dz*nv)
    g    = zero

    mid_idx_x = int(real(nxp+1)/2)
    mid_idx_y = int(real(nyp+1)/2)
    mid_idx_z = int(real(nv +1)/2)

    kabs(mid_idx_z, mid_idx_x, :) = 1/dz
    ksca(mid_idx_z, mid_idx_x, :) = 1/dz
    g   (mid_idx_z, mid_idx_x, :) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call set_angles(solver, sundir)

    call solve_pprts(solver, &
      & lthermal=.False., &
      & lsolar=.True., &
      & edirTOA=incSolar )

    call pprts_get_result(solver, fdn,fup,fdiv,fdir)

    print *,'edir', fdir(:, mid_idx_x, mid_idx_y)
    print *,'edn:', fdn (:, mid_idx_x, mid_idx_y)
    print *,'eup:', fup (:, mid_idx_x, mid_idx_y)
    print *,'divE', fdiv(:, mid_idx_x, mid_idx_y)
    call destroy_pprts(solver, .True.)
end subroutine

end module


program main
  use m_pprts_ex1
  integer(mpiint) :: myid, ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

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
  call mpi_finalize(ierr)
end program
