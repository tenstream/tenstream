module m_py_boxmc
  use m_data_parameters, only : ireals, mpiint, iintegers, &
                                init_mpi_data_parameters
  use m_boxmc, only : t_boxmc_8_10
  use mpi

  implicit none

  private
  public :: get_coeff_8_10

  contains

    subroutine get_coeff_8_10( &
        comm, op_bg, src,      &
        ldir, phi0, theta0,    &
        dx, dy, dz,            &
        S, T, S_tol,T_tol,     &
        inp_atol, inp_rtol)

    real(ireals),intent(in)       :: op_bg(3)     !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
    real(ireals),intent(in)       :: phi0         !< @param[in] phi0 solar azimuth angle
    real(ireals),intent(in)       :: theta0       !< @param[in] theta0 solar zenith angle
    integer(iintegers),intent(in) :: src          !< @param[in] src stream from which to start photons - see init_photon routines
    integer(mpiint),intent(in)    :: comm         !< @param[in] comm MPI Communicator
    logical,intent(in)            :: ldir         !< @param[in] ldir determines if photons should be started with a fixed incidence angle
    real(ireals),intent(in)       :: dx,dy,dz     !< @param[in] dx,dy,dz box with dimensions in [m]
    real(ireals),intent(out)      :: S(10)        !< @param[out] S_out diffuse streams transfer coefficients
    real(ireals),intent(out)      :: T(8)         !< @param[out] T_out direct streams transfer coefficients
    real(ireals),intent(out)      :: S_tol(10)    !< @param[out] absolute tolerances of results
    real(ireals),intent(out)      :: T_tol(8)     !< @param[out] absolute tolerances of results
    real(ireals),intent(in),optional :: inp_atol  !< @param[in] inp_atol if given, determines targeted absolute stddeviation
    real(ireals),intent(in),optional :: inp_rtol  !< @param[in] inp_rtol if given, determines targeted relative stddeviation

    type(t_boxmc_8_10) :: bmc

    call init_mpi_data_parameters(comm)
    call bmc%init(comm)

    call bmc%get_coeff(comm,op_bg,src,ldir,phi0,theta0,dx,dy,dz, &
      S, T, S_tol, T_tol, inp_atol, inp_rtol )

    end subroutine

end module

program main
  use m_py_boxmc
  use m_data_parameters, only : ireals, iintegers, mpiint
  use mpi

  real(ireals) :: op_bg(3) = [1e-3, 1e-3, .5]
  real(ireals) :: phi0 = 0
  real(ireals) :: theta0 = 0
  integer(iintegers) :: src = 1
  logical :: ldir=.True.
  real(ireals) :: dx=100, dy=100, dz=50
  real(ireals) :: S(10)
  real(ireals) :: T(8)
  real(ireals) :: S_tol(10)
  real(ireals) :: T_tol(8)
  real(ireals) :: atol=1e-2_ireals, rtol=1e-1_ireals

  integer(mpiint) :: ierr

  call mpi_init(ierr)

  call get_coeff_8_10( MPI_COMM_WORLD, &
        op_bg, src,            &
        ldir, phi0, theta0,    &
        dx, dy, dz,            &
        S, T, S_tol,T_tol,     &
        atol, rtol)

  print *,'T', T
  print *,'S', S

  call mpi_finalize(ierr)
end program
