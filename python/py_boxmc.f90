module m_py_boxmc
  use m_data_parameters, only : ireals, ireal_dp, mpiint, iintegers, &
                                init_mpi_data_parameters
  use m_boxmc, only : t_boxmc, t_boxmc_3_10, t_boxmc_8_10, t_boxmc_wedge_5_8
  use m_boxmc_geometry, only:  setup_default_unit_cube_geometry, setup_default_unit_wedge_geometry
  use mpi

  implicit none

  private
  public :: py_get_coeff_3_10, py_get_coeff_8_10, py_get_coeff_wedge_5_8

contains

  subroutine py_get_coeff_8_10( &
      comm, dx, dy, dz,      &
      op_bg, src,            &
      ldir, phi0, theta0,    &
      ret_S, ret_T,          &
      ret_S_tol,ret_T_tol,   &
      inp_atol, inp_rtol)
    double precision,intent(in)  :: dx, dy, dz   !< @param[in] size of cube
    double precision,intent(in)  :: op_bg(3)     !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
    double precision,intent(in)  :: phi0         !< @param[in] phi0 solar azimuth angle
    double precision,intent(in)  :: theta0       !< @param[in] theta0 solar zenith angle
    integer,intent(in)           :: src                         !< @param[in] src stream from which to start photons - see init_photon routines
    integer(mpiint),intent(in)   :: comm         !< @param[in] comm MPI Communicator
    logical,intent(in)           :: ldir         !< @param[in] ldir determines if photons should be started with a fixed incidence angle
    double precision,intent(out) :: ret_S(10)    !< @param[out] S_out diffuse streams transfer coefficients
    double precision,intent(out) :: ret_T(8)     !< @param[out] T_out direct streams transfer coefficients
    double precision,intent(out) :: ret_S_tol(10)!< @param[out] absolute tolerances of results
    double precision,intent(out) :: ret_T_tol(8) !< @param[out] absolute tolerances of results
    double precision,intent(in)  :: inp_atol     !< @param[in] inp_atol if given, determines targeted absolute stddeviation
    double precision,intent(in)  :: inp_rtol     !< @param[in] inp_rtol if given, determines targeted relative stddeviation

    type(t_boxmc_8_10) :: bmc

    real(ireals)       :: S(size(ret_S))
    real(ireals)       :: T(size(ret_T))
    real(ireals)       :: S_tol(size(ret_S_tol))
    real(ireals)       :: T_tol(size(ret_T_tol))

    real(ireals), allocatable :: vertices(:)

    call init_mpi_data_parameters(comm)
    call bmc%init(comm)

    call setup_default_unit_cube_geometry(&
      real(dx, ireals), real(dy, ireals), real(dz,ireals), &
      vertices)

    call bmc%get_coeff(comm,      &
      real(op_bg, kind=ireals),   &
      int(src,kind=iintegers),    &
      ldir,                       &
      real(phi0, kind=ireals),    &
      real(theta0, kind=ireals),  &
      real(vertices, kind=ireals), &
      S, T, S_tol, T_tol,         &
      real(inp_atol, kind=ireals),&
      real(inp_rtol, kind=ireals))

    ret_S = S
    ret_T = T
    ret_S_tol = S_tol
    ret_T_tol = T_tol
  end subroutine

  subroutine py_get_coeff_wedge_5_8( &
      comm, dx, dy, dz,      &
      op_bg, src,            &
      ldir, phi0, theta0,    &
      ret_S, ret_T,          &
      ret_S_tol,ret_T_tol,   &
      inp_atol, inp_rtol)
    double precision,intent(in)  :: dx, dy, dz   !< @param[in] size of cube
    double precision,intent(in)  :: op_bg(3)     !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
    double precision,intent(in)  :: phi0         !< @param[in] phi0 solar azimuth angle
    double precision,intent(in)  :: theta0       !< @param[in] theta0 solar zenith angle
    integer,intent(in)           :: src                         !< @param[in] src stream from which to start photons - see init_photon routines
    integer,intent(in)           :: comm         !< @param[in] comm MPI Communicator
    logical,intent(in)           :: ldir         !< @param[in] ldir determines if photons should be started with a fixed incidence angle
    double precision,intent(out) :: ret_S(8)     !< @param[out] S_out diffuse streams transfer coefficients
    double precision,intent(out) :: ret_T(5)     !< @param[out] T_out direct streams transfer coefficients
    double precision,intent(out) :: ret_S_tol(8) !< @param[out] absolute tolerances of results
    double precision,intent(out) :: ret_T_tol(5) !< @param[out] absolute tolerances of results
    double precision,intent(in)  :: inp_atol     !< @param[in] inp_atol if given, determines targeted absolute stddeviation
    double precision,intent(in)  :: inp_rtol     !< @param[in] inp_rtol if given, determines targeted relative stddeviation

    type(t_boxmc_wedge_5_8) :: bmc

    real(ireals)       :: S(size(ret_S))
    real(ireals)       :: T(size(ret_T))
    real(ireals)       :: S_tol(size(ret_S_tol))
    real(ireals)       :: T_tol(size(ret_T_tol))

    real(ireals), allocatable :: vertices(:)

    call init_mpi_data_parameters(comm)
    call bmc%init(comm)

    call setup_default_unit_wedge_geometry( &
      real(dx, ireals), real(dy, ireals), real(dz,ireals), &
      vertices)

    call bmc%get_coeff(comm,      &
      real(op_bg, kind=ireals),   &
      int(src,kind=iintegers),    &
      ldir,                       &
      real(phi0, kind=ireals),    &
      real(theta0, kind=ireals),  &
      real(vertices, kind=ireals), &
      S, T, S_tol, T_tol,         &
      real(inp_atol, kind=ireals),&
      real(inp_rtol, kind=ireals))

    ret_S = S
    ret_T = T
    ret_S_tol = S_tol
    ret_T_tol = T_tol
  end subroutine

  subroutine py_get_coeff_3_10( &
      comm, dx, dy, dz,      &
      op_bg, src,            &
      ldir, phi0, theta0,    &
      ret_S, ret_T,          &
      ret_S_tol,ret_T_tol,   &
      inp_atol, inp_rtol)
    double precision,intent(in)  :: dx, dy, dz   !< @param[in] size of cube
    double precision,intent(in)  :: op_bg(3)     !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
    double precision,intent(in)  :: phi0         !< @param[in] phi0 solar azimuth angle
    double precision,intent(in)  :: theta0       !< @param[in] theta0 solar zenith angle
    integer,intent(in)           :: src                         !< @param[in] src stream from which to start photons - see init_photon routines
    integer(mpiint),intent(in)   :: comm         !< @param[in] comm MPI Communicator
    logical,intent(in)           :: ldir         !< @param[in] ldir determines if photons should be started with a fixed incidence angle
    double precision,intent(out) :: ret_S(10)    !< @param[out] S_out diffuse streams transfer coefficients
    double precision,intent(out) :: ret_T(3)     !< @param[out] T_out direct streams transfer coefficients
    double precision,intent(out) :: ret_S_tol(10)!< @param[out] absolute tolerances of results
    double precision,intent(out) :: ret_T_tol(3) !< @param[out] absolute tolerances of results
    double precision,intent(in)  :: inp_atol     !< @param[in] inp_atol if given, determines targeted absolute stddeviation
    double precision,intent(in)  :: inp_rtol     !< @param[in] inp_rtol if given, determines targeted relative stddeviation

    type(t_boxmc_3_10) :: bmc

    real(ireals)       :: S(size(ret_S))
    real(ireals)       :: T(size(ret_T))
    real(ireals)       :: S_tol(size(ret_S_tol))
    real(ireals)       :: T_tol(size(ret_T_tol))

    real(ireals), allocatable :: vertices(:)

    call init_mpi_data_parameters(comm)
    call bmc%init(comm)

    call setup_default_unit_cube_geometry(&
      real(dx, ireals), real(dy, ireals), real(dz,ireals), &
      vertices)

    call bmc%get_coeff(comm,      &
      real(op_bg, kind=ireals),   &
      int(src,kind=iintegers),    &
      ldir,                       &
      real(phi0, kind=ireals),    &
      real(theta0, kind=ireals),  &
      real(vertices, kind=ireals), &
      S, T, S_tol, T_tol,         &
      real(inp_atol, kind=ireals),&
      real(inp_rtol, kind=ireals))

    ret_S = S
    ret_T = T
    ret_S_tol = S_tol
    ret_T_tol = T_tol
  end subroutine

end module

program main
  use m_py_boxmc
  use m_data_parameters, only : ireal_dp, iintegers, mpiint
  use mpi

  real(ireal_dp) :: dx=100, dy=100, dz=100
  real(ireal_dp) :: op_bg(3) = [1e-3, 1e-3, .5]
  real(ireal_dp) :: phi0 = 0
  real(ireal_dp) :: theta0 = 0
  integer(iintegers) :: src = 1
  logical :: ldir=.True.
  real(ireal_dp) :: S(10)
  real(ireal_dp) :: T3(3), T(8)
  real(ireal_dp) :: S_tol(10)
  real(ireal_dp) :: T_tol3(3), T_tol(8)
  real(ireal_dp) :: atol=1e-2_ireal_dp, rtol=1e-1_ireal_dp

  integer(mpiint) :: ierr

  call mpi_init(ierr)

  call py_get_coeff_8_10(MPI_COMM_WORLD, &
        dx, dy, dz,            &
        op_bg, src,            &
        ldir, phi0, theta0,    &
        S, T, S_tol,T_tol,     &
        atol, rtol)

  print *,'T', T
  print *,'S', S

  call py_get_coeff_3_10(MPI_COMM_WORLD, &
        dx, dy, dz,            &
        op_bg, src,            &
        ldir, phi0, theta0,    &
        S, T3, S_tol,T_tol3,   &
        atol, rtol)

  print *,'T3_10', T3
  print *,'S3_10', S

  call mpi_finalize(ierr)
end program
