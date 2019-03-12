module m_py_param_phi
  use m_data_parameters, only : ireals, irealLUT, mpiint, iintegers, &
                                init_mpi_data_parameters
  use mpi
  use m_LUT_param_phi, only: iterative_phi_theta_from_param_phi_and_param_theta
  use m_boxmc_geometry, only: setup_default_wedge_geometry

  implicit none

  private
  public :: py_iterative_phi_theta_from_param_phi_and_param_theta

  contains

    subroutine py_iterative_phi_theta_from_param_phi_and_param_theta( &
        sphere_radius, dz, Cx, Cy, &
        param_phi, param_theta, pyphi, pytheta)

    double precision, intent(in)  :: param_phi, param_theta, dz, sphere_radius, Cx, Cy
    double precision, intent(out) :: pyphi, pytheta

    real(ireals), allocatable :: vertices(:)
    real(irealLUT) :: phi, theta
    integer(mpiint) :: ierr

    call setup_default_wedge_geometry(&
      [0._ireals, 0._ireals], &
      [1._ireals, 0._ireals], &
      real([Cx, Cy], ireals), &
      real(dz, ireals), vertices, &
      sphere_radius=real(sphere_radius, ireals))

    call iterative_phi_theta_from_param_phi_and_param_theta(&
      real(vertices, irealLUT), &
      real(param_phi, irealLUT), real(param_theta, irealLUT), &
      phi, theta, ierr)

    pyphi = phi
    pytheta = theta
    end subroutine

end module
