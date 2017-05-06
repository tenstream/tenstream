module m_py_eddington
  use m_data_parameters, only : ireals, iintegers
  use m_eddington, only : eddington_coeff_zdun

  implicit none

  private
  public :: py_eddington_coeff_zdun

  contains

     subroutine py_eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33,g1,g2)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: a11,a12,a13,a23,a33,g1,g2

          call eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33,g1,g2)
    end subroutine

end module
