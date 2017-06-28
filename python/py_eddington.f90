module m_py_eddington
  use m_data_parameters, only : ireals, iintegers
  use m_eddington, only : eddington_coeff_zdun

  implicit none

  private
  public :: py_eddington_coeff_zdun

  contains

     subroutine py_eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33,g1,g2)
       double precision,intent(in) :: dtau_in,g_in,omega_0_in,mu_0
       double precision,intent(out) :: a11,a12,a13,a23,a33,g1,g2
       real(ireals) :: a11i,a12i,a13i,a23i,a33i,g1i,g2i

       call eddington_coeff_zdun(real(dtau_in, kind=ireals),        &
                                 real(omega_0_in, kind=ireals),     &
                                 real(g_in, kind=ireals),           &
                                 real(mu_0, kind=ireals),           &
                                 a11i,a12i,a13i,a23i,a33i,g1i,g2i)

      a11 = a11i
      a12 = a12i
      a13 = a13i
      a23 = a23i
      a33 = a33i
      g1  = g1i
      g2  = g2i
    end subroutine

end module
