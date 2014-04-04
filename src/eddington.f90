module eddington
      use data_parameters, only: ireals,iintegers
      implicit none

      private
      public :: rodents

      logical,parameter :: ldelta_scale=.False.

      contains

      pure subroutine delta_scale(dtau,w0,g, dtau_d,g_d,w0_d)
        real(ireals),intent(in) :: g,w0,dtau
        real(ireals),intent(out) :: dtau_d,g_d,w0_d
        real(ireals) :: f

        if(ldelta_scale) then
          f = g**2
          dtau_d = dtau * ( 1. - w0 * f )
          g_d    = ( g - f ) / ( 1. - f )
          w0_d   = w0 * ( 1. - f ) / ( 1. - f * w0 )
        else
          dtau_d = dtau
          w0_d   = w0
          g_d    = g
        endif
      end subroutine

      pure subroutine rodents(dtau,w0,g,mu0, coeff)
        real(ireals),intent(out) :: coeff(5)
        real(ireals),intent(in) :: g,w0,dtau,mu0
        real(ireals) :: g_d,w0_d,dtau_d

        real(ireals) :: mu_0_inv,b_mmu_0,bscr,alpha_1,alpha_2,lambd, &
          exp1,term1,A,a11,a12,alpha_3,alpha_4,den1,   &
          alpha_5,alpha_6,a33,a13,a23

        call delta_scale(dtau,w0,g, dtau_d,w0_d,g_d)

        mu_0_inv = 1.0_ireals/mu0
        b_mmu_0 = 0.5_ireals - 0.75_ireals * g_d * mu0

        bscr = 0.5_ireals - 0.375_ireals * g_d
        alpha_1 = 2._ireals * ( 1._ireals - w0_d * ( 1._ireals - bscr ) ) - 0.25_ireals
        alpha_2 = 2._ireals * w0_d * bscr - 0.25_ireals

        lambd = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 )

        exp1  = exp( lambd * dtau_d )
        term1 = alpha_2 / ( alpha_1 - lambd ) * exp1
        A = 1.0_ireals / ( term1 - 1._ireals / term1 )
        a11 = A * 2.0_ireals * lambd / alpha_2
        a12 = A * ( exp1 - 1._ireals / exp1 )

        if(dtau_d.gt.100_ireals) then
          a11 = 0._ireals
          a12 = ( alpha_1 - lambd ) / alpha_2
        else
          exp1  = exp( lambd * dtau_d );
          term1 = alpha_2 / ( alpha_1 - lambd ) * exp1;

          A = 1.0_ireals / ( term1 - 1._ireals / term1 );

          a11 = A * 2.0_ireals * lambd / alpha_2;
          a12 = A * ( exp1 - 1._ireals / exp1 );
        endif

        alpha_3 = - w0_d * b_mmu_0
        alpha_4 = w0_d + alpha_3

        den1 = 1._ireals / ( mu_0_inv**2 - lambd * lambd )
        alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1
        alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1

        a33 = exp ( - dtau_d  * mu_0_inv )
        a13 = + alpha_5 * ( 1.0_ireals - a33 * a11 ) - alpha_6 * a12
        a23 = - alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 )

        coeff = [a11,a12,a13,a23,a33]
      end subroutine

  end module
