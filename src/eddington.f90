module eddington
      use data_parameters, only: ireals,iintegers,zero,one
      implicit none

      private
      public :: rodents

      logical,parameter :: ldelta_scale=.True.

      contains

      pure subroutine delta_scale(dtau,w0,g, dtau_d,g_d,w0_d)
        real(ireals),intent(in) :: g,w0,dtau
        real(ireals),intent(out) :: dtau_d,g_d,w0_d
        real(ireals) :: f

        if(ldelta_scale) then
          f = g**2
          dtau_d = dtau * ( one - w0 * f )
          g_d    = ( g - f ) / ( one - f )
          w0_d   = w0 * ( one - f ) / ( one - f * w0 )
        else
          dtau_d = dtau
          w0_d   = w0
          g_d    = g
        endif
      end subroutine

     subroutine rodents(dtau,w0,g,mu0, coeff)
        real(ireals),intent(out) :: coeff(5)
        real(ireals),intent(in) :: g,w0,dtau,mu0
        real(ireals) :: g_d,w0_d,dtau_d

        real(ireals) :: mu_0_inv,b_mmu_0,bscr,alpha_1,alpha_2,lambd, &
          exp1,term1,A,a11,a12,alpha_3,alpha_4,den1,   &
          alpha_5,alpha_6,a33,a13,a23


        
        call delta_scale(dtau,min( w0,one-1e-6_ireals ),g, dtau_d, w0_d, g_d)

        mu_0_inv = one/mu0
        b_mmu_0 = 0.5_ireals - 0.75_ireals * g_d * mu0

        bscr = 0.5_ireals - 0.375_ireals * g_d
        alpha_1 = 2._ireals * ( one - w0_d * ( one - bscr ) ) - 0.25_ireals
        alpha_2 = 2._ireals * w0_d * bscr - 0.25_ireals

        lambd = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 )

        exp1  = exp( lambd * dtau_d )
        term1 = alpha_2 / ( alpha_1 - lambd ) * exp1
        A = one / ( term1 - one / term1 )
        a11 = A * 2.0_ireals * lambd / alpha_2
        a12 = A * ( exp1 - one / exp1 )

        if(dtau_d.gt.30_ireals) then
          a11 = zero
          a12 = ( alpha_1 - lambd ) / alpha_2
        else
          exp1  = exp( lambd * dtau_d );
          term1 = alpha_2 / ( alpha_1 - lambd ) * exp1;

          A = one / ( term1 - one / term1 );

          a11 = A * 2.0_ireals * lambd / alpha_2;
          a12 = A * ( exp1 - one / exp1 );
        endif

        alpha_3 = - w0_d * b_mmu_0
        alpha_4 = w0_d + alpha_3

        den1 = one / ( mu_0_inv**2 - lambd * lambd )
        alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1
        alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1

        a33 = exp ( - dtau_d  * mu_0_inv )
        a13 = + alpha_5 * ( one - a33 * a11 ) - alpha_6 * a12
        a23 = - alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 )

        coeff = [a11,a12,a13,a23,a33]
!        call eddington_coeffc (dtau, g, w0, mu0, coeff)
        !todo twostream coefficients are neither energy conservant nor are they
        !restricted to positive values :(
        coeff = min( max(zero, coeff),one)

        if(dtau.lt.1e-4_ireals) coeff = [exp(-dtau),zero,zero,zero,exp(-dtau)] ! if theres no optical depth theres not much happening anyway
        if(sum(coeff(3:4) ).gt.one-coeff(5)) coeff(3:4) = coeff(3:4)/sum(coeff(3:4)) * (one-exp(-dtau*w0/mu0) ) ! dirty hack to conserve energy for dir2diff coeffs i.e. norm the coefficients to the direct scattering path extinction
        if(sum(coeff(1:2) ).gt.one) coeff(1:2) = coeff(1:2)/sum(coeff(1:2)) *.99_ireals
        
!        where (coeff.lt.zero .and. coeff.gt.-1e-3_ireals )
!          coeff = zero
!        end where
        if(any(coeff.lt.zero).or.any(coeff.gt.one)) then
          print *,'dtau,w0,g,mu0',dtau,w0,g,mu0,'==>', coeff
          stop 'error eddington coeffs - bigger one or smaller zero!'
        endif
        if(sum(coeff(1:2) ).gt.one) then
          print *,'dtau,w0,g,mu0',dtau,w0,g,mu0,'==>', coeff
          stop 'error eddington coeffs - diffuse streams sum bigger one - this is not compliant to energy conservation!'
        endif
!        if(sum(coeff(3:4) ).gt.one-coeff(5)) then
!          print *,'dtau,w0,g,mu0',dtau,w0,g,mu0,'==>', coeff
!          stop 'error eddington coeffs - dir2diff streams sum bigger 1-Transmission - this is not compliant to energy conservation!'
!        endif
      end subroutine
      subroutine eddington_coeffc (dtau, g, omega0, mu0, coeff)
          real(ireals),intent(out) :: coeff(5)
          real(ireals),intent(in) :: g,omega0,dtau,mu0

          real(ireals) :: a11,a12,a13,a23,a33
          real(ireals) :: alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
          real(ireals) :: lambda, b, A
          real(ireals) :: denom

          alpha1= (one-omega0)+0.75*(one-omega0*g);
          alpha2=-(one-omega0)+0.75*(one-omega0*g);

          lambda=sqrt(alpha1*alpha1-alpha2*alpha2);

          if(dtau.gt.100._ireals) then
            a11 = zero
            a12 = ( alpha1 - lambda ) / alpha2
          else
            A=one/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));

            a11=A*2.0*lambda/alpha2;
            a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));
          endif

          b=0.5-0.75*g*mu0;
          alpha3=-omega0*b;
          alpha4=omega0*(1-b);

          denom = (one/mu0/mu0-lambda*lambda);
          alpha5=((alpha1-one/mu0)*alpha3-alpha2*alpha4)/denom;
          alpha6=(alpha2*alpha3-(alpha1+one/mu0)*alpha4)/denom;

          a33=exp(-dtau/mu0);

          a13=alpha5*(one-(a11)*(a33))-alpha6*(a12);
          a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11));
          coeff = [a11,a12,a13,a23,a33]
      end subroutine


  end module
