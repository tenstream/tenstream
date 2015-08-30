@test
subroutine test_eddington_zdun()

   use m_eddington
   use m_data_parameters, only: ireals,iintegers,zero,one,pi

   use pfunit_mod


   implicit none

   real(ireals),parameter :: tol=1e-6

   real(ireals) :: tau, w0, g, mu0
   real(ireals) :: a11,a12,a13,a23,a33,g1,g2

   tau = zero
   w0  = one
   g   = one
   mu0 = one

   call eddington_coeff_zdun(tau,w0,g,mu0,a11,a12,a13,a23,a33,g1,g2)
   print *,a11,a12,a13,a23,a33,g1,g2

!   @assertAll([one,zero,zero],[a11,a12,a13], tol)

   @assertEqual(one , a11,  tol)
   @assertEqual(zero, a12,  tol)
   @assertEqual(zero, a13,  tol)
   @assertEqual(zero, a23,  tol)
   @assertEqual(one , a33,  tol)

!   @assertEqual(g1 , one , tol) ! TODO: dont know what exactly this value should be?!
!   @assertEqual(g2 , zero, tol) ! TODO: dont know what exactly this value should be?! 

end subroutine
