@test(npes =[1])
subroutine test_eddington(this)

    use m_schwarzschild, only: schwarzschild, B_eff
    use m_data_parameters, only: ireals, iintegers, pi
    use m_helper_functions, only: is_between

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: Nlay=10
    integer(iintegers) :: k, itau
    real(ireals) :: Blev(Nlay+1), B, tau

    do k=1, Nlay+1
      Blev(k) = 5.67e-8_ireals * real(280+k, ireals)**4 / pi
    enddo

    tau = 1000._ireals
    do k=1, Nlay
      call B_eff(Blev(k), Blev(k+1), tau, B)
      !print *,k,Blev(k), Blev(k+1), B, is_between(B,Blev(k), Blev(k+1))
      @assertEqual(Blev(k+1), B, .1_ireals, 'for very high optical depth, the values should be approx the same')
      @assertTrue(is_between(B, Blev(k+1), Blev(k)), 'resulting B_eff has to be between the far and near value')
    enddo

    tau = 0._ireals
    do k=1, Nlay
      call B_eff(Blev(k), Blev(k+1), tau, B)
      !print *,k,Blev(k), Blev(k+1), B, is_between(B,Blev(k), Blev(k+1))
      @assertEqual((Blev(k)+Blev(k+1))/2, B, .1_ireals, 'for zero optical depth, should be approx the mean')
      @assertTrue(is_between(B, Blev(k+1), Blev(k)), 'resulting B_eff has to be between the far and near value')
    enddo

    do itau=-100,40
      tau = 10._ireals**(real(itau, ireals)/10)
      do k=1, Nlay
        call B_eff(Blev(k), Blev(k+1), tau, B)
        !print *,tau, k,Blev(k), Blev(k+1), B, is_between(B,Blev(k), Blev(k+1))
        @assertTrue(is_between(B, Blev(k+1), Blev(k)), 'resulting B_eff has to be between the far and near value')
      enddo
    enddo

    ! make sure that this also works in reverse
    do itau=-100,40
      tau = 10._ireals**(real(itau, ireals)/10)
      do k=1, Nlay
        call B_eff(Blev(k+1), Blev(k), tau, B)
        !print *,tau, k,Blev(k), Blev(k+1), B, is_between(B,Blev(k), Blev(k+1))
        @assertTrue(is_between(B, Blev(k+1), Blev(k)), 'resulting B_eff has to be between the far and near value')
      enddo
    enddo
end subroutine

