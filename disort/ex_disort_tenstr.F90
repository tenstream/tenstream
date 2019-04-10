program main
  use m_tenstr_disort, only: default_flx_computation
  implicit none
  integer, parameter :: nlyr=5, nstreams=16
  real, parameter :: mu0=.5, S0=1000, Ag=.2
  real, dimension(nlyr) :: Bfracs, dtau, ssalb, gasym
  real, dimension(nlyr+1) :: temper, RFLDIR, RFLDN, FLUP, DFDT, UAVG
  integer :: k

  dtau = 1
  ssalb = 0
  gasym = .5
  temper = 288
  Bfracs = .5

  call solar()
  call thermal()
contains
  subroutine thermal()
    Bfracs = .1
    call default_flx_computation(&
      0., 0., Ag, &
      .True., [1., 10000.], Bfracs, &
      dtau, ssalb, gasym, temper, &
      RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
      nstreams, lverbose=.False.)
    print *,'----------- THERMAL COMPUTATION 10% source term -----------'
    print *,' z-idx         '// &
      'RFLDIR           RFLDN            FLUP              DFDT              UAVG'
    do k=1,nlyr+1
      print *, k, RFLDIR(k), RFLDN(k), FLUP(k), DFDT(k), UAVG(k)
    enddo

    Bfracs = 1
    call default_flx_computation(&
      0., 0., Ag, &
      .True., [1., 10000.], Bfracs, &
      dtau, ssalb, gasym, temper, &
      RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
      nstreams, lverbose=.False.)

    print *,'----------- THERMAL COMPUTATION 100% source term ----------'
    print *,' z-idx         '// &
      'RFLDIR           RFLDN            FLUP              DFDT              UAVG'
    do k=1,nlyr+1
      print *, k, RFLDIR(k), RFLDN(k), FLUP(k), DFDT(k), UAVG(k)
    enddo
  end subroutine
  subroutine solar()
    real, dimension(nlyr+1) :: Transmission
    call default_flx_computation(&
      mu0, S0, Ag, &
      .False., [0., 0.], Bfracs, &
      dtau, ssalb, gasym, temper, &
      RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
      nstreams, lverbose=.False.)

    Transmission(1) = 1.
    do k=1,nlyr
      Transmission(k+1) = Transmission(k) * exp(-dtau(k)/mu0)
    enddo

    print *,'----------- SOLAR COMPUTATION -------------'
    print *,' z-idx         Transmission       '// &
      'RFLDIR           RFLDN            FLUP              DFDT              UAVG'
    do k=1,nlyr+1
      print *, k, Transmission(k)*S0*mu0, RFLDIR(k), RFLDN(k), FLUP(k), DFDT(k), UAVG(k)
    enddo
  end subroutine
end program
