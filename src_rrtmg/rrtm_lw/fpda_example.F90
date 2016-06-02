program main

      use parkind, only: im => kind_im, rb => kind_rb

      use m_fpda_rrtm_lw
      implicit none

      integer(im),parameter :: nlay=20
      real(rb),dimension(1,nlay+1) :: plev
      real(rb),dimension(1,nlay)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
      real(rb),dimension(1,nlay)   :: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr

      integer(im) :: nbands
      real(rb),allocatable,target,dimension(:)   :: band_lbound,band_ubound
      real(rb),allocatable,target,dimension(:,:) :: tau_gas,planck_fracs!,planck_lev

      integer(im) :: k,iv
      real(rb) :: tmp,tmp_tot, planck

      do k=1,nlay+1
        plev(1,k)=1e3 - (k-1)*1e3/nlay
      enddo
      do k=1,nlay
        tlay(1,k)=300! - (k-1)*50/nlay
      enddo

      h2ovmr   = 9e-6
      o3vmr    = 5e-9
      co2vmr   = 400e-6
      ch4vmr   = 10e-6
      n2ovmr   = 320e-9
      o2vmr    = .209
      cfc11vmr = 0
      cfc12vmr = 0
      cfc22vmr = 0
      ccl4vmr  = 0

      call fpda_rrtm_lw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
              cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,                                     &
              nbands, band_lbound,band_ubound, tau_gas, planck_fracs)

      tmp_tot=0
      planck=0
      do iv=1,nbands
        tmp=0
        do k=nlay,nlay
          tmp = tmp + planck_fracs(k,iv)
!          planck = planck + planck_lev(k,iv)
        enddo
        tmp_tot=tmp_tot+tmp
        print *,iv,'::',1e4/band_lbound(iv),1e4/band_ubound(iv),tmp,':: sum ',tmp_tot,'planck',planck
      enddo

end program
