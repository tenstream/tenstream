module m_fpda_rrtm_lw
      use rrtmg_lw_init, only: rrtmg_lw_ini
      use parkind, only: im => kind_im, rb => kind_rb
      use rrlw_wvn
      use parrrtm, only: nbndlw
      use rrtmg_lw_rad, only: rrtmg_lw, fpda_taugas,fpda_planck_fracs,fpda_planck_lev
      
      use iso_c_binding
      implicit none

      logical :: linit=.False.


    contains
      subroutine f2c_fpda_rrtm_lw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,&
              cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,                                              &
              c_nbands,c_band_lbound,c_band_ubound,c_planck_fracs, c_tau_gas) bind(C,name="fpda_rrtm_lw")

          integer(c_int),intent(in) :: nlay

          real(c_double),dimension(1,nlay+1),intent(in) :: plev
          real(c_double),dimension(1,nlay),intent(in)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
          real(c_double),dimension(1,nlay),intent(in)   :: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr

          integer(c_int),intent(out) :: c_nbands
          type(c_ptr),intent(out) :: c_band_lbound,c_band_ubound
          type(c_ptr),intent(out) :: c_tau_gas,c_planck_fracs

          real(rb),dimension(1,nlay+1) :: plev_r
          real(rb),dimension(1,nlay)   :: tlay_r, h2ovmr_r, o3vmr_r, co2vmr_r, ch4vmr_r, n2ovmr_r, o2vmr_r
          real(rb),dimension(1,nlay)   :: cfc11vmr_r,cfc12vmr_r,cfc22vmr_r,ccl4vmr_r

          integer(im) :: nbands
          real(rb),allocatable,save,target,dimension(:)   :: band_lbound,band_ubound         ! [      nbands ]
          real(rb),allocatable,save,target,dimension(:,:) :: tau_gas,planck_fracs            ! [ nlay,nbands ]

          integer(im) :: ib

          plev_r         (1,:)= rev( plev         (1,:) )
          tlay_r         (1,:)= rev( tlay         (1,:) )
          h2ovmr_r       (1,:)= rev( h2ovmr       (1,:) )
          o3vmr_r        (1,:)= rev( o3vmr        (1,:) )
          co2vmr_r       (1,:)= rev( co2vmr       (1,:) )
          ch4vmr_r       (1,:)= rev( ch4vmr       (1,:) )
          n2ovmr_r       (1,:)= rev( n2ovmr       (1,:) )
          o2vmr_r        (1,:)= rev( o2vmr        (1,:) )
          cfc11vmr_r     (1,:)= rev( cfc11vmr     (1,:) )
          cfc12vmr_r     (1,:)= rev( cfc12vmr     (1,:) )
          cfc22vmr_r     (1,:)= rev( cfc22vmr     (1,:) )
          ccl4vmr_r      (1,:)= rev( ccl4vmr      (1,:) )

          call fpda_rrtm_lw(nlay,plev_r,tlay_r,                           &
              h2ovmr_r, o3vmr_r, co2vmr_r, ch4vmr_r, n2ovmr_r, o2vmr_r,   &
              cfc11vmr_r,cfc12vmr_r,cfc22vmr_r,ccl4vmr_r,                 &
              nbands, band_lbound,band_ubound, tau_gas, planck_fracs)

          do ib=1,nbands
            tau_gas        (:,ib)= rev( tau_gas      (:,ib) )   
            planck_fracs   (:,ib)= rev( planck_fracs (:,ib) )        
          enddo

          c_nbands = ngptlw
          c_band_lbound = c_loc(band_lbound)
          c_band_ubound = c_loc(band_ubound)
          c_tau_gas     = c_loc(tau_gas)
          c_planck_fracs= c_loc(planck_fracs)
      end subroutine

      subroutine fpda_rrtm_lw(nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
              cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,                                           &
              nbands, band_lbound,band_ubound, tau_gas, planck_fracs)

        integer(im),intent(in) :: nlay

        real(rb),dimension(:,:),intent(in) :: plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
        real(rb),dimension(:,:),intent(in) :: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr     

        integer(im),intent(out) :: nbands
        real(rb),allocatable,target,dimension(:),intent(out)   :: band_lbound,band_ubound
        real(rb),allocatable,target,dimension(:,:),intent(out) :: tau_gas,planck_fracs


        integer(im),parameter :: ncol=1,inflglw=0,iceflglw=0,liqflglw=0
        integer(kind=im) :: icld=0         ! Cloud overlap method
        integer(kind=im) :: idrv=0         ! Flag for calculation of dFdT

        real(rb),dimension(ncol,nlay+1) :: tlev,hhl

        real(rb),dimension(ncol,nlay  ) :: play, cldfr, cicewp, cliqwp, reice, reliq 

        real(rb),dimension(nbndlw, ncol,nlay  ) :: taucld
        real(rb),dimension(ncol, nlay, nbndlw ) :: tauaer
        real(rb),dimension(ncol, nbndlw ) :: emis


        real(rb),dimension(ncol     ) :: tsfc 

        real(rb),dimension(ncol,nlay+1) :: uflx,dflx,uflxc,dflxc
        real(rb),dimension(ncol,nlay  ) :: hr,hrc

        integer(im) :: iv,ig,k,icol,iw


        do icol=1,ubound(plev,1)
          call hydrostat_lev(plev(icol,:),tlay(icol,:), hhl(icol,:))
          play(icol,:)      = .5_rb*(plev(icol,1:nlay)+plev(icol,2:nlay+1))

          tlev(icol,nlay+1) = tlay(icol,nlay)
          tlev(icol,1)      = tlay(icol,1)
          tsfc(icol)        = tlay(icol,1)

          do k = 1,nlay-1
            tlev(icol,k+1)  = .5_rb * (tlay(icol,k+1)+tlay(icol,k) )
          enddo

        enddo

        cldfr  = 0; taucld = 0; emis   = 1
        cicewp = 0; cliqwp = 0; reice  = 0;
        reliq  = 0; tauaer = 0;
        

        if(.not.linit) then
          call rrtmg_lw_ini(1006._rb)
          linit = .True.

!          print *,'hhl',hhl(1,:)
!          print *,'tlev',tlev(1,:)
!          print *,'tlay',tlay(1,:)
!          print *,'plev',plev(1,:)
!          print *,'play',play(1,:)
        endif

      call rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfr   , &
             taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc)

        nbands = ngptlw
        if(.not.allocated(band_lbound)) allocate(band_lbound (ngptlw) )
        if(.not.allocated(band_ubound)) allocate(band_ubound (ngptlw) )
        if(.not.allocated(tau_gas    )) allocate(tau_gas     (nlay,ngptlw) )
        if(.not.allocated(planck_fracs))allocate(planck_fracs(nlay,ngptlw) )
!        if(.not.allocated(planck_lev  ))allocate(planck_lev  (nlay+1,ngptlw) )

        iw=0
        do iv=1,nbndlw
          do ig=1,ngc(iv)
            iw = iw+1
            band_lbound(iw) = wavenum1(iv)
            band_ubound(iw) = wavenum2(iv)
            tau_gas(:,iw)   = fpda_taugas(:,iw)
            planck_fracs(:,iw)   = fpda_planck_fracs(:,iw)
!            planck_lev  (:,iw)   = fpda_planck_lev  (:,iv)
          enddo
        enddo

!        do k=nlay+1,1,-1
!          print *,'flx',k,'::',uflxc(1,k),':',dflxc(1,k)
!        enddo
      end subroutine

      subroutine hydrostat_lev(plev,tlay, hhl)
        real(rb),intent(in) :: plev(:),tlay(:)
        real(rb),intent(out) :: hhl(size(plev))
        real(rb) :: dp,dz,rho
        integer(im) :: k
        hhl(1) = 0._rb
        do k=1,size(tlay)
          dp  = abs( plev(k)-plev(k+1) ) 
          rho = ( plev(k)+dp/2._rb  ) / 287.058_rb / tlay(k)
          dz  = dp / rho / 9.8065_rb
          hhl(k+1) = hhl(k) + dz
        enddo
      end subroutine

      function rev(inp) 
          real(rb),intent(in) :: inp(:)
          real(rb) :: rev(size(inp))
          rev = inp(ubound(inp,1):lbound(inp,1):-1)
      end function

end module

