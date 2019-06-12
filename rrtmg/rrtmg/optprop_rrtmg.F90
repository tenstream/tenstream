!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

!> \page Routines to call tenstream with optical properties from RRTM
!! Routine in this file call a hacked rrtmg radiative transfer code
!! which returns the optical properties in each band
!! or may also use the rrtmg solver to compute radiative fluxes

module m_optprop_rrtmg
  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb
  use m_tenstr_rrtmg_lw_init, only: rrtmg_lw_ini
  use m_tenstr_parrrtm, only: nbndlw
  use m_tenstr_rrtmg_lw_rad, only: rrtmg_lw

  use m_tenstr_rrtmg_sw_init, only: rrtmg_sw_ini
  use m_tenstr_parrrsw, only: nbndsw, naerec
  use m_tenstr_rrtmg_sw_rad, only: rrtmg_sw
  use m_tenstr_rrtmg_sw_spcvrt, only: tenstr_solsrc

  use m_data_parameters, only: iintegers, ireals, one, mpiint
  use m_helper_functions, only: deg2rad, CHKERR

  implicit none

  private
  public :: optprop_rrtm_lw, optprop_rrtm_sw, linit_rrtmg_lw

  logical,parameter :: ldebug=.False.
  logical,save :: linit_rrtmg_lw=.False.

contains
  subroutine optprop_rrtm_lw(ncol_in, nlay_in, &
      albedo, plev, tlev, tlay, tsrfc, &
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      lwp, reliq, iwp, reice, tau, Bfrac, opt_tau_f, &
      opt_lwuflx, opt_lwdflx, opt_lwhr, opt_cldfr)
    ! RRTM needs the arrays to start at the surface

    integer(iintegers),intent(in) :: ncol_in, nlay_in
    real(ireals), intent(in) :: albedo

    real(ireals),dimension(ncol_in,nlay_in+1), intent(in) :: plev, tlev
    real(ireals),dimension(ncol_in,nlay_in)  , intent(in) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
    real(ireals),dimension(ncol_in,nlay_in)  , intent(in) :: lwp, reliq, iwp, reice
    real(ireals),dimension(ncol_in)          , intent(in) :: tsrfc

    real(ireals), dimension(:,:,:), intent(out) :: tau, Bfrac ! [nlay, ncol, ngptlw]
    real(ireals), dimension(:,:,:), intent(out), optional :: opt_tau_f ! [nlay, ncol, ngptlw]
    real(ireals), dimension(:,:), intent(out), optional :: opt_lwuflx, opt_lwdflx, opt_lwhr ! [nlay+1, ncol]
    real(ireals), dimension(:,:), intent(in), optional :: opt_cldfr ! [ncol, nlay]

    real(rb),dimension(ncol_in,nlay_in) :: play, cldfr

    real(rb),dimension(nbndlw, ncol_in, nlay_in) :: taucld
    real(rb),dimension(ncol_in, nlay_in, nbndlw ) :: tauaer
    real(rb),dimension(ncol_in, nbndlw ) :: emis

    real(rb),dimension(ncol_in, nlay_in)   :: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr

    real(rb),dimension(ncol_in,nlay_in+1) :: lwuflx,lwdflx,lwuflxc,lwdflxc
    real(rb),dimension(ncol_in,nlay_in  ) :: lwhr,lwhrc

    integer(im) :: icol, ilay, ncol, nlay

    integer(im),parameter :: inflglw=2,liqflglw=1,iceflglw=3
    integer(kind=im) :: icld=2         ! Cloud overlap method
    integer(kind=im) :: idrv=0         ! Flag for calculation of dFdT

    ! copy from TenStream to RRTM precision:
    ncol   = int(ncol_in, kind=im)
    nlay   = int(nlay_in, kind=im)

    ! Take average pressure and temperature as mean values for voxels --
    ! should probably use log interpolation for pressure...
    do ilay=1,nlay
      do icol=1,ncol
        play(icol,ilay) = .5_rb*(plev(icol,ilay)+plev(icol,ilay+1))
      enddo
    enddo

    taucld   = 0;
    tauaer   = 0; lwdflxc  = 0; lwuflxc  = 0;
    cfc11vmr = 0; cfc12vmr = 0; cfc22vmr = 0; ccl4vmr = 0;

    emis = one - albedo

    if(present(opt_cldfr)) then
        cldfr = real(opt_cldfr, rb)
    else
      where(lwp.gt.0)
        cldfr = 1
      elsewhere
        cldfr = 0
      end where
    endif

    if(.not.linit_rrtmg_lw) then
      call rrtmg_lw_ini(1006._rb)
      linit_rrtmg_lw = .True.

      !if(ldebug .and. myid.eq.0) then
      !  do k=nlay,1,-1
      !    print *,'rrtm_optprop_lw',k,'tlev',tlev(1,k),'tlay',tlay(1,k),'plev',plev(1,k),'play',play(1,k), &
      !      'lwp',lwp(1,k), 'reliq',reliq(1,k), 'iwp',iwp(1,k), 'reice',reice(1,k),&
      !      'h2o',h2ovmr(1,k), 'o3' , o3vmr(1,k), 'co2', co2vmr(1,k), 'ch4', ch4vmr(1,k),    &
      !      'n2o', n2ovmr(1,k), 'o2' , o2vmr(1,k)
      !  enddo
      !endif
    endif

    if (present(opt_lwuflx).and.present(opt_lwdflx).and.present(opt_lwhr)) then
      call rrtmg_lw &
        (ncol, nlay, icld, idrv, &
        real(play,rb), real(plev,rb), &
        real(tlay,rb), real(tlev,rb), real(tsrfc, rb), &
        real(h2ovmr,rb), real(o3vmr,rb), real(co2vmr,rb), &
        real(ch4vmr,rb), real(n2ovmr,rb), real(o2vmr,rb), &
        cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
        inflglw, iceflglw, liqflglw, cldfr, taucld, &
        real(iwp, rb), real(lwp, rb), real(reice, rb), real(reliq, rb), &
        tauaer, &
        lwuflx, lwdflx  ,lwhr ,lwuflxc ,lwdflxc ,lwhrc, &
        tau, Bfrac, loptprop_only=.False., tenstr_tau_f=opt_tau_f)
      opt_lwuflx = transpose(real(lwuflx, ireals))
      opt_lwdflx = transpose(real(lwdflx, ireals))
      opt_lwhr   = transpose(real(lwhr, ireals))
    else
      call rrtmg_lw &
        (ncol, nlay, icld, idrv, &
        real(play,rb), real(plev,rb), &
        real(tlay,rb), real(tlev,rb), real(tsrfc, rb), &
        real(h2ovmr,rb), real(o3vmr,rb), real(co2vmr,rb), &
        real(ch4vmr,rb), real(n2ovmr,rb), real(o2vmr,rb), &
        cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
        inflglw, iceflglw, liqflglw, cldfr, taucld, &
        real(iwp, rb), real(lwp, rb), real(reice, rb), real(reliq, rb), &
        tauaer, &
        lwuflx, lwdflx  ,lwhr ,lwuflxc ,lwdflxc ,lwhrc, &
        tau, Bfrac, loptprop_only=.True., tenstr_tau_f=opt_tau_f)
    endif
  end subroutine

  subroutine optprop_rrtm_sw(ncol_in, nlay_in, &
      theta0, albedo, &
      plev, tlev, tlay, &
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      lwp, reliq, iwp, reice, tau, w0, g, &
      opt_solar_constant, opt_cldfr, &
      opt_swdirflx, opt_swuflx, opt_swdflx, opt_swhr, &
      opt_tau_f, opt_w0_f, opt_g_f)
    ! RRTM needs the arrays to start at the surface

    integer(iintegers),intent(in)          :: ncol_in, nlay_in
    real(ireals), intent(in) :: theta0, albedo

    real(ireals),dimension(ncol_in,nlay_in+1), intent(in) :: plev, tlev
    real(ireals),dimension(ncol_in,nlay_in)  , intent(in) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr
    real(ireals),dimension(ncol_in,nlay_in)  , intent(in) :: lwp, reliq, iwp, reice

    real(ireals), dimension(:,:,:), intent(out) :: tau, w0, g ! [nlay, ncol, ngptsw]

    real(ireals), intent(in), optional :: opt_solar_constant
    real(ireals),dimension(:,:), intent(in), optional :: opt_cldfr ! [ncol, nlay]

    real(ireals), dimension(:,:), intent(out), optional :: opt_swdirflx, opt_swuflx, opt_swdflx ! [nlay+1, ncol]
    real(ireals), dimension(:,:), intent(out), optional :: opt_swhr ! [nlay, ncol]
    real(ireals), dimension(:,:,:), intent(out),optional :: opt_tau_f, opt_w0_f, opt_g_f ! [nlay, ncol, ngptsw]

    real(rb),dimension(ncol_in,nlay_in) :: play, cldfr

    real(rb),dimension(nbndsw, ncol_in, nlay_in) :: taucld, ssacld, asmcld, fsfcld
    real(rb),dimension(ncol_in, nlay_in, nbndsw ) :: tauaer, ssaaer, asmaer
    real(rb),dimension(ncol_in, nlay_in, naerec ) :: ecaer

    real(rb),dimension(ncol_in) :: tsfc, asdir, aldir, asdif, aldif, coszen

    real(rb),dimension(ncol_in,nlay_in+1) :: swuflx,swdflx,swuflxc,swdflxc
    real(rb),dimension(ncol_in,nlay_in  ) :: swhr,swhrc

    integer(im) :: icol, ncol, nlay, k

    integer(im),parameter :: dyofyr=0,inflgsw=2,iceflgsw=3,liqflgsw=1
    real(rb)   ,parameter :: adjes=1, scon=1.36822e+03
    real(rb) :: solar_const
    integer(kind=im) :: icld=2         ! Cloud overlap method
    integer(kind=im) :: iaer=0         ! Aerosol option flag

    logical,save :: linit_rrtmg=.False.

    if(present(opt_solar_constant)) then
      solar_const = opt_solar_constant
    else
      solar_const = scon
    endif

    ! copy from TenStream to RRTM precision:
    ncol   = int(ncol_in, kind=im)
    nlay   = int(nlay_in, kind=im)

    ! Take average pressure and temperature as mean values for voxels --
    ! Todo: should we use log interpolation for pressure...?
    do icol=1,ncol
      play(icol,:) = .5_rb*(plev(icol,1:nlay)+plev(icol,2:nlay+1))

      tsfc(icol)   = tlev(icol,1)
    enddo

    taucld = 0; ssacld = 0; asmcld  = 0;
    fsfcld = 0;
    tauaer = 0; ssaaer  = 0; asmaer  = 0;
    ecaer  = 0; asdir   = albedo; aldir   = albedo;
    asdif  = albedo; aldif  = albedo; swdflxc = 0; swuflxc = 0;
    coszen = cos(deg2rad(theta0));

    if(present(opt_cldfr)) then
      cldfr = real(opt_cldfr, rb)
    else
      where(lwp.gt.0)
        cldfr = 1
      elsewhere
        cldfr = 0
      end where
    endif

    if(.not.linit_rrtmg) then
      call rrtmg_sw_ini(1006._rb)
      linit_rrtmg = .True.

      if(ldebug) then
        do k=nlay,1,-1
          print *,'rrtm_optprop_sw',k,'tlev',tlev(1,k),'tlay',tlay(1,k),'plev',plev(1,k),'play',play(1,k), &
            'lwp',lwp(1,k), 'reliq',reliq(1,k), 'iwp',iwp(1,k), 'reice',reice(1,k),&
            'h2o',h2ovmr(1,k), 'o3' , o3vmr(1,k), 'co2', co2vmr(1,k), 'ch4', ch4vmr(1,k),    &
            'n2o', n2ovmr(1,k), 'o2' , o2vmr(1,k)
        enddo
      endif
    endif

    if ( all([ present(opt_swdirflx), present(opt_swuflx), present(opt_swdflx), present(opt_swhr) ])) then
      call rrtmg_sw &
        (ncol, nlay, icld, iaer, real(play,rb), real(plev,rb), &
        real(tlay,rb), real(tlev,rb), tsfc,                    &
        real(h2ovmr,rb), real(o3vmr,rb), real(co2vmr,rb),      &
        real(ch4vmr,rb), real(n2ovmr,rb), real(o2vmr,rb),      &
        asdir, asdif, aldir, aldif, &
        coszen, adjes, dyofyr, solar_const, &
        inflgsw, iceflgsw, liqflgsw, cldfr, &
        taucld, ssacld, asmcld, fsfcld, &
        real(iwp, rb), real(lwp, rb), real(reice, rb), real(reliq, rb), &
        tauaer, ssaaer, asmaer, ecaer, &
        swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
        tau, w0, g, loptprop_only=.False., &
        tenstr_tau_f=opt_tau_f, tenstr_w_f=opt_w0_f, tenstr_g_f=opt_g_f)

      opt_swuflx = transpose(real(swuflx, ireals))
      opt_swdflx = transpose(real(swdflx, ireals))
      opt_swhr   = transpose(real(swhr, ireals))

      if(present(opt_solar_constant)) then
        call Edir_lambert_beer(&
          theta0, &
          tenstr_solsrc(:) / sum(tenstr_solsrc) * opt_solar_constant, &
          tau, opt_swdirflx)
      else
        call Edir_lambert_beer(&
          theta0, &
          tenstr_solsrc(:), &
          tau, opt_swdirflx)
      endif
      opt_swdflx = opt_swdflx - opt_swdirflx ! remove direct radiation part from edn

    else
      if ( any([ present(opt_swdirflx), present(opt_swuflx), present(opt_swdflx), present(opt_swhr) ])) &
        call CHKERR(1_mpiint, 'guess you wanted to call rrtmg and compute fluxes '// &
                              'but in this case you have to present all arguments!')

      call rrtmg_sw &
        (ncol, nlay, icld, iaer, real(play,rb), real(plev,rb), &
        real(tlay,rb), real(tlev,rb), tsfc,                    &
        real(h2ovmr,rb), real(o3vmr,rb), real(co2vmr,rb),      &
        real(ch4vmr,rb), real(n2ovmr,rb), real(o2vmr,rb),      &
        asdir, asdif, aldir, aldif, &
        coszen, adjes, dyofyr, solar_const, &
        inflgsw, iceflgsw, liqflgsw, cldfr, &
        taucld, ssacld, asmcld, fsfcld, &
        real(iwp, rb), real(lwp, rb), real(reice, rb), real(reliq, rb), &
        tauaer, ssaaer, asmaer, ecaer, &
        swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
        tau, w0, g, loptprop_only=.True., &
        tenstr_tau_f=opt_tau_f, tenstr_w_f=opt_w0_f, tenstr_g_f=opt_g_f)
    endif
  end subroutine

  ! Compute direct radiation from lambert beers law.
  ! This is good to split the edn flx which rrtmg returns into a direct component and a diffuse component
  ! ordering of tau and flx here start from surface and go up to TOA
  subroutine Edir_lambert_beer(theta0, E0, dtau, edir)
    real(ireals), intent(in) :: theta0     ! solar angle [deg]
    real(ireals), intent(in) :: E0(:)      ! solar incident irradiance at TOA dim(nbands)
    real(ireals), intent(in) :: dtau(:,:,:)  ! vertical optical depth dim(layers, ncol, nbands)
    real(ireals), intent(out):: edir(:,:)  ! direct solar irradiance dim(levels, ncol)
    real(ireals) :: mu, mu_inv
    real(ireals) :: E
    integer(iintegers) :: k, ib, i
    mu = cos(deg2rad(theta0))
    mu_inv = 1._ireals / mu
    edir = 0

    do ib=1,size(E0)
      do i=1,size(dtau,2)
        E = E0(ib) * mu
        edir(size(edir),i) = edir(size(edir),i) + E

        do k=size(dtau,1),1,-1
          E = E * exp(- dtau(k,i,ib) * mu_inv )
          edir(k,i) = edir(k,i) + E
        enddo
      enddo
    enddo
  end subroutine

end module
