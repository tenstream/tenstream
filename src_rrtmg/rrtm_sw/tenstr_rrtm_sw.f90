module m_tenstr_rrtm_sw
    use rrtmg_sw_init, only: rrtmg_sw_ini
    use parkind, only: im => kind_im, rb => kind_rb
    use rrsw_wvn
    use parrrsw, only: ngptsw, nbndsw,naerec,jpb1, jpb2
    use rrtmg_sw_rad, only: rrtmg_sw
    use rrtmg_sw_spcvrt, only: tenstr_solsrc      
    use iso_c_binding

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, myid, zero, one, i0, i1

    use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
        tenstream_get_result, C_one
    use m_tenstream_options, only: read_commandline_options

    use m_helper_functions, only: imp_reduce_sum

    implicit none

#include "petsc/finclude/petsc.h90"
    PetscErrorCode :: ierr

    logical :: linit=.False.

contains

    subroutine tenstr_rrtm_sw(comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir,edn,eup,abso)
        integer(iintegers), intent(in) :: comm
        integer(iintegers), intent(in) :: nlay, nxp, nyp

        real(ireals), intent(in) :: dx, dy, phi0, theta0, albedo
        
        real(ireals),intent(in) :: plev   (:,:,:) ! dim(nlay+1, nxp, nyp)

        real(ireals),intent(in) :: tlay   (:,:,:) ! dim(nlay  , nxp, nyp)
        real(ireals),intent(in) :: h2ovmr (:,:,:)
        real(ireals),intent(in) :: o3vmr  (:,:,:)
        real(ireals),intent(in) :: co2vmr (:,:,:)
        real(ireals),intent(in) :: ch4vmr (:,:,:)
        real(ireals),intent(in) :: n2ovmr (:,:,:)
        real(ireals),intent(in) :: o2vmr  (:,:,:)
        real(ireals),intent(in) :: lwc    (:,:,:)
        real(ireals),intent(in) :: reliq  (:,:,:)

        real(ireals) :: col_plev   (nxp*nyp, nlay+1)
        real(ireals) :: col_tlay   (nxp*nyp, nlay)
        real(ireals) :: col_h2ovmr (nxp*nyp, nlay)
        real(ireals) :: col_o3vmr  (nxp*nyp, nlay)
        real(ireals) :: col_co2vmr (nxp*nyp, nlay)
        real(ireals) :: col_ch4vmr (nxp*nyp, nlay)
        real(ireals) :: col_n2ovmr (nxp*nyp, nlay)
        real(ireals) :: col_o2vmr  (nxp*nyp, nlay)
        real(ireals) :: col_lwc    (nxp*nyp, nlay)
        real(ireals) :: col_reliq  (nxp*nyp, nlay)

        real(ireals) :: hhl(nlay+1, nxp, nyp)
        real(ireals) :: dz (nlay  , nxp, nyp)

        integer(iintegers) :: i, j, k, icol, ib
        integer(iintegers) :: is,ie,js,je

        integer(iintegers) :: nbands
        real(ireals),allocatable,dimension(:)   :: band_lbound,band_ubound,weights ! [nbands]
        real(ireals),allocatable,dimension(:, :, :) :: col_tau, col_w0, col_g      ! [ncol, nlyr, nbands]
        real(ireals),allocatable, dimension(:,:,:,:) :: kabs,ksca,g                ! [nlyr, local_nx, local_ny, nbands]
        real(ireals),allocatable, dimension(:,:,:), intent(out) :: edir,edn,eup,abso          ! [nlyr(+1), local_nx, local_ny ]
        real(ireals),allocatable, dimension(:,:,:) :: spec_edir,spec_edn,spec_eup,spec_abso ! [nlyr(+1), local_nx, local_ny ]

        do j=1,nyp
            do i=1,nxp
                call hydrostat_lev(plev(:,i,j),tlay(:,i,j), hhl(:,i,j))
            enddo
        enddo
        dz = hhl(2:nlay+1,:,:)-hhl(1:nlay,:,:)

        call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, albedo, dz1d=dz(:,1,1))
        if(myid.eq.0) print *,'dz1d',dz(:,1,1)

        is = C_one%xs +1
        ie = C_one%xe +1
        js = C_one%ys +1
        je = C_one%ye +1

        call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, albedo, dz3d=dz(:,is:ie,js:je))

        print *,'domain has sizes:',is,ie,js,je

        icol=0
        do j=js,je
            do i=is,ie
                icol = icol+1
                col_plev  (icol, :) = plev  (:, i, j)  
                col_tlay  (icol, :) = tlay  (:, i, j)
                col_h2ovmr(icol, :) = h2ovmr(:, i, j)
                col_o3vmr (icol, :) = o3vmr (:, i, j)
                col_co2vmr(icol, :) = co2vmr(:, i, j)
                col_ch4vmr(icol, :) = ch4vmr(:, i, j)
                col_n2ovmr(icol, :) = n2ovmr(:, i, j)
                col_o2vmr (icol, :) = o2vmr (:, i, j)
                col_lwc   (icol, :) = lwc   (:, i, j)*dz(:,i,j)
                col_reliq (icol, :) = reliq (:, i, j)
            enddo
        enddo

        call optprop_rrtm_sw(icol, nlay, &
            col_plev, col_tlay, &
            col_h2ovmr, col_o3vmr , col_co2vmr,&
            col_ch4vmr, col_n2ovmr, col_o2vmr ,&
            col_lwc, col_reliq,        &
            nbands, band_lbound, band_ubound, weights, &
            col_tau, col_w0, col_g)


        allocate(kabs(nlay, is:ie, js:je, ngptsw))
        allocate(ksca(nlay, is:ie, js:je, ngptsw))
        allocate(g   (nlay, is:ie, js:je, ngptsw))

        allocate(edir(nlay+1, is:ie, js:je), source=zero)
        allocate(edn (nlay+1, is:ie, js:je), source=zero)
        allocate(eup (nlay+1, is:ie, js:je), source=zero)
        allocate(abso(nlay  , is:ie, js:je), source=zero)

        allocate(spec_edir(nlay+1, is:ie, js:je))
        allocate(spec_edn (nlay+1, is:ie, js:je))
        allocate(spec_eup (nlay+1, is:ie, js:je))
        allocate(spec_abso(nlay  , is:ie, js:je))

        icol=0
        do j=js,je
            do i=is,ie
                icol = icol+1
                ! copy from number columns of rrtm interface back onto regular grid
                kabs(:,i,j,:) = max(zero, col_tau(icol,:,:) * (one-col_w0 (icol,:,:))) 
                ksca(:,i,j,:) = max(zero, col_tau(icol,:,:) *      col_w0 (icol,:,:) ) 
                g   (:,i,j,:) = col_g  (icol,:,:)

                ! reverse output from rrtm to start at TOA again
                kabs(:,i,j,:) = kabs(nlay:1:-1,i,j,:)
                ksca(:,i,j,:) = ksca(nlay:1:-1,i,j,:)
                g   (:,i,j,:) = g   (nlay:1:-1,i,j,:)

                ! divide by thickness to convert from tau to coefficients per meter
                do ib=1, nbands
                    kabs(:,i,j,ib) = kabs(:,i,j,ib) / dz(:,i,j)
                    ksca(:,i,j,ib) = ksca(:,i,j,ib) / dz(:,i,j)
                enddo
            enddo
        enddo

        do ib=1, nbands
            call set_optical_properties(kabs(:,:,:,ib), ksca(:,:,:,ib), g(:,:,:,ib))
            call solve_tenstream(weights(ib))
            call tenstream_get_result(spec_edir, spec_edn, spec_eup, spec_abso)
            edir = edir + spec_edir
            edn  = edn  + spec_edn 
            eup  = eup  + spec_eup 
        enddo
    end subroutine

    subroutine optprop_rrtm_sw(ncol_in, nlay_in, plev_in, tlay_in, h2ovmr_in, o3vmr_in, co2vmr_in, ch4vmr_in, n2ovmr_in, o2vmr_in, lwp_in, reliq_in, nbands, band_lbound, band_ubound, weights, tau, w0, g)
        integer(im),parameter :: dyofyr=0,inflgsw=2,iceflgsw=3,liqflgsw=1
        real(rb)   ,parameter :: adjes=1, scon=1.36822e+03

        ! RRTM needs the arrays to start at the surface
        ! Input is however given starting from top
        ! Copy input variables from tenstream precision to rrtm precision
        integer(iintegers),intent(in) :: ncol_in, nlay_in

        real(ireals),dimension(:,:),intent(in) :: plev_in, tlay_in, h2ovmr_in, o3vmr_in, co2vmr_in, ch4vmr_in, n2ovmr_in, o2vmr_in ! [ncol_in, nlay]
        real(ireals),dimension(:,:),intent(in) :: lwp_in, reliq_in ! [ncol_in, nlay]

        ! and use those without '_in' as before
        integer(im) :: ncol, nlay


        real(rb),dimension(ncol_in,nlay_in+1) :: plev 
        real(rb),dimension(ncol_in,nlay_in)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr 
        real(rb),dimension(ncol_in,nlay_in)   :: lwp, reliq 

        integer(iintegers),intent(out) :: nbands
        real(ireals), allocatable, dimension(:),intent(out)   :: band_lbound, band_ubound, weights ! [nbands]
        real(ireals), allocatable, dimension(:,:,:),intent(out) :: tau, w0, g ! [ncol, nlay, nbands]

        real(rb),dimension(ncol_in,nlay_in+1) :: tlev

        real(rb),dimension(ncol_in,nlay_in) :: play, cldfr, cicewp, reice

        real(rb),dimension(nbndsw, ncol_in,nlay_in) :: taucld, ssacld, asmcld, fsfcld
        real(rb),dimension(ncol_in, nlay_in, nbndsw ) :: tauaer, ssaaer, asmaer
        real(rb),dimension(ncol_in, nlay_in, naerec ) :: ecaer

        real(rb),dimension(ncol_in) :: tsfc, asdir, aldir, asdif, aldif, coszen

        real(rb),dimension(ncol_in,nlay_in+1) :: swuflx,swdflx,swuflxc,swdflxc
        real(rb),dimension(ncol_in,nlay_in  ) :: swhr,swhrc

        integer(im) :: iv,ig,k,icol,iw

        integer(kind=im) :: icld=2         ! Cloud overlap method
        integer(kind=im) :: iaer=0         ! Aerosol option flag

        if(.not.allocated(band_lbound)) allocate(band_lbound(ngptsw))
        if(.not.allocated(band_ubound)) allocate(band_ubound(ngptsw))
        if(.not.allocated(weights    )) allocate(weights    (ngptsw))

        if(.not.allocated(tau)) allocate(tau(ncol_in, nlay_in, ngptsw))
        if(.not.allocated(w0 )) allocate(w0 (ncol_in, nlay_in, ngptsw))
        if(.not.allocated(g  )) allocate(g  (ncol_in, nlay_in, ngptsw))

        ! copy to correct precision:
        ncol   = ncol_in
        nlay   = nlay_in
        plev   = rev(plev_in)
        tlay   = rev(tlay_in)
        h2ovmr = rev(h2ovmr_in)
        o3vmr  = rev(o3vmr_in)
        co2vmr = rev(co2vmr_in)
        ch4vmr = rev(ch4vmr_in)
        n2ovmr = rev(n2ovmr_in)
        o2vmr  = rev(o2vmr_in)
        lwp    = rev(lwp_in)
        reliq  = rev(reliq_in)

        do icol=1,ncol
            play(icol,:)      = .5_rb*(plev(icol,1:nlay)+plev(icol,2:nlay+1))

            tlev(icol,nlay+1) = tlay(icol,nlay)
            tsfc(icol)        = tlev(icol,nlay+1)
            tlev(icol,1)      = tlay(icol,1)

            do k = 1,nlay-1
                tlev(icol,k+1)  = .5_rb * (tlay(icol,k+1)+tlay(icol,k) )
            enddo
        enddo

        taucld = 0; ssacld = 0; asmcld  = 0;
        fsfcld = 0; cicewp = 0; reice   = 0;
        tauaer = 0; ssaaer  = 0; asmaer  = 0;
        ecaer  = 0; coszen = 1; asdir   = 0; aldir   = 0;
        asdif  = 0; aldif  = 0; swdflxc = 0; swuflxc = 0;

        where ( lwp.gt.0 )
            cldfr = 1
        elsewhere
            cldfr = 0
        endwhere

        if(.not.linit) then
            call rrtmg_sw_ini(1006._rb)
            linit = .True.

!            if(myid.eq.0) then
!                do k=1,nlay
!                    print *,k,'tlev',tlev(1,k),'tlay',tlay(1,k),'plev',plev(1,k),'play',play(1,k),'lwp',lwp(1,k)
!                enddo
!            endif
        endif

        call rrtmg_sw &
            (ncol    ,nlay    ,icld    ,iaer    , &
            play    ,plev    ,tlay    ,tlev    ,tsfc    , &
            h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
            asdir   ,asdif   ,aldir   ,aldif   , &
            coszen  ,adjes   ,dyofyr  ,scon    , &
            inflgsw ,iceflgsw,liqflgsw,cldfr   , &
            taucld  ,ssacld  ,asmcld  ,fsfcld  , &
            cicewp  ,lwp  ,reice   ,reliq      , &
            tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
            swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
            tau, w0, g)

        nbands = ngptsw

        iw=0
        do iv=1,nbndsw
            do ig=1,ngc(iv)
                iw = iw+1
                band_lbound(iw) = wavenum1(iv+jpb1-1)
                band_ubound(iw) = wavenum2(iv+jpb1-1)
                weights(iw)     = tenstr_solsrc(iw)
            enddo
        enddo
    end subroutine

    subroutine hydrostat_lev(plev,tlay, hhl)
        real(ireals),intent(in) :: plev(:),tlay(:)
        real(ireals),intent(out) :: hhl(size(plev))
        real(ireals) :: dp,dz,rho
        integer(im) :: k
        hhl(1) = zero
        do k=1,size(tlay)
            dp  = abs( plev(k)-plev(k+1) ) 
            rho = ( plev(k)+dp/2._ireals  ) / 287.058_ireals / tlay(k)
            dz  = dp / rho / 9.8065_ireals
            hhl(k+1) = hhl(k) + dz
        enddo
    end subroutine

    function rev(inp) ! reverse second dimension
        real(rb),intent(in) :: inp(:,:)
        real(rb) :: rev(size(inp,1),size(inp,2))
        rev = inp(:,ubound(inp,2):lbound(inp,2):-1)
    end function

end module

