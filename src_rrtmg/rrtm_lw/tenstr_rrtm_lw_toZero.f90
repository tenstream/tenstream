! Routines to call tenstream with optical properties from RRTM -- to be called from rank zero
! for functions callable from unique rank, please refer to m_tenstr_rrtm_lw

module m_tenstr_rrtm_lw_toZero

      use m_tenstr_rrtmg_lw_init, only: rrtmg_lw_ini
      use m_tenstr_parkind_lw, only: im => kind_im, rb => kind_rb
      use m_tenstr_rrlw_wvn, only : ngb, wavenum1, wavenum2
      use m_tenstr_parrrtm, only: ngptlw, nbndlw
      use m_tenstr_rrtmg_lw_rad, only: rrtmg_lw

      use m_data_parameters, only : init_mpi_data_parameters, &
        iintegers, ireals, myid, zero, one, i0, i1,           &
        mpiint, pi, default_str_len

      use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
          tenstream_get_result, tenstream_get_result_toZero, C_one

      use m_netcdfIO, only : ncwrite

    implicit none

    private
    public :: tenstream_rrtm_lw_toZero

    logical :: linit=.False.

    interface
      real function PLKINT(WVLLO, WVLHI, T)
        real :: WVLLO, WVLHI, T
      end function 
    end interface


contains

    subroutine tenstream_rrtm_lw_toZero(comm, nlay, nxp, nyp, dx, dy, phi0, theta0, albedo, plev, tlay, tsrfc, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, lwc, reliq, edir,edn,eup,abso, icollapse)
        integer(mpiint), intent(in) :: comm
        integer(iintegers), intent(in) :: nlay, nxp, nyp

        real(ireals), intent(in) :: dx, dy, phi0, theta0, albedo

        real(ireals),intent(in) :: plev   (:,:,:) ! dim(nlay+1, nxp, nyp)

        real(ireals),intent(in) :: tlay   (:,:,:) ! all have
        real(ireals),intent(in) :: h2ovmr (:,:,:) ! dim(nlay , nxp, nyp)
        real(ireals),intent(in) :: o3vmr  (:,:,:) !
        real(ireals),intent(in) :: co2vmr (:,:,:) !
        real(ireals),intent(in) :: ch4vmr (:,:,:) !
        real(ireals),intent(in) :: n2ovmr (:,:,:) !
        real(ireals),intent(in) :: o2vmr  (:,:,:) !
        real(ireals),intent(in) :: lwc    (:,:,:) !
        real(ireals),intent(in) :: reliq  (:,:,:) !
        real(ireals),intent(in) :: tsrfc  (:,:)   ! dim(nxp, nyp)

        integer(iintegers),optional :: icollapse

        real(ireals),allocatable :: col_plev   (:,:)
        real(ireals),allocatable :: col_tlay   (:,:)
        real(ireals),allocatable :: col_h2ovmr (:,:)
        real(ireals),allocatable :: col_o3vmr  (:,:)
        real(ireals),allocatable :: col_co2vmr (:,:)
        real(ireals),allocatable :: col_ch4vmr (:,:)
        real(ireals),allocatable :: col_n2ovmr (:,:)
        real(ireals),allocatable :: col_o2vmr  (:,:)
        real(ireals),allocatable :: col_lwc    (:,:)
        real(ireals),allocatable :: col_reliq  (:,:)

        real(ireals) :: hhl(nlay+1, nxp, nyp)
        real(ireals) :: dz (nlay  , nxp, nyp)

        integer(iintegers) :: i, j, k, icol, ib, ig
        integer(iintegers) :: is,ie,js,je
        real(ireals) :: global_maxheight

        real(ireals),allocatable, dimension(:, :, :) :: col_tau, col_Bfrac                    ! [ncol, nlyr, ngptlw]
        real(ireals),allocatable, dimension(:,:,:)   :: ksca,g                                ! [nlyr, local_nx, local_ny, ngptlw]
        real(ireals),allocatable, dimension(:,:,:,:) :: kabs,Bfrac                            ! [nlyr, local_nx, local_ny, ngptlw]
        real(ireals),allocatable, dimension(:,:,:,:) :: Blev                                  ! [nlyr+1, local_nx, local_ny, nbndlw]
        real(ireals),allocatable, dimension(:,:,:)   :: spec_edir,spec_edn,spec_eup,spec_abso ! [nlyr(+1), local_nx, local_ny ]

        real(ireals),allocatable, dimension(:,:,:), intent(out) :: edir,edn,eup,abso          ! [nlyr(+1), local_nx, local_ny ]

        character(default_str_len) :: output_path(2) ! [ filename, varname ]

        do j=1,nyp
            do i=1,nxp
                call hydrostat_lev(plev(:,i,j),tlay(:,i,j), zero, hhl(:,i,j), dz(:,i,j))
            enddo
        enddo
        global_maxheight = maxval(hhl)
        do j=1,nyp
            do i=1,nxp
                hhl(:, i, j) = hhl(:, i, j) + global_maxheight - hhl(1, i, j)
            enddo
        enddo
        output_path(1) = 'output.nc'
        if(myid.eq.0) then
            output_path(2) = 'dz3d' ; call ncwrite(output_path, dz, i)
            output_path(2) = 'hhl'  ; call ncwrite(output_path, hhl, i)
            output_path(2) = 'hsrfc'; call ncwrite(output_path, hhl(ubound(hhl,1),:,:), i)
        endif

        if(present(icollapse)) then 
            call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, dz1d=dz(:,1,1), collapseindex=icollapse)
            is = C_one%xs +1; ie = C_one%xe +1; js = C_one%ys +1; je = C_one%ye +1
            call destroy_tenstream(.True.)
            call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, dz3d=dz(:,is:ie,js:je), collapseindex=icollapse)
        else
            call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, dz1d=dz(:,1,1))
            is = C_one%xs +1; ie = C_one%xe +1; js = C_one%ys +1; je = C_one%ye +1
            call destroy_tenstream(.True.)
            call init_tenstream(comm, nlay, nxp, nyp, dx,dy,phi0, theta0, dz3d=dz(:,is:ie,js:je))
        endif

        allocate(col_plev   (C_one%xm*C_one%ym, nlay+1))
        allocate(col_tlay   (C_one%xm*C_one%ym, nlay))
        allocate(col_h2ovmr (C_one%xm*C_one%ym, nlay))
        allocate(col_o3vmr  (C_one%xm*C_one%ym, nlay))
        allocate(col_co2vmr (C_one%xm*C_one%ym, nlay))
        allocate(col_ch4vmr (C_one%xm*C_one%ym, nlay))
        allocate(col_n2ovmr (C_one%xm*C_one%ym, nlay))
        allocate(col_o2vmr  (C_one%xm*C_one%ym, nlay))
        allocate(col_lwc    (C_one%xm*C_one%ym, nlay))
        allocate(col_reliq  (C_one%xm*C_one%ym, nlay))

        allocate(col_tau  (C_one%xm*C_one%ym, nlay, ngptlw))
        allocate(col_Bfrac(C_one%xm*C_one%ym, nlay, ngptlw))

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

        call optprop_rrtm_lw(icol, nlay, albedo,&
            col_plev, col_tlay, &
            col_h2ovmr, col_o3vmr , col_co2vmr,&
            col_ch4vmr, col_n2ovmr, col_o2vmr ,&
            col_lwc, col_reliq,       &
            col_tau, col_Bfrac)


        allocate(kabs (nlay,  is:ie, js:je, ngptlw))
        allocate(ksca (nlay,  is:ie, js:je))
        allocate(g    (nlay,  is:ie, js:je))
        allocate(Bfrac(nlay+1,is:ie, js:je, ngptlw))

        ! rrtmg_lw does not support thermal scattering... set to zero
        ksca = zero
        g    = zero

        icol=0
        do j=js,je
            do i=is,ie
                icol = icol+1
                ! copy from number columns of rrtm interface back onto regular grid
                kabs(:,i,j,:) = max(zero, col_tau(icol,:,:))
                ! reverse output from rrtm to start at TOA again
                kabs(:,i,j,:) = kabs(nlay:1:-1,i,j,:)

                ! divide by thickness to convert from tau to coefficients per meter
                do ib=1, ngptlw
                    kabs(:,i,j,ib) = kabs(:,i,j,ib) / dz(:,i,j)
                enddo

                Bfrac(2:nlay+1,i,j,:) = col_Bfrac(icol,:,:)
                Bfrac(1,i,j,:) = Bfrac(2,i,j,:) ! surface weights are the same as lowest layer
                Bfrac(:,i,j,:) = Bfrac(nlay+1:1:-1,i,j,:)
            enddo
        enddo

        ! Free up some intermediate memory
        deallocate(col_plev  )
        deallocate(col_tlay  )
        deallocate(col_h2ovmr)
        deallocate(col_o3vmr )
        deallocate(col_co2vmr)
        deallocate(col_ch4vmr)
        deallocate(col_n2ovmr)
        deallocate(col_o2vmr )
        deallocate(col_lwc   )
        deallocate(col_reliq )
        deallocate(col_tau   )
        deallocate(col_Bfrac )

        ! Allocate space for results -- for integrated values and for temporary spectral integration...
        if(myid.eq.0) then
            allocate(edn (nlay+1, C_one%glob_xm, C_one%glob_ym ), source=zero)
            allocate(eup (nlay+1, C_one%glob_xm, C_one%glob_ym ), source=zero)
            allocate(abso(nlay  , C_one%glob_xm, C_one%glob_ym ), source=zero)

            allocate(spec_edn (nlay+1, C_one%glob_xm, C_one%glob_ym))
            allocate(spec_eup (nlay+1, C_one%glob_xm, C_one%glob_ym))
            allocate(spec_abso(nlay  , C_one%glob_xm, C_one%glob_ym))
        endif

        ! Compute source term(planck function)
        allocate(Blev (nlay+1,is:ie, js:je, nbndlw))
        do ib=1,nbndlw
          do j=js,je
            do i=is,ie
              do k=1,nlay
                Blev(k,i,j,ib) = plkint(real(wavenum1(ib)), real(wavenum2(ib)), real(tlay(k,i,j)))
              enddo
              Blev(nlay+1,i,j,ib) = plkint(real(wavenum1(ib)), real(wavenum2(ib)), real(tsrfc(i,j)))
            enddo
          enddo
        enddo

        ! Loop over spectral intervals and call solver
        do ib=1,ngptlw 
            call set_optical_properties(albedo, kabs(:,:,:,ib), ksca(:,:,:), g(:,:,:), Blev(:,:,:,ngb(ib))*Bfrac(:,:,:,ib))
            call solve_tenstream(zero)
            call tenstream_get_result_toZero(spec_edir, spec_edn, spec_eup, spec_abso)

            if(myid.eq.0) then
                edn  = edn  + spec_edn 
                eup  = eup  + spec_eup 
                abso = abso + spec_abso
            endif
        enddo

        ! Tidy up the solver
        call destroy_tenstream(lfinalizepetsc=.True.)
    end subroutine

    subroutine optprop_rrtm_lw(ncol_in, nlay_in, albedo, plev_in, tlay_in, h2ovmr_in, o3vmr_in, co2vmr_in, ch4vmr_in, n2ovmr_in, o2vmr_in, lwp_in, reliq_in, tau, Bfrac)
        ! RRTM needs the arrays to start at the surface
        ! Input is however given starting from top
        ! Copy input variables from tenstream precision to rrtm precision

        integer(iintegers),intent(in)          :: ncol_in, nlay_in
        real(ireals), intent(in) :: albedo
        real(ireals),dimension(:,:),intent(in) :: plev_in, tlay_in, h2ovmr_in, o3vmr_in, co2vmr_in, ch4vmr_in, n2ovmr_in, o2vmr_in ! [ncol_in, nlay]
        real(ireals),dimension(:,:),intent(in) :: lwp_in, reliq_in ! [ncol_in, nlay]

        ! and use those without '_in' as before

        real(rb),dimension(ncol_in,nlay_in+1) :: plev 
        real(rb),dimension(ncol_in,nlay_in)   :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr 
        real(rb),dimension(ncol_in,nlay_in)   :: lwp, reliq 

        real(ireals), dimension(:,:,:), intent(out) :: tau, Bfrac ! [ncol, nlay, ngptlw]

        real(rb),dimension(ncol_in,nlay_in+1) :: tlev
        real(rb),dimension(ncol_in,nlay_in) :: play, cldfr, cicewp, reice

        real(rb),dimension(nbndlw, ncol_in, nlay_in) :: taucld
        real(rb),dimension(ncol_in, nlay_in, nbndlw ) :: tauaer
        real(rb),dimension(ncol_in, nbndlw ) :: emis

        real(rb),dimension(ncol_in, nlay_in)   :: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr

        real(rb),dimension(ncol_in) :: tsfc 

        real(rb),dimension(ncol_in,nlay_in+1) :: lwuflx,lwdflx,lwuflxc,lwdflxc
        real(rb),dimension(ncol_in,nlay_in  ) :: lwhr,lwhrc

        integer(im) :: iv,ig,k,icol,iw
        integer(im) :: ncol, nlay

        integer(im),parameter :: inflglw=2,iceflglw=3,liqflglw=1
        integer(kind=im) :: icld=2         ! Cloud overlap method
        integer(kind=im) :: iaer=0         ! Aerosol option flag
        integer(kind=im) :: idrv=0         ! Flag for calculation of dFdT

        ! copy from TenStream to RRTM precision:
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

        ! Take average pressure and temperature as mean values for voxels --
        ! should probably use log interpolation for pressure...
        do icol=1,ncol
            play(icol,:)      = .5_rb*(plev(icol,1:nlay)+plev(icol,2:nlay+1))

            tlev(icol,nlay+1) = tlay(icol,nlay)
            tsfc(icol)        = tlev(icol,nlay+1)
            tlev(icol,1)      = tlay(icol,1)

            do k = 1,nlay-1
                tlev(icol,k+1)  = .5_rb * (tlay(icol,k+1)+tlay(icol,k) )
            enddo
        enddo

        taucld   = 0; cicewp   = 0; reice    = 0;
        tauaer   = 0; lwdflxc  = 0; lwuflxc  = 0;
        cfc11vmr = 0; cfc12vmr = 0; cfc22vmr = 0; ccl4vmr = 0;

        emis = one - albedo

        where ( lwp.gt.0 )
            cldfr = 1
        elsewhere
            cldfr = 0
        endwhere

        if(.not.linit) then
            call rrtmg_lw_ini(1006._rb)
            linit = .True.

!            if(myid.eq.0) then
!                do k=1,nlay
!                    print *,k,'tlev',tlev(1,k),'tlay',tlay(1,k),'plev',plev(1,k),'play',play(1,k),'lwp',lwp(1,k)
!                enddo
!            endif
        endif

        call rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv   , &
            play    ,plev    ,tlay    ,tlev    ,tsfc    , &
            h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
            cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
            inflglw, iceflglw, liqflglw, cldfr, &
            taucld , cicewp  , lwp  ,reice   ,reliq      , &
            tauaer , &
            lwuflx , lwdflx  ,lwhr    ,lwuflxc ,lwdflxc ,lwhrc, &
            tau, Bfrac)
    end subroutine

    subroutine hydrostat_lev(plev,tlay, hsrfc, hhl, dz)
        ! Integrate vertical height profile hydrostatically.
        real(ireals),intent(in) :: plev(:),tlay(:)
        real(ireals),intent(in) :: hsrfc
        real(ireals),intent(out) :: hhl(size(plev))
        real(ireals),intent(out) :: dz(size(tlay))
        real(ireals) :: dp,rho
        integer(im) :: k
        hhl(size(plev)) = hsrfc
        do k=size(tlay),1,-1
            dp  = abs( plev(k)-plev(k+1) ) 
            rho = ( plev(k)+dp/2._ireals  ) / 287.058_ireals / tlay(k)
            dz(k) = dp / rho / 9.8065_ireals
            hhl(k) = hhl(k+1) + dz(k)
        enddo
        if(any(dz.le.zero)) then
            print *,'plev',plev
            print *,'tlay',tlay
            print *,'dz',dz
            print *,'hhl',hhl
            stop 'error in dz'
        endif
    end subroutine

    function rev(inp) ! reverse second dimension
        real(ireals),intent(in) :: inp(:,:)
        real(rb) :: rev(size(inp,1),size(inp,2))
        rev = inp(:,ubound(inp,2):lbound(inp,2):-1)
    end function

end module

