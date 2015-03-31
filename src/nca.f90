!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2015, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------
module m_ts_nca
  use m_data_parameters, only : ireals,iintegers
  implicit none
  !global vars

contains
  subroutine ts_nca (nlyr, Nx, Ny, dz, B, kabs_3d, hr_nca_3d_tmp, dx, dy, Ldn, Lup)

    integer(iintegers), intent(in) :: Nx, Ny, nlyr
    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in),dimension(nlyr  ,0:Nx+1,0:Ny+1) :: dz, kabs_3d
    real(ireals), intent(in),dimension(nlyr+1,0:Nx+1,0:Ny+1) :: B,Ldn,Lup

    real(ireals), intent(out) :: hr_nca_3d_tmp(nlyr,Nx,Ny)

    ! ############## Definition of variables ##################

    integer(iintegers)        :: ilyr = 0         ! counter z levels
    integer(iintegers)        :: ixx = 0          ! counter x grid boxes
    integer(iintegers)        :: iyy = 0          ! counter y grid boxes
    integer(iintegers)        :: iface = 0        ! counter side faces/directions of grid box

    real(ireals)           :: Absup = 0        ! Upwelling, Absorption, lower/upper face
    real(ireals)           :: Absdn = 0        ! Downwelling Absorption, lower/upper face
    real(ireals)           :: Absup1 = 0       ! Upwelling, Absorption, face 1
    real(ireals)           :: Absup2 = 0       ! Upwelling, Absorption, face 2
    real(ireals)           :: Absup3 = 0       ! Upwelling, Absorption, face 3
    real(ireals)           :: Absup4 = 0       ! Upwelling, Absorption, face 4
    real(ireals)           :: Absdn1 = 0       ! Downwelling Absorption, face 1  
    real(ireals)           :: Absdn2 = 0       ! Downwelling Absorption, face 2
    real(ireals)           :: Absdn3 = 0       ! Downwelling Absorption, face 3
    real(ireals)           :: Absdn4 = 0       ! Downwelling Absorption, face 4

    real(ireals)           :: Emup = 0         ! Upwelling, Emission, lower/upper face
    real(ireals)           :: Emdn = 0         ! Downwelling Emission, lower/upper face
    real(ireals)           :: Emup1 = 0        ! Upwelling, Emission, face 1
    real(ireals)           :: Emup2 = 0        ! Upwelling, Emission, face 2
    real(ireals)           :: Emup3 = 0        ! Upwelling, Emission, face 3
    real(ireals)           :: Emup4 = 0        ! Upwelling, Emission, face 4
    real(ireals)           :: Emdn1 = 0        ! Downwelling Emission, face 1  
    real(ireals)           :: Emdn2 = 0        ! Downwelling Emission, face 2  
    real(ireals)           :: Emdn3 = 0        ! Downwelling Emission, face 3  
    real(ireals)           :: Emdn4 = 0        ! Downwelling Emission, face 4  
   
    real(ireals)           :: HR_up = 0        ! Downwelling Emission, face 1  
    real(ireals)           :: HR_dn = 0        ! Downwelling Emission, face 2  
    real(ireals)           :: HR_up_s = 0      ! Downwelling Emission, face 3  
    real(ireals)           :: HR_dn_s = 0      ! Downwelling Emission, face 4  

    real(ireals)           :: mu = 0           ! zenith angle

    real(ireals)           :: areaweight = 0   ! area weight for non-cubic grid boxes

    real(ireals)           :: ax = 0           ! integration boarder for side contributions
    real(ireals)           :: bx = 0           ! integration boarder for side contributions 
    real(ireals)           :: cx = 0           ! integration boarder for side contributions

    real(ireals)           :: az = 0           ! integration boarder for side contributions
    real(ireals)           :: bz = 0           ! integration boarder for side contributions 
    real(ireals)           :: cz = 0           ! integration boarder for side contributions 

    real(ireals)           :: pi=3.141592653589793
 
    real(ireals)           :: factor = 0
    real(ireals)           :: l = 0.0
    real(ireals)           :: Trans = 0.0

    real(ireals)           :: L_up_3d(nlyr+1,0:Nx+1,0:Ny+1)
    real(ireals)           :: L_dn_3d(nlyr+1,0:Nx+1,0:Ny+1)
    real(ireals)           :: B1 = 0
    real(ireals)           :: B2 = 0
    real(ireals)           :: B1_1 = 0
    real(ireals)           :: B2_1 = 0

    integer(iintegers) :: is,ie,js,je

    is = 1
    ie = Nx
    js = 1
    je = Ny

    ! ### get radiance from flux (back converted at the end)
    L_dn_3d=Ldn/pi
    L_up_3d=Lup/pi

    ! ################# start 3d calculation ####################
    ! ###########################################################
    do ilyr=1,nlyr    !  loop over all height levels 
       do ixx=is,ie           ! loop over all x gridboxes 
          do iyy=js,je         !  loop over all y gridboxes  


           !  print*,ilyr, iyy, ixx, dz(ilyr,ixx,iyy), dx, dy
             ! set and reset boundary conditions 
             Emdn=0.
             Emdn1=0.
             Emdn2=0.
             Emdn3=0.
             Emdn4=0.
             Absdn=0.
             Absdn1=0.
             Absdn2=0.
             Absdn3=0.
             Absdn4=0.
             Emup=0.
             Emup1=0.
             Emup2=0.
             Emup3=0.
             Emup4=0.
             Absup=0.
             Absup1=0.
             Absup2=0.
             Absup3=0.
             Absup4=0.
             areaweight=1.
    
             do iface=1,5     ! loop over the faces of the grid box! 

                mu=cos(45.*pi/180.) 

                if(mu.lt.0) then
                   mu=mu*(-1.)
                   print *,' Some error occured, mu should not be negative. Resetted to positive value, but result might be wrong! Please check!'
                endif

                
                ax=0
                !ay=0
                az=0
                cx=dx
                !  cy=dy
                cz=dz(ilyr,ixx,iyy)

                ! ## set parameters for integration 
                bx=cx-(dz(ilyr,ixx,iyy)/mu)*sqrt(1.-mu*mu)
              !  by=cy-(dz(ilyr,ixx,iyy)/mu)*sqrt(1.-mu*mu)
                bz=cz-(dx/sqrt(1.-mu*mu))*mu 


                if (bx.lt.ax) then
                   bx=ax
                endif
             !   if (bz.lt.az) then
             !      by=ay
             !   endif
                if (bz.lt.az) then
                   bz=az
                endif

                B2=B(ilyr+1,ixx,iyy)
                B1=B(ilyr,ixx,iyy)
                   
                ! if lower face of gridbox 
                if (iface .eq. 1) then
                   Trans=0
                   
                   if(ilyr.gt.1) then
                      B2_1=B(ilyr,ixx,iyy)
                      B1_1=B(ilyr-1,ixx,iyy)
                      Trans = integrate_flux((B1_1+B2_1)/2., L_dn_3d(ilyr-1,modulo(ixx,Nx)+1,iyy), &
                           kabs_3d(ilyr-1,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr-1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr-1,ixx,iyy)), dy,  mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_dn_3d(ilyr-1,modulo(ixx-2,Nx)+1,iyy),&
                           kabs_3d(ilyr-1,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr-1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr-1,ixx,iyy)), dy, mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_dn_3d(ilyr-1,ixx,modulo(iyy-2,Ny)+1),&
                           kabs_3d(ilyr-1,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr-1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr-1,ixx,iyy)), dx, mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_dn_3d(ilyr-1,ixx,modulo(iyy,Ny)+1),&
                           kabs_3d(ilyr-1,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr-1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr-1,ixx,iyy)), dx, mu)            
                      l = (Trans+L_dn_3d(ilyr,ixx,iyy))/5.
                   else
                      l = L_dn_3d(ilyr,ixx,iyy)               
                   end if
   
            !       if(l.gt.B(ilyr+1,ixx,iyy)) print*, l, B(ilyr+1, ixx,iyy)
                   Absdn =  l*integrate_emis(kabs_3d(ilyr,ixx,iyy),&
                        ax, bx, cx, dx, abs(dz(ilyr,ixx,iyy)), mu)
                   Emdn = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, abs(dz(ilyr,ixx,iyy)), mu)
                endif


                ! #### go to left ####
                if (iface .eq. 2) then

                   areaweight=dz(ilyr,ixx,iyy)/dy

                   Absdn1 = integrate_abs ((B1+B2)/2., L_dn_3d(ilyr,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight
                   Emdn1 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight
                endif ! iface .eq. 2 


                ! #### go to right #### 
                if (iface .eq. 3)  then

                   areaweight=dz(ilyr,ixx,iyy)/dy

                   Absdn2 = integrate_abs((B1+B2)/2.,  L_dn_3d(ilyr,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dy,  mu)*areaweight
                   Emdn2 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight

                endif ! end face 2 


                ! #### go to front #### 
                if (iface .eq. 4) then

                   areaweight=dz(ilyr,ixx,iyy)/dx

                   Absdn3 = integrate_abs((B1+B2)/2.,  L_dn_3d(ilyr,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                   Emdn3 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight    

                endif ! ifacd .eq. 4


                ! ##### go to back #### 
                if (iface .eq. 5) then

                   areaweight=dz(ilyr,ixx,iyy)/dx

                   Absdn4 = integrate_abs((B1+B2)/2.,  L_dn_3d(ilyr,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                   Emdn4 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                endif ! iface .eq. 5


              ! ## upwelling           

                ! if lower face of gridbox 
                if (iface .eq. 1) then  
                   Trans=0
                   if(ilyr.lt.nlyr-1) then
                      B2_1=B(ilyr+2,ixx,iyy)
                      B1_1=B(ilyr+1,ixx,iyy)
                       Trans = integrate_flux((B1_1+B2_1)/2., L_up_3d(ilyr+2,modulo(ixx,Nx)+1,iyy), &
                           kabs_3d(ilyr+1,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr+1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr+1,ixx,iyy)), dy,  mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_up_3d(ilyr+2,modulo(ixx-2,Nx)+1,iyy), &
                           kabs_3d(ilyr+1,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr+1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr+1,ixx,iyy)), dy, mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_up_3d(ilyr+2,ixx,modulo(iyy,Ny)+1), &
                           kabs_3d(ilyr+1,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr+1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr+1,ixx,iyy)), dx, mu) + &
                           integrate_flux((B1_1+B2_1)/2., L_up_3d(ilyr+2,ixx,modulo(iyy-2,Ny)+1),&
                           kabs_3d(ilyr+1,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr+1,ixx,iyy),&
                           az, bz, cz, abs(dz(ilyr+1,ixx,iyy)), dx, mu)         
                
                      l = (Trans+L_up_3d(ilyr+1,ixx,iyy))/5.
     
                   else
                      l=L_up_3d(nlyr-1,ixx,iyy)
                   end if
                   
                   Absup =l*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, abs(dz(ilyr,ixx,iyy)), mu)
                   Emup = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, abs(dz(ilyr,ixx,iyy)), mu) 
                endif ! iface .eq. 1


                ! #### go to left ####
                if (iface .eq. 2) then

                   areaweight=dz(ilyr,ixx,iyy)/dy

                   Absup1=integrate_abs((B1+B2)/2., L_up_3d(ilyr+1,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr,modulo(ixx-2,Nx)+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dy,  mu)*areaweight
                   Emup1 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight 
                endif !face .eq. 2 


                ! #### go to right #### 
                if (iface .eq. 3) then

                   areaweight=dz(ilyr,ixx,iyy)/dy

                   Absup2=integrate_abs((B1+B2)/2., L_up_3d(ilyr+1,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr,modulo(ixx,Nx)+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight
                   Emup2 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)*areaweight
                endif ! iface .eq. 3


                ! #### go to front #### 
                if (iface .eq. 4) then

                   areaweight=dz(ilyr,ixx,iyy)/dx

                   Absup3=integrate_abs((B1+B2)/2., L_up_3d(ilyr+1,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr,ixx,modulo(iyy-2,Ny)+1), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                   Emup3 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                endif ! iface .eq. 4


                ! ##### go to back #### 
                if (iface .eq. 5) then

                   areaweight=dz(ilyr,ixx,iyy)/dx

                   Absup4=integrate_abs((B1+B2)/2., L_up_3d(ilyr+1,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr,ixx,modulo(iyy,Ny)+1), kabs_3d(ilyr,ixx,iyy),&
                        az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                   Emup4 = -((B1+B2)/2.)*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)*areaweight
                endif ! iface.eq.5

             enddo ! end iface


             if(dx.lt.abs(dz(ilyr,ixx,iyy))) then
                factor=atan(sqrt(kabs_3d(ilyr,ixx,iyy)*dx)-1.4*1.1)*0.105*1.1+(1.-1.57*0.105*1.1)
             end if
             if(dx.gt.abs(dz(ilyr,ixx,iyy))) then
                factor=atan(sqrt(kabs_3d(ilyr,ixx,iyy)*abs(dz(ilyr,ixx,iyy)))-1.4*1.)*0.105*abs(dz(ilyr,ixx,iyy))/dx+(1.-1.57*0.105*abs(dz(ilyr,ixx,iyy))/dx)
             end if
             if(dx.eq.abs(dz(ilyr,ixx,iyy))) then
                factor=atan(sqrt(kabs_3d(ilyr,ixx,iyy)*dx)-1.4)*0.105+(1.-1.57*0.105)
             end if
            
             HR_up = Absup+Emup
             HR_dn = Absdn+Emdn
             HR_up_s = (Absup1+Absup2+Absup3+Absup4+Emup1+Emup2+Emup3+Emup4)/2.*factor
             HR_dn_s = (Absdn1+Absdn2+Absdn3+Absdn4+Emdn1+Emdn2+Emdn3+Emdn4)/2.*factor     
             
             hr_nca_3d_tmp(ilyr,ixx,iyy)  = (HR_up + HR_dn + HR_up_s + HR_dn_s)*pi

#ifndef _XLF
             if(isnan(hr_nca_3d_tmp(ilyr,ixx,iyy))) print *, 'nca shows nan', ilyr, ixx, iyy, hr_nca_3d_tmp(ilyr,ixx,iyy)
#endif
          enddo ! iyy
       enddo ! ixx   
    enddo ! end ilyr
  end subroutine


  ! ################################################################################ 
  ! ####################### Function - integrate_abs ############################### 
  ! # This function integrates the contributiong emission/absorption of a grid box #
  ! ################################################################################ 

  elemental function integrate_abs(B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu)

    real(ireals), intent(in) ::B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu

    real(ireals) :: integrate_abs
    real(ireals) :: factor 
    real(ireals) :: factor1_ab
    real(ireals) :: factor1_bc
    real(ireals) :: factor2_bc
    real(ireals) :: factor3_bc
    real(ireals) :: factor4_bc
    real(ireals) :: sin_mu 
    sin_mu = sqrt(1.-mu*mu)


    if(kabs1*(delta_2/mu).lt.1e-4 .or. kabs1*(delta_1/mu).lt.1e-4) then !Taylor solution
       if (kabs2*(delta_2/mu).lt.1e-4 .or. kabs2*(delta_1/mu).lt.1e-4) then
          integrate_abs =  L*(1. - 1./delta_1 *((b-a) * (1-kabs2*delta_2/mu)&
               + sqrt(1.-mu*mu)/kabs2*((c-b)*kabs2/sqrt(1.-mu*mu))) )
       else
          integrate_abs = L*(1. - 1./delta_1 *((b-a) * exp(-kabs2*delta_2/mu)&
               + sqrt(1.-mu*mu)/kabs2*(1.-exp(-(c-b)*kabs2/sqrt(1.-mu*mu)))) )
       endif
    else
       if (kabs2*(delta_2/mu).lt.1e-4 .or. kabs2*(delta_1/mu).lt.1e-4) then !Taylor solution, 
          integrate_abs = 0.0 
       else 
          factor1_ab = (1-exp(-kabs2*delta_2/mu))*((b-a)*B_planck-( sin_mu/kabs1 * (L-B_planck)&
               * (exp(-kabs1*b/sin_mu)-exp(-kabs1*a/sin_mu))))
          factor1_bc = (c-b)*B_planck
          factor2_bc = -(L-B_planck) * sin_mu/kabs1 * (exp(-kabs1*c/sin_mu)-exp(-kabs1*b/sin_mu))
          if(kabs1.gt.kabs2) then
             factor3_bc = (L-B_planck) * exp(-kabs2*c/sin_mu) * sin_mu/(kabs1-kabs2) *&
                  (exp(-(kabs1-kabs2)*c/sin_mu) - exp(-(kabs1-kabs2)*b/sin_mu))
          else if (kabs1.lt.kabs2) then
             factor3_bc = (L-B_planck) * sin_mu/(kabs1-kabs2) * exp(-kabs1*c/sin_mu)*&
                  (1 -exp((kabs1-kabs2)*(c-b)/sin_mu)) 
          else 
             factor3_bc = - (L-B_planck) * exp(-kabs2*c/sin_mu) * (c-b)
          endif
          factor4_bc = -B_planck* sin_mu/kabs2 *  (1-exp(-(c-b)*kabs2/sin_mu))
          
          integrate_abs = (factor1_ab+factor1_bc+factor2_bc+factor3_bc+factor4_bc)/delta_1
       endif
    endif
  end function integrate_abs



 ! ################################################################################ 
  ! ####################### Function - integrate_flux ############################### 
  ! # This function integrates the contributiong emission/absorption of a grid box #
  ! ################################################################################ 
  elemental function integrate_flux(B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu)

    real(ireals), intent(in) ::B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu

    real(ireals) :: integrate_flux
    real(ireals) :: factor1_ab
    real(ireals) :: factor1_bc
    real(ireals) :: factor3_bc
    real(ireals) :: sin_mu 
    sin_mu = sqrt(1.-mu*mu)

 !   if(kabs1*(delta_2/mu).lt.1e-4 .or. kabs1*(delta_1/mu).lt.1e-4) then !Taylor solution
 !      if (kabs2*(delta_2/mu).lt.1e-4 .or. kabs2*(delta_1/mu).lt.1e-4) then
 !         integrate_abs =  L*(1. - 1./delta_1 *((b-a) * (1-kabs2*delta_2/mu)&
 !              + sqrt(1.-mu*mu)/kabs2*((c-b)*kabs2/sqrt(1.-mu*mu))) )
 !      else
 !         integrate_abs = L*(1. - 1./delta_1 *((b-a) * exp(-kabs2*delta_2/mu)&
 !              + sqrt(1.-mu*mu)/kabs2*(1.-exp(-(c-b)*kabs2/sqrt(1.-mu*mu)))) )
 !      endif
!    else
   if (kabs2*(delta_2/mu).lt.1e-4 .or. kabs2*(delta_1/mu).lt.1e-4) then !Taylor solution, 
      integrate_flux = L 
    else 
       factor1_ab = (b-a)*B_planck-exp(-kabs2*delta_2/mu)*sin_mu/kabs1*(L-B_planck)*&
             (exp(-kabs1*b/sin_mu)-exp(-kabs1*a/sin_mu))
       factor1_bc = (c-b)*B_planck
       if(kabs1.gt.kabs2) then
          factor3_bc = -(L-B_planck) * exp(-kabs2*c/sin_mu) * sin_mu/(kabs1-kabs2) *&
               (exp(-(kabs1-kabs2)*c/sin_mu) - exp(-(kabs1-kabs2)*b/sin_mu))
       else if (kabs1.lt.kabs2) then
          factor3_bc = -(L-B_planck) * sin_mu/(kabs1-kabs2) * exp(-kabs1*c/sin_mu)*&
               (1 -exp((kabs1-kabs2)*(c-b)/sin_mu)) 
       else 
          factor3_bc = + (L-B_planck) * exp(-kabs2*c/sin_mu) * (c-b)
       endif

          
       integrate_flux = (factor1_ab+factor1_bc+factor3_bc)/delta_1
       endif
 !   endif
  end function integrate_flux

  ! ################################################################################ 
  ! ####################### Function - integrate_emis ############################## 
  ! # This function integrates the contributiong emission/absorption of a grid box # 
  ! ################################################################################ 

  elemental function integrate_emis (kabs, a, b, c, delta_1, delta_2, mu) 

    real(ireals), intent(in) :: kabs, a, b, c, delta_1, delta_2, mu
    real(ireals) :: integrate_emis

    if (kabs*(delta_2/mu).lt.1e-4) then  ! Taylor solution
       integrate_emis =  1. - 1./delta_1 *((b-a) * (1-kabs*delta_2/mu)&
            + sqrt(1.-mu*mu)/kabs*((c-b)*kabs/sqrt(1.-mu*mu))) 
    else
       integrate_emis = 1. - 1./delta_1 *((b-a) * exp(-kabs*delta_2/mu)&
            + sqrt(1.-mu*mu)/kabs*(1.-exp(-(c-b)*kabs/sqrt(1.-mu*mu)))) 

    endif
  end function integrate_emis
end module m_ts_nca
