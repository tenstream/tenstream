!----------------------------------------------------------------------------
! This file is part of the tenstream solver.
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
! Copyright (C) 2010-2015  Carolin Klinger, <carolin.klinger@physik.lmu.de>
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2015, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------

module m_ts_nca
  use m_data_parameters, only : ireals,iintegers
  implicit none

  private
  public :: ts_nca

contains
  subroutine ts_nca (dx, dy, dz, B, kabs_3d, Edn, Eup, hr)

    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in), dimension(:,:,:) :: dz, kabs_3d   ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(in), dimension(:,:,:) :: B,Edn,Eup     ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(out),dimension(:,:,:) :: hr            ! dimensions including ghost values (Nlevel,Nx,Ny)

    ! ############## Definition of variables ##################

    integer(iintegers)     :: ilyr          ! counter z levels
    integer(iintegers)     :: ixx           ! counter x grid boxes
    integer(iintegers)     :: iyy           ! counter y grid boxes

    real(ireals)           :: Absup         ! Upwelling, Absorption, bottom
    real(ireals)           :: Absdn         ! Downwelling Absorption, top
    real(ireals)           :: Absup_s       ! Upwelling, Absorption, side
    real(ireals)           :: Absdn_s       ! Downwelling Absorption, side  
    real(ireals)           :: Emup          ! Upwelling, Emission, bottom
    real(ireals)           :: Emdn          ! Downwelling Emission, top
    real(ireals)           :: Emup_s        ! Upwelling, Emission, side
    real(ireals)           :: Emdn_s        ! Downwelling Emission, side  

    real(ireals)           :: HR_T          ! Heating Rate Top Contribution  
    real(ireals)           :: HR_s          ! Heating Rate Side Contribution 

    real(ireals)           :: mu            ! zenith angle

    real(ireals)           :: ax,ay         ! integration boarder for side contributions
    real(ireals)           :: bx,by         ! integration boarder for side contributions 
    real(ireals)           :: cx,cy         ! integration boarder for side contributions
    real(ireals)           :: az1,az        ! integration boarder for side contributions
    real(ireals)           :: bz1,bz        ! integration boarder for side contributions 
    real(ireals)           :: cz1,cz        ! integration boarder for side contributions 


    real(ireals)           :: l      ! Averaged Transmitted Flux, Top Contribution 
    real(ireals)           :: Trans  ! Transmitted Flux, Top Contribution

    real(ireals)           :: L_up_3d(lbound(Eup,1):ubound(Eup,1),lbound(Eup,2):ubound(Eup,2),lbound(Eup,3):ubound(Eup,3))
    real(ireals)           :: L_dn_3d(lbound(Edn,1):ubound(Edn,1),lbound(Edn,2):ubound(Edn,2),lbound(Edn,3):ubound(Edn,3))
    real(ireals)           :: B2   ! Planck General
    real(ireals)           :: B2_1 ! Planck Top Contribution

    integer(iintegers) :: is,ie,js,je,nlyr

    ! ## Fitting Parameters
    real(ireals)           :: afit1 
    real(ireals)           :: bfit1 
    real(ireals)           :: bfit2 
    real(ireals)           :: afit2 
    real(ireals)           :: cfit2 
    real(ireals)           :: factor2 
    real(ireals)           :: factor1 

    real(ireals),parameter :: eps = 0.0001
    real(ireals),parameter :: pi=3.141592653589793

    nlyr = ubound(hr,1)-1
    is = lbound(hr,2) +1
    ie = ubound(hr,2) -1
    js = lbound(hr,3) +1
    je = ubound(hr,3) -1
    

 !   hr(1:nlyr,is:ie,js:je) = Edn(1:nlyr,is:ie,js:je) + Eup(2:nlyr+1,is:ie,js:je) - Eup(1:nlyr,is:ie,js:je) - Edn(2:nlyr+1,is:ie,js:je)

!    do ilyr=1,Nlyr+1
!      print *,'Flux in NCA:',Edn(ilyr,is,js),Eup(ilyr,is,js)
!    enddo

!    return

!    stop 'NCA is not freely available.... please consult Carolin Klinger for an implementation.'


    ! ### get radiance from flux (back converted at the end)
    L_dn_3d=Edn/pi
    L_up_3d=Eup/pi

    ! ### average angle for heating rate calculation
    mu=cos(45._ireals*pi/180._ireals) 
    if(mu.lt.0) then
       mu=mu*(-1._ireals)
    endif
    
    ! ### setting boundaries for radiation integration
    ax=0
    ay=0
    az=0   
    cx=dx
    cy=dy
    az1=0

    ! ################# start 3d calculation ####################
    ! ###########################################################
    do ilyr=1,nlyr    !  loop over all height levels 
       do ixx=is,ie           ! loop over all x gridboxes 
          do iyy=js,je         !  loop over all y gridboxes  

             ! set and reset boundary conditions 
             Emdn=0
             Emdn_s=0
             Absdn=0
             Absdn_s=0
             Emup=0
             Emup_s=0
             Absup=0
             Absup_s=0

             ! ### setting boundaries for radiation integration
             cz=dz(ilyr,ixx,iyy)

             bx=cx-(dz(ilyr,ixx,iyy)/mu)*sqrt(1._ireals-mu*mu) 
             by=cy-(dz(ilyr,ixx,iyy)/mu)*sqrt(1._ireals-mu*mu) 
             bz=cz-(dx/sqrt(1.-mu*mu))*mu 

             if (abs(bx-cx).le.eps .or. abs(by-cy).le.eps) then
                bx=ax
                by=ay
             endif
             if (abs(bz-cz).le.eps) then
                bz=az
             endif
             if(bx.lt.0._ireals .or. by.lt.0._ireals)then
                bx=ax
                by=ay
             endif
             if(bz.lt.0)then
                bz=az
             endif
             if(bx.lt.eps .or. by.lt.eps)then
                bx=0._ireals
                by=0._ireals
             endif
             if(bz.lt.eps)then
                bz=0._ireals
             endif


             ! ## Average Planck of layer
             B2=(B(ilyr+1,ixx,iyy)+B(ilyr,ixx,iyy))/2._ireals

             ! if lower face of gridbox 
             if(ilyr.gt.1) then

               ! ## Special treatment for boundaries and Planck, one layer up/down
                B2_1=(B(ilyr,ixx,iyy)+B(ilyr-1,ixx,iyy))/2._ireals
      
                az1=0._ireals
                cz1=dz(ilyr-1,ixx,iyy)
                bz1=cz1-(dx/sqrt(1.-mu*mu))*mu 

                if (abs(bz1-cz1).le.eps) then
                   bz1=az1
                endif
                if(bz1.lt.0._ireals)then
                   bz1=az1
                endif
                if(bz1.lt.eps)then
                   bz1=0._ireals
                endif
                
                !### averaged transmitted flux
              Trans = integrate_flux((B(ilyr,ixx+1,iyy)+B(ilyr-1,ixx+1,iyy))/2._ireals, B2_1, L_dn_3d(ilyr-1,ixx+1,iyy), &
                     kabs_3d(ilyr-1,ixx+1,iyy), kabs_3d(ilyr-1,ixx,iyy),&
                     bz1, cz1, dz(ilyr-1,ixx,iyy), dy,  mu) + &
                     integrate_flux((B(ilyr,ixx-1,iyy)+B(ilyr-1,ixx-1,iyy))/2._ireals, B2_1, L_dn_3d(ilyr-1,ixx-1,iyy),&
                     kabs_3d(ilyr-1,ixx-1,iyy), kabs_3d(ilyr-1,ixx,iyy),&
                     bz1, cz1, dz(ilyr-1,ixx,iyy), dy, mu) + &                
                     integrate_flux((B(ilyr,ixx,iyy-1)+B(ilyr-1,ixx,iyy-1))/2._ireals, B2_1, L_dn_3d(ilyr-1,ixx,iyy-1),&
                     kabs_3d(ilyr-1,ixx,iyy-1), kabs_3d(ilyr-1,ixx,iyy),&
                     bz1, cz1, dz(ilyr-1,ixx,iyy), dx, mu) + &
                     integrate_flux((B(ilyr,ixx,iyy+1)+B(ilyr-1,ixx,iyy+1))/2._ireals, B2_1, L_dn_3d(ilyr-1,ixx,iyy+1),&
                     kabs_3d(ilyr-1,ixx,iyy+1), kabs_3d(ilyr-1,ixx,iyy),&
                     bz1, cz1, dz(ilyr-1,ixx,iyy), dx, mu)


                l = (Trans+L_dn_3d(ilyr,ixx,iyy))/5._ireals

             else
                l = L_dn_3d(ilyr,ixx,iyy)     
             end if

             !### top/bottom contribution
             Absdn =  l*integrate_emis(kabs_3d(ilyr,ixx,iyy),&
                  ax, bx, cx, dx, dz(ilyr,ixx,iyy), mu) 
             Emdn = -B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, dz(ilyr,ixx,iyy), mu)
             Absdn =  Absdn + l*integrate_emis(kabs_3d(ilyr,ixx,iyy),&
                  ay, by, cy, dy, dz(ilyr,ixx,iyy), mu) 
             Emdn = Emdn - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), ay, by, cy, dy, dz(ilyr,ixx,iyy), mu)

             ! #### side contribution - go to left ####     
             Absdn_s = Absdn_s+integrate_abs ((B(ilyr+1,ixx+1,iyy)+B(ilyr,ixx+1,iyy))/2, L_dn_3d(ilyr,ixx+1,iyy), &
                  kabs_3d(ilyr,ixx+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dy, mu) 
             Emdn_s = -B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)


             ! #### side contribution - go to right #### 
             Absdn_s = Absdn_s + integrate_abs((B(ilyr+1,ixx-1,iyy)+B(ilyr,ixx-1,iyy))/2,  L_dn_3d(ilyr,ixx-1,iyy), &
                  kabs_3d(ilyr,ixx-1,iyy), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dy,  mu) 
             Emdn_s = Emdn_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu) 


             ! #### side contribution - go to front ####    
             Absdn_s = Absdn_s + integrate_abs((B(ilyr+1,ixx,iyy-1)+B(ilyr,ixx,iyy-1))/2,  L_dn_3d(ilyr,ixx,iyy-1), &
                  kabs_3d(ilyr,ixx,iyy-1), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)
             Emdn_s = Emdn_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)    


             ! ##### side contribution - go to back ####           
             Absdn_s = Absdn_s + integrate_abs((B(ilyr+1,ixx,iyy+1)+B(ilyr,ixx,iyy+1))/2,  L_dn_3d(ilyr,ixx,iyy+1), &
                  kabs_3d(ilyr,ixx,iyy+1), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)
             Emdn_s = Emdn_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)



             ! ## upwelling           
             ! if lower face of gridbox 
             if(ilyr.lt.nlyr) then

                ! ## Special treatment for boundaries and Planck, one layer up/down
                B2_1=(B(ilyr+1,ixx,iyy)+B(ilyr+2,ixx,iyy))/2._ireals

                az1=0._ireals
                cz1=dz(ilyr+1,ixx,iyy)

                ! ## set parameters for integration 
                bz1=cz1-(dx/sqrt(1.-mu*mu))*mu 

                if (abs(bz1-cz1).le.eps) then
                   bz1=az1
                endif
                if(bz1.lt.0._ireals)then
                   bz1=az1
                endif
                if(bz1.lt.eps)then
                   bz1=0._ireals
                endif

                !### averaged transmitted flux
                Trans = integrate_flux((B(ilyr+1,ixx+1,iyy)+B(ilyr+2,ixx+1,iyy))/2._ireals, B2_1, L_up_3d(ilyr+2,ixx+1,iyy), &
                     kabs_3d(ilyr+1,ixx+1,iyy), kabs_3d(ilyr+1,ixx,iyy),&
                     bz1, cz1, dz(ilyr+1,ixx,iyy), dy,  mu) + &
                     integrate_flux((B(ilyr+1,ixx-1,iyy)+B(ilyr+2,ixx-1,iyy))/2._ireals, B2_1, L_up_3d(ilyr+2,ixx-1,iyy), &
                     kabs_3d(ilyr+1,ixx-1,iyy), kabs_3d(ilyr+1,ixx,iyy),&
                     bz1, cz1, dz(ilyr+1,ixx,iyy), dy, mu) + &
                     integrate_flux((B(ilyr+1,ixx,iyy+1)+B(ilyr+2,ixx,iyy+1))/2._ireals, B2_1, L_up_3d(ilyr+2,ixx,iyy+1), &
                     kabs_3d(ilyr+1,ixx,iyy+1), kabs_3d(ilyr+1,ixx,iyy),&
                     bz1, cz1, dz(ilyr+1,ixx,iyy), dx, mu) + &
                     integrate_flux((B(ilyr+1,ixx,iyy-1)+B(ilyr+2,ixx,iyy-1))/2._ireals, B2_1, L_up_3d(ilyr+2,ixx,iyy-1),&
                     kabs_3d(ilyr+1,ixx,iyy-1), kabs_3d(ilyr+1,ixx,iyy),&
                     bz1, cz1, dz(ilyr+1,ixx,iyy), dx, mu)
                l = (Trans+L_up_3d(ilyr+1,ixx,iyy))/5._ireals
             else
                l=L_up_3d(nlyr+1,ixx,iyy)
             end if

             !### top/bottom contribution
             Absup = l*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, dz(ilyr,ixx,iyy), mu)
             Emup  = -B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), ax, bx, cx, dx, dz(ilyr,ixx,iyy), mu) 
             Absup = Absup + l*integrate_emis(kabs_3d(ilyr,ixx,iyy), ay, by, cy, dy, dz(ilyr,ixx,iyy), mu)
             Emup  = Emup - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), ay, by, cy, dy, dz(ilyr,ixx,iyy), mu) 


             ! #### side contribution - go to left ####
            Absup_s = integrate_abs((B(ilyr+1,ixx+1,iyy)+B(ilyr,ixx+1,iyy))/2, L_up_3d(ilyr+1,ixx+1,iyy), &
                  kabs_3d(ilyr,ixx+1,iyy), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)
             Emup_s = -B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz,  dz(ilyr,ixx,iyy), dy, mu)
             
             ! #### side contribution - go to right #### 
             Absup_s = Absup_s + integrate_abs((B(ilyr+1,ixx-1,iyy)+B(ilyr,ixx-1,iyy))/2, L_up_3d(ilyr+1,ixx-1,iyy), &
                  kabs_3d(ilyr,ixx-1,iyy), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)
             Emup_s = Emup_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dy, mu)
            
             ! #### side contribution - go to front #### 
             Absup_s = Absup_s + integrate_abs((B(ilyr+1,ixx,iyy+1)+B(ilyr,ixx,iyy+1))/2, L_up_3d(ilyr+1,ixx,iyy+1), &
                  kabs_3d(ilyr,ixx,iyy+1), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz,dz(ilyr,ixx,iyy), dx, mu)
             Emup_s = Emup_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)
            
             ! ##### side contribution - go to back #### 
             Absup_s = Absup_s + integrate_abs((B(ilyr+1,ixx,iyy-1)+B(ilyr,ixx,iyy-1))/2, L_up_3d(ilyr+1,ixx,iyy-1),&
                  kabs_3d(ilyr,ixx,iyy-1), kabs_3d(ilyr,ixx,iyy),&
                  az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)
             Emup_s = Emup_s - B2*integrate_emis(kabs_3d(ilyr,ixx,iyy), az, bz, cz, dz(ilyr,ixx,iyy), dx, mu)

             ! ### Correction fits
             afit1 = atan(log(dx/dz(ilyr,ixx,iyy))-0.24_ireals)*(0.39_ireals) 
             bfit1 = atan(log(dx/dz(ilyr,ixx,iyy))-0.3_ireals)*(-0.67_ireals)+1.28_ireals 
             afit2 = atan(dx/dz(ilyr,ixx,iyy)*7.35_ireals)*0.6_ireals-0.35_ireals  
             cfit2 = atan(dx/dz(ilyr,ixx,iyy)-0.36_ireals)*1.14_ireals+(2.02_ireals)/exp(dx/dz(ilyr,ixx,iyy)*(0.48_ireals))-0.84_ireals  
             bfit2 = (1.1-cfit2)/(pi/2._ireals)

             factor2=atan(kabs_3d(ilyr,ixx,iyy)*dx*afit2)*bfit2+cfit2

             if(dx.le.dz(ilyr,ixx,iyy)) then
                factor1=atan(kabs_3d(ilyr,ixx,iyy)*dx)*afit1+bfit1
             else 
                factor1=atan(kabs_3d(ilyr,ixx,iyy)*dz(ilyr,ixx,iyy))*afit1+bfit1
             end if

             !### scale and summarize heating
             HR_T = (Absup+Emup+Absdn+Emdn)*pi/2._ireals*dy/dy/dx/dz(ilyr,ixx,iyy)/factor2
             HR_s = (Absup_s+Emup_s+Absdn_s+Emdn_s)*pi/2._ireals*dy/dy/dx/dz(ilyr,ixx,iyy)/factor1 

             hr(ilyr,ixx,iyy)  = (HR_T + HR_s)*dz(ilyr,ixx,iyy)
 
#ifndef _XLF
             if(isnan(hr(ilyr,ixx,iyy))) print *, 'nca shows nan', ilyr, ixx, iyy, hr(ilyr,ixx,iyy)
#endif

          enddo ! iyy
       enddo ! ixx   
    enddo ! end ilyr


    !   hr_nca_3d_tmp = -9999999
    !   stop 'NCA not freely available.... please consult Carolin Klinger for an implementation.'
  end subroutine ts_nca



  ! ################################################################################ 
  ! ####################### Function - integrate_abs ############################### 
  ! # This function integrates the contributiong emission/absorption of a grid box #
  ! ################################################################################ 

  elemental function integrate_abs(B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu)

    real(ireals), intent(in) ::B_planck, L, kabs1, kabs2, a, b, c, delta_1, delta_2, mu

    real(ireals) :: integrate_abs
    real(ireals) :: factor 
    real(ireals) :: factor_ab
    real(ireals) :: factor1_ab
    real(ireals) :: factor1_bc
    real(ireals) :: factor2_bc
    real(ireals) :: factor3_bc
    real(ireals) :: factor4_bc
    real(ireals) :: sin_mu 
    real(ireals) :: exp1, exp2, exp3, exp4, exp5
    real(ireals) :: LB

    sin_mu = sqrt(1._ireals-mu*mu)
    factor_ab = 0
    factor1_bc = 0
    factor2_bc = 0
    factor3_bc = 0
    factor4_bc = 0
    integrate_abs = 0 
    LB = L - B_planck

    exp1 = exp(-kabs2*delta_2/sin_mu)
    exp2 = exp(-kabs2*delta_1/mu)
    exp3 = exp(-kabs1*b/mu)
    exp4 = exp(-kabs1*a/mu)
    exp5 = exp(-kabs1*c/mu)

    if (kabs1*c.lt.1e-4_ireals .or. kabs1*c.lt.1e-4_ireals)   then !Taylor solution, 
       if (kabs2*c.lt.1e-4_ireals .or. kabs2*c.lt.1e-4_ireals)   then !Taylor solution, 
          integrate_abs = 0.0 
       else
          integrate_abs = L*((1._ireals-exp1)*(b-a) -(c-b)-exp2*(c-b))
       endif
    else
       if(kabs1*c.lt.1e-4_ireals)then
          factor_ab = (1._ireals-exp1)*L*(b-a)
          factor1_bc = -(LB)*(c-b)
       else
          factor_ab = (1._ireals-exp1)*(B_planck*(b-a) - (LB)*mu/kabs1* (exp3-exp4))
          factor1_bc = -(LB)*mu/kabs1*(exp5-exp3)
       endif

       factor2_bc = B_planck*(c-b)

       if (kabs2*c.lt.1e-4_ireals) then
          factor3_bc = -B_planck*exp2*(c-b)
       else
          factor3_bc = -B_planck*mu/kabs2*(exp(kabs2*(c-delta_1)/mu)-exp(kabs2*(b-delta_1)/mu))  
       endif

       if ((kabs1*c.lt.1e-4_ireals).or.(abs(kabs2-kabs1).lt.1e-4_ireals)) then
          factor4_bc = -(LB)*exp2*(c-b)
       else
          factor4_bc = -(LB)*mu/(kabs2-kabs1)*&
               (exp(((kabs2-kabs1)*c-kabs2*delta_1)/mu)-exp(((kabs2-kabs1)*b-kabs2*delta_1)/mu))  
       endif

    endif

    integrate_abs = (factor_ab+factor1_bc+factor2_bc+factor3_bc+factor4_bc)

end function integrate_abs


! ################################################################################ 
! ####################### Function - integrate_flux ############################### 
! # This function integrates the contributiong emission/absorption of a grid box #
! ################################################################################ 
elemental function integrate_flux(B_planck1, B_planck2, L, kabs1, kabs2, b, c, delta_1, delta_2, mu)
  
  real(ireals), intent(in) ::B_planck1,B_planck2, L, kabs1, kabs2, b, c, delta_1, delta_2, mu
  
  real(ireals) :: integrate_flux
  real(ireals) :: factor1_bc
  real(ireals) :: factor2_bc
  real(ireals) :: factor3_bc
  real(ireals) :: sin_mu
  real(ireals) :: exp1
  real(ireals) :: LB

  factor1_bc = 0
  factor2_bc = 0
  factor3_bc = 0
  integrate_flux = 0 
  sin_mu = sqrt(1._ireals-mu*mu)
  
  LB = L - B_planck1
  
  exp1 = exp(-kabs2*delta_2/sin_mu)
  
  factor1_bc = B_planck2 * (c-b)

  if(kabs1.lt.1e-4_ireals)then
     factor2_bc =  (B_planck1 - B_planck2) * exp1 * (c-b)
  else
     factor2_bc =  (B_planck1 - B_planck2) * exp1 * mu/kabs2 * &
          (exp(-kabs2*c/mu) - exp(-kabs2*b/mu))
  endif
  if ((kabs1.lt.1e-4_ireals).or.(abs(kabs2-kabs1).lt.1e-4_ireals)) then
     factor3_bc = (LB) * (c-b) * exp(-kabs2*delta_1/mu)
  else
     factor3_bc = (LB) * mu/(kabs2-kabs1) *&
          (exp(((kabs2-kabs1)*c-(kabs2*delta_1))/mu) - exp(((kabs2-kabs1)*b-(kabs2*delta_1))/mu))
  endif
  
  integrate_flux = (factor1_bc+factor2_bc+factor3_bc)/(c-b)

end function integrate_flux

! ################################################################################ 
! ####################### Function - integrate_emis ############################## 
! # This function integrates the contributiong emission/absorption of a grid box # 
! ################################################################################ 

elemental function integrate_emis (kabs, a, b, c, delta_1, delta_2, mu) 
  
  real(ireals), intent(in) :: kabs, a, b, c, delta_1, delta_2, mu
  real(ireals) :: integrate_emis
  real(ireals) :: exp1,sin_mu
  
  integrate_emis = 0
  sin_mu=sqrt(1._ireals-mu*mu)
  exp1 = exp(-kabs*delta_2/sin_mu)
  
  if(kabs*c.lt.1e-4_ireals)then
     integrate_emis = (b-a)*(1._ireals-exp1) +(c-b)*(1.-exp(-kabs*delta_1/mu))
  else
     integrate_emis = ((b-a)*(1._ireals-exp1) + (c-b) - mu/kabs*&
          (exp(kabs*(c-delta_1)/mu)-exp(kabs*(b-delta_1)/mu)))
  endif
  
end function integrate_emis

end module m_ts_nca
