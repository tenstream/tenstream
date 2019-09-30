!----------------------------------------------------------------------------
! Copyright (C) 2010-2019  Carolin Klinger, <carolin.klinger@physik.lmu.de>
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2019, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------

module m_plexrt_nca
#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR
  use m_data_parameters, only : ireals, iintegers, mpiint, default_str_len
  implicit none

  private
  public :: plexrt_nca_init, plexrt_nca

  real(ireals), dimension(:,:), allocatable :: eps_tab_side, eps_tab_top
  real(ireals), dimension(:,:), allocatable :: corr_tab_side, corr_tab_top
  real(ireals), dimension(:)  , allocatable :: tau_hx, tau_z
  real(ireals), dimension(:)  , allocatable :: var_1, var_2

contains
  subroutine plexrt_nca_init()
    integer(iintegers) :: i,j
    integer :: funit, ierr
    logical :: lflg

    integer(iintegers), parameter :: ntau=16, n1=9, n2=36
    character(len=*), parameter :: eps_tab_side_fname = "lookup_side_triangle.dat"
    character(len=*), parameter :: eps_tab_top_fname  = "lookup_top_triangle.dat"
    character(len=*), parameter :: cor_tab_side_fname = "lookup_correct_triangle_side.dat"
    character(len=*), parameter :: cor_tab_top_fname  = "lookup_correct_triangle_top.dat"
    character(len=default_str_len) :: lut_dir

    lut_dir = '.'
    call PetscInitialized(lflg, ierr)
    if(lflg) then
      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-nca_lut_dir', lut_dir, lflg, ierr); call CHKERR(ierr)
    endif

    ! lookup table for emissivity
    call load_table(trim(lut_dir)//'/'//eps_tab_side_fname, ntau, ntau, eps_tab_side)
    call load_table(trim(lut_dir)//'/'//eps_tab_top_fname, ntau, ntau, eps_tab_top)

    ! lookup table for correction
    call load_table(trim(lut_dir)//'/'//cor_tab_side_fname, n1, n2, corr_tab_side)
    call load_table(trim(lut_dir)//'/'//cor_tab_top_fname, n1, n2, corr_tab_top)

    allocate(tau_hx(ntau))
    allocate(tau_z (ntau))

    tau_hx(:) = [0.0001, 0.00025, 0.000625, 0.0015625, 0.00390625, &
         0.00976562, 0.0244141, 0.0610352, 0.152588, 0.38147, 0.953674, &
         2.38419, 5.96046, 14.9012, 37.2529, 93.1323]

    tau_z(:) = [0.0001, 0.00025, 0.000625, 0.0015625, 0.00390625, &
         0.00976562, 0.0244141, 0.0610352, 0.152588, 0.38147, 0.953674, &
         2.38419, 5.96046, 14.9012, 37.2529, 93.1323]

    allocate(var_1(n1+2))
    allocate(var_2(n2+2))
    !var_1(:) = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
    var_1(:) = [0.11547, 0.23094, 0.34641, 0.57735, 1.1547, 2.3094, &
         3.4641, 5.7735, 11.547]

    var_2(:) = [0.00015, 0.000225, 0.000338, 0.000506, 0.000759, &
         0.001139, 0.001709, 0.002563, 0.003844, 0.005766, 0.00865, &
         0.012975, 0.019462, 0.029193, 0.043789, 0.065684, 0.098526, &
         0.147789, 0.221684, 0.332526, 0.498789, 0.748183, 1.12227, &
         1.68341, 2.52512, 3.78768, 5.68151, 8.52227, 12.7834, 19.1751, &
         28.7627, 43.144, 64.716, 97.074, 145.611, 218.416]
  contains
    subroutine load_table(fname, n1, n2, arr)
      character(len=*), intent(in) :: fname
      integer(iintegers), intent(in) :: n1, n2
      real(ireals), allocatable, intent(inout) :: arr(:,:)
      logical :: lexists
      real(ireals) :: tmp

      if(.not.allocated(arr)) then
        allocate(arr(n1,n2))

        inquire(file=trim(fname), exist=lexists)
        if(.not.lexists) call CHKERR(1_mpiint, 'Could not find NCA LUT : '//trim(fname))

        open (newunit=funit, file=trim(fname), status="old", action="read")
        do i = 1, n1
          do j = 1, n2
            read(funit,*) tmp, tmp, arr(i,j)
          end do
        end do
        close(funit)
      endif
    end subroutine
  end subroutine


    ! NCA geometry for Wedges with vertices (ABC, DEF):
    !
    !
    !         F
    !        / \
    !    dx2/   \dx3
    !      /     \
    !     /  atop \
    !    /   dx1   \
    !   D __________E
    !   |           |
    !   |     |     |
    !   | a2  |  a3 |
    !   |     |     |
    !   |           |
    !   |     C     |
    !   |    / \    |
    !   |   /   \   |
    !   |  /     \  |
    !   | /       \ |
    !   |/  abot   \|
    !   A _________ B

    ! NCA flux/optprop info for Wedges with vertices (ABC, DEF), box on top has vertices (DEF, GHI):
    !      Base Wedge(1)      Neighbor Wedge
    !
    !         I                         F
    !        / \                       / \
    !       /   \                     /   \
    !      /     \                   /     \
    !     /       \                 /       \
    !    /         \               /         \
    !   G __________H             D __________E
    !
    !   |  kabs_top |
    !   |           |
    !   |     F     |                   F
    !   |    / \    |                  / \
    !   |   /   \   |                 /   \
    !   |  /     \  |                /     \
    !     /       \                 /       \
    !    /   Etop  \               /  Etop_s \
    !   D __________E             D __________E
    !   |           |             |           |
    !   |   kabs    |             |   kabs_s  |
    !   |     |     |             |     |     |
    !   |     |     |             |     |     |
    !   |           |             |           |
    !   |     C     |             |     C     |
    !   |    / \    |             |    / \    |
    !   |   /   \   |             |   /   \   |
    !   |  /     \  |             |  /     \  |
    !   | /       \ |             | /       \ |
    !   |/   Ebot  \|             |/  Ebot_s \|
    !   A _________ B             A _________ B


    subroutine plexrt_nca (dx1, dx2, dx3, dz, atop, abot, a1, a2, a3, v, &
        base_info, side_info, hr)
      real(ireals), intent(in) :: dx1, dx2, dx3          ! edge lengths of triangle: dx1, dx2, dx3
      real(ireals), intent(in) :: a1, a2, a3, atop, abot ! area of side faces and top and bot faces
      real(ireals), intent(in) :: dz, v                  ! height of grid box and volume

      ! info on this voxel dim(7) ( Edn_top, Btop, Eup_bot, Bbot, kabs, kabs_top, kabs_bot)
      real(ireals), intent(in), dimension(:) :: base_info

      ! info for each of the side voxels dim(3 * 5) ( Edn_top, Eup_top, Edn_bot, Eup_bot, kabs)
      real(ireals), intent(in), dimension(:) :: side_info
      real(ireals), intent(out) :: hr        ! new 3D heating rate in the voxel

      ! ############## Definition of variables ##################

      real(ireals)           :: Absup     ! Upwelling, Absorption, bottom
      real(ireals)           :: Absdn     ! Downwelling Absorption, top
      real(ireals)           :: Absups   ! Upwelling, Absorption, side
      real(ireals)           :: Absdns   ! Downwelling Absorption, side
      real(ireals)           :: Emup      ! Upwelling, Emission, bottom
      real(ireals)           :: Emdn      ! Downwelling Emission, top
      real(ireals)           :: Emups    ! Upwelling, Emission, side
      real(ireals)           :: Emdns    ! Downwelling Emission, side
      real(ireals)           :: B         ! Average Planck of Layer

      real(ireals)           :: asp       ! Aspect Ratio

      real(ireals)           :: tauz      ! Vertical optical depth
      real(ireals)           :: tauhx     ! horizontal optical depth (along height of triangle)

      real(ireals)           :: eps_top   ! Emissivity from Top Direction
      real(ireals)           :: eps_bot   ! Emissivity from Base Direction
      real(ireals)           :: eps_s1    ! Emissivity from Side Face 1
      real(ireals)           :: eps_s2    ! Emissivity from Side Face 2
      real(ireals)           :: eps_s3    ! Emissivity from Side Face 3

      real(ireals)           :: f_final_t, f_final_b, f_final_s1, f_final_s2, f_final_s3  ! Correction factors
      real(ireals)           :: l, Trans

      ! factors and weights
      real(ireals)           :: f1
      real(ireals)           :: f2
      real(ireals)           :: w1,w2
      real(ireals)           :: wa,wb,wc

      ! Variables from main program
      real(ireals),parameter :: pi=3.141592653589793

      if(.not. all( [ &
        allocated(eps_tab_side), &
        allocated(eps_tab_top), &
        allocated(corr_tab_side), &
        allocated(corr_tab_top) ] ) ) &
          call CHKERR(1_mpiint, 'called NCA but seems the NCA LUT`s havent been loaded ...'// &
          ' Try to call plexrt_nca_init() first')

      ! Dummy computations to get rid of unused compiler warnings
      !hr = dx1 + dx2 + dx3 + a1 + a2 + a3 + atop + abot + dz + v + base_info(1) + side_info(1)


      ! Get variables form input array and convert flux to radiance
      ! (necessary because of planck - will be converted back at the end)

      associate( kabs        => base_info(1),     &
                 kabs_top    => base_info(2),     &
                 Ldn_top     => base_info(3)/pi,  &
                 Btop        => base_info(4),     &
                 kabs_bot    => base_info(5),     &
                 Lup_bot     => base_info(6)/pi,  &
                 Bbot        => base_info(7),     &
                 kabs_s1     => side_info( 1),    &
                 Ldn_top_s1  => side_info( 2)/pi, &
                 Lup_top_s1  => side_info( 3)/pi, &
                 Ldn_bot_s1  => side_info( 4)/pi, &
                 Lup_bot_s1  => side_info( 5)/pi, &
                 kabs_s2     => side_info( 6),    &
                 Ldn_top_s2  => side_info( 7)/pi, &
                 Lup_top_s2  => side_info( 8)/pi, &
                 Ldn_bot_s2  => side_info( 9)/pi, &
                 Lup_bot_s2  => side_info(10)/pi, &
                 kabs_s3     => side_info(11)   , &
                 Ldn_top_s3  => side_info(12)/pi, &
                 Lup_top_s3  => side_info(13)/pi, &
                 Ldn_bot_s3  => side_info(14)/pi, &
                 Lup_bot_s3  => side_info(15)/pi )



      ! ###########################################################
      ! ################# start 3d calculation ####################
      ! ###########################################################

      ! set and reset boundary conditions
      Emdn   = 0._ireals
      Emdns  = 0._ireals
      Absdn  = 0._ireals
      Absdns = 0._ireals
      Emup   = 0._ireals
      Emups  = 0._ireals
      Absup  = 0._ireals
      Absups = 0._ireals


      ! ## Average Planck of layer
      B = ( Btop + Bbot ) / 2._ireals

      asp = dz / (( 2._ireals / (( dx1 + dx2 + dx3 ) / 3._ireals )) * atop )
      ! Lookup table from  0.11 to 11, for asp=dz/hc
      if( asp > 11._ireals )then
         asp = 11._ireals
      else if( asp < 0.11 )then
         asp = 0.11
      end if

      !! emissivity from lookup table - linear interpolation
      ! get tau in vertical dimension (same for all considerations)
      tauz = kabs * dz

      
      ! find emissivity for top face
      ! get tau of average side edge length
      tauhx     = kabs * ( dx1 + dx2 + dx3 ) / 3 * 0.86603
      eps_top   = interpol_emis(tauhx, tauz,  eps_tab_top)
      f_final_t = interpol_2d(asp, tauz, corr_tab_top)
 
      
      asp = dz / (( 2._ireals / (( dx1 + dx2 + dx3 ) / 3._ireals )) * abot )
      ! Lookup table from  0.11 to 11, for asp=dz/hc
      if( asp > 11._ireals )then
         asp = 11._ireals
      else if( asp < 0.11 )then
         asp = 0.11
      end if

      ! find emissivity for base face
      ! get tau of average side edge length
      tauhx     = kabs * ( dx1 + dx2 + dx3 ) / 3 * 0.86603
      eps_bot   = interpol_emis(tauhx, tauz,  eps_tab_top)
      f_final_b = interpol_2d(asp, tauz, corr_tab_top)

      ! find emissivity for side face 1 (index 3)
      ! get tau for specific side edge length
      tauhx      = kabs * dx1 * 0.86603
      eps_s1     = interpol_emis(tauhx, tauz, eps_tab_side)
      f_final_s1 = interpol_2d(asp, tauhx, corr_tab_side)

      ! find emissivity for side face 2 (index 4)
      ! get tau for specific side edge length
      tauhx      = kabs * dx2 * 0.86603
      eps_s2     = interpol_emis(tauhx, tauz, eps_tab_side)
      f_final_s2 = interpol_2d(asp, tauhx, corr_tab_side)

      ! find emissivity for side face 3 (index 5)
      ! get tau for specific side edge length
      tauhx      = kabs * dx3 * 0.86603
      eps_s3     = interpol_emis(tauhx, tauz, eps_tab_side)
      f_final_s3 = interpol_2d(asp, tauhx, corr_tab_side)


      !#############################################################################################################

      ! #### downwelling ####

      ! #### top ####
      ! Fit was made form 0.1 to 10, for asp = dz / hc
      ! get aspect ratio of

      asp = dz / (( 2._ireals / (( dx1 + dx2 + dx3 ) / 3._ireals )) * atop )
      if( asp >  10._ireals )then
         asp = 10._ireals
      else if( asp < 0.1 )then
         asp = 0.1
      end if

      ! get weights for incoming top flux
      wa = atan( asp * 1.29 ) * ( -0.75 ) + 1.21
      wb = ( asp**0.027 ) * ( -7.98 ) + asp * ( -0.01 ) + atan( asp * 0.11 ) + 7.36
      wc = ( asp**0.49 ) * ( 1.46 )   + asp * ( -0.25 ) + atan( asp * (-0.29)) -0.12
      w1 = atan( kabs_top * dz * wa ) * wb + wc
      w2 = 1._ireals - w1

      Trans = Ldn_top_s1 + Ldn_top_s2 + Ldn_top_s3  ! side face
      l = w1 * Trans / 3._ireals + Ldn_top * w2     ! weight average of side fluxes and top flux


      ! Get Emission and Absorption from each face
      Absdn = l  * atop * eps_top * f_final_t
      Emdn  = -B * atop * eps_top * f_final_t

      ! #### First side face ####
      f1 = atan( kabs_s1 * dz *( -2.08 / ( dz / dx1 ) )) * 0.31192 + 0.49
      if( f1 < 0._ireals )then
         f1 = 0._ireals
      end if
      f2 = 1._ireals-f1

      Absdns = ( f1 * Ldn_top_s1 + f2 * Ldn_bot_s1 ) * a1 * eps_s1 * f_final_s1
      Emdns  = -B * a1 * eps_s1 * f_final_s1

      
      ! #### Second side face ####
      f1 = atan( kabs_s2 * dz * ( -2.08 / ( dz / dx2 ) )) * 0.31192 + 0.49
      if(f1 < 0._ireals)then
         f1 = 0._ireals
      end if
      f2 = 1._ireals-f1

      Absdns = Absdns + ( f1 * Ldn_top_s2 + f2 * Ldn_bot_s2 ) * a2 * eps_s2 * f_final_s2
      Emdns  = Emdns - B * a2 * eps_s2 * f_final_s2


      ! #### Third side face ####
      f1 = atan( kabs_s3 * dz * ( -2.08  / ( dz / dx3 ) )) * 0.31192 + 0.49
      if( f1 < 0._ireals )then
         f1 = 0._ireals
      end if
      f2 = 1._ireals-f1

      Absdns = Absdns + ( f1 * Ldn_top_s3 + f2 * Ldn_bot_s3 ) * a3 * eps_s3 * f_final_s3
      Emdns  = Emdns  - B * a3 * eps_s3 * f_final_s3
 

      !#############################################################################################################


      ! #### upwelling ####

      ! #### bot ####
      ! Fit was made form 0.1 to 10, for asp=dz/hc

      asp = dz / (( 2._ireals / (( dx1 + dx2 + dx3 ) / 3._ireals )) *abot )
      if( asp > 10._ireals )then
         asp = 10._ireals
      else if( asp < 0.1 )then
         asp = 0.1
      end if

      ! get weights for incoming bot flux
      wa = atan( asp * 1.29 ) * ( -0.75 ) + 1.21
      wb = ( asp**0.027 ) * ( -7.98 ) + asp * ( -0.01 ) + atan( asp * 0.11 ) + 7.36
      wc = ( asp**0.49 ) * ( 1.46 ) + asp * ( -0.25 ) + atan( asp * (-0.29) ) - 0.12
      w1 = atan( kabs_bot * dz * wa) * wb + wc     ! kabs of grid box (index 0) is needed
      w2 = 1._ireals - w1

      Trans = Lup_bot_s1 + Lup_bot_s2 + Lup_bot_s3 ! side face
      l = w1 * Trans / 3._ireals + Lup_bot * w2


      ! Get Emission and Absorption from each face
      ! ## Face 1
      Absup = l  * abot * eps_bot * f_final_b
      Emup  = -B * abot * eps_bot * f_final_b

      ! ####  First side face ####
      f1 = atan( kabs_s1 * dz * ( -2.08 / ( dz / dx1 ) )) * 0.31192 + 0.49
      if( f1 < 0._ireals )then
         f1 = 0._ireals
      end if
      f2 = 1._ireals - f1

      Absups = ( f1 *Lup_bot_s1 + f2 * Lup_top_s1 ) * a1 * eps_s1 * f_final_s1
      Emups  = -B * a1 * eps_s1 * f_final_s1


      ! #### Second side face ####
      f1 = atan( kabs_s2 * dz *( -2.08 / ( dz / dx2 ) )) * 0.31192 + 0.49
      if( f1 < 0._ireals )then
         f1 = 0._ireals
      end if
      f2 = 1._ireals - f1

      Absups = Absups + ( f1 * Lup_bot_s2 + f2 * Lup_top_s2 ) * a2 * eps_s2 * f_final_s2
      Emups  = Emups - B * a2 * eps_s2 * f_final_s2


      ! #### Third side face ####
      f1 = atan( kabs_s3 * dz * ( -2.08 / ( dz / dx3 ) )) * 0.31192 + 0.49
      if( f1 < 0._ireals )then
         f1 = 0._ireals
      end if
      f2 = 1._ireals - f1

      Absups = Absups + ( f1 * Lup_bot_s3 + f2 * Lup_top_s3 ) * a3 * eps_s3 * f_final_s3
      Emups  = Emups - B * a3 * eps_s3 * f_final_s3

      ! Calculate total heating rate from single contributions
      hr = ( Absup + Emup + Absdn + Emdn &
           + ( Absups + Emups+ Absdns + Emdns ) / 2._ireals ) / v * pi

     
      
#ifndef _XLF
      if(isnan(hr)) print *, 'nca shows nan', hr
#endif

      end associate
    end subroutine plexrt_nca




  ! ################################################################################
  ! ################## Function - interpolate emissivity ###########################
  ! # This function interpolates in 3D space between different emissivities ########
  ! ################################################################################

  function interpol_emis (tauxx, tauzz, eps_tab) result(emis)

    real(ireals), intent(in) :: tauxx, tauzz
    real(ireals), intent(in) :: eps_tab(:,:)
    integer(iintegers) :: ix, iy, i, ntau
    real(ireals) :: emis, tauhx, tauz
    real(ireals) :: f1,  f2

    tauhx=tauxx
    tauz =tauzz

    ntau = size(eps_tab, dim=1)

    ! find indeces of tauhx and tauz
    ! at the table limits: set to lower or upper boundary
    do i=1,ntau
       if ((tauhx).gt.tau_hx(i) .and. (tauhx).lt.tau_hx(ntau) .and. (tauhx).gt.tau_hx(1)) then
          ix = i
       else if((tauxx).le.tau_hx(1)) then
          ix=1
       else if ((tauxx).gt.tau_hx(ntau))then
          ix = ntau
          tauhx=tau_hx(ntau)
       end if
    enddo

    do i=1,ntau
       if ((tauz).gt.tau_z(i) .and. (tauz).lt.tau_z(ntau) .and. (tauz).gt.tau_z(1)) then
          iy = i
       else if((tauz).le.tau_z(1)) then
          iy=1
       else if ((tauz).gt.tau_z(ntau))then
          iy = ntau
          tauz=tau_z(ntau)
       end if
    enddo

   !  if optical depth is lower than lookup table, set to 1-exp(-tau) and exit program
    if(tauxx.lt.tau_hx(1) .or. tauzz.lt.tau_z(1)) then
       if(tauzz.lt.tauxx)then
          emis=1-(exp(-tauzz))
       else
          emis=1-(exp(-tauxx))
       end if
    else if(tauxx.lt.tau_hx(ntau).and. tauxx.gt.tau_hx(1))then
       if(tauzz.lt.tau_z(ntau).and.tauzz.gt.tau_z(1))then
          f1=(tau_hx(ix+1)-tauhx)/(tau_hx(ix+1)-tau_hx(ix)) *eps_tab(ix,iy) + (tauhx-tau_hx(ix)) &
               / (tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix+1,iy)
          f2=(tau_hx(ix+1)-tauhx)/(tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix,iy+1) + (tauhx-tau_hx(ix)) &
               / (tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix+1,iy+1)
          emis = (tau_z(iy+1)-tauz)/(tau_z(iy+1)-tau_z(iy)) *f1 + (tauz-tau_z(iy)) &
               / (tau_z(iy+1)-tau_z(iy))*f2
       else if(tauzz.ge.tau_z(ntau)) then
          f1=(tau_hx(ix+1)-tauhx)/(tau_hx(ix+1)-tau_hx(ix)) *eps_tab(ix,iy) + (tauhx-tau_hx(ix)) &
               / (tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix+1,iy)
          f2=(tau_hx(ix+1)-tauhx)/(tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix,iy) + (tauhx-tau_hx(ix)) &
               / (tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix+1,iy)
          emis = (tau_hx(ix+1)-tauhx)/(tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix,iy) + (tauhx-tau_hx(ix)) &
               / (tau_hx(ix+1)-tau_hx(ix))*eps_tab(ix+1,iy)
       end if
    else if(tauxx.ge.tau_hx(ntau-1))then
       f1=1
       f2=1
       if(tauzz.lt.tau_z(ntau).and.tauzz.gt.tau_z(1))then
          emis = (tau_z(iy+1)-tauz)/(tau_z(iy+1)-tau_z(iy)) * eps_tab(ix,iy)+ (tauz-tau_z(iy)) &
               / (tau_z(iy+1)-tau_z(iy))*eps_tab(ix,iy+1)
       else if(tauzz.ge.tau_z(ntau-1))then
          emis=eps_tab(ix,iy)
       end if
    end if

    !correct monte carlo noise in lookup table for high optical thickness (must be 1 in the limit)
    emis = min(emis,real(1._ireals))
  end function interpol_emis

  ! ################################################################################
  ! ################## Function - interpolate emissivity ###########################
  ! # This function interpolates in 3D space between different emissivities ########
  ! ################################################################################

  function interpol_2d (var1, var2, tab1) result(res)

    real(ireals), intent(in) :: var1, var2
    real(ireals), intent(in) :: tab1(:,:)
    integer(iintegers) :: ix, iy, i, n1, n2
    real(ireals) :: res
    real(ireals) :: f1,  f2, tmp1
    n1=size(tab1,dim=1)
    n2=size(tab1,dim=2)

    f1=0
    f2=0
    res=0

    ! find indeces of var1 and var2
    ! at the table limits: set to lower or upper boundary
    do i=1,n1,1
       if ((var1).gt.var_1(i) .and. (var1).lt.var_1(n1) .and. (var1).gt.var_1(1)) then
          ix = i
       else if((var1).le.var_1(1)) then
          ix=1
       else if ((var1).ge.var_1(n1))then
          ix = n1
       end if
    enddo

    do i=1,n2,1
       if ((var2).gt.var_2(i) .and. (var2).lt.var_2(n2) .and. (var2).gt.var_2(1)) then
          iy = i
       else if((var2).le.var_2(1)) then
          iy=1
       else if ((var2).ge.var_2(n2))then
          iy = n2

       end if
    enddo

    !if optical depth is lower than lookup table, set to 1-exp(-tau) and exit program
    if(var1.lt.var_1(n1) .and. var1.ge.var_1(1)) then
       if(var2.lt.var_2(n2).and.var2.ge.var_2(1))then
          f1=(var_1(ix+1)-var1)/(var_1(ix+1)-var_1(ix)) *tab1(ix,iy) + (var1-var_1(ix)) &
               / (var_1(ix+1)-var_1(ix))*tab1(ix+1,iy)
          f2=(var_1(ix+1)-var1)/(var_1(ix+1)-var_1(ix))*tab1(ix,iy+1) + (var1-var_1(ix)) &
               / (var_1(ix+1)-var_1(ix))*tab1(ix+1,iy+1)
          res = (var_2(iy+1)-var2)/(var_2(iy+1)-var_2(iy)) *f1 + (var2-var_2(iy)) &
               / (var_2(iy+1)-var_2(iy))*f2
       else if(var2.ge.var_2(n2)) then
          f1=(var_1(ix+1)-var1)/(var_1(ix+1)-var_1(ix)) *tab1(ix,iy) + (var1-var_1(ix)) &
               / (var_1(ix+1)-var_1(ix))*tab1(ix+1,iy)
          f2=(var_1(ix+1)-var1)/(var_1(ix+1)-var_1(ix))*tab1(ix,iy) + (var1-var_1(ix)) &
               / (var_1(ix+1)-var_1(ix))*tab1(ix+1,iy)
          res=(var_1(ix+1)-var1)/(var_1(ix+1)-var_1(ix))*tab1(ix,iy) + (var1-var_1(ix)) &
               / (var_1(ix+1)-var_1(ix))*tab1(ix+1,iy)
       else if (var2.lt.var_2(1)) then
          tmp1= tab1(ix, iy)
          res=tmp1
       end if

    else if (var1.ge.var_1(n1))then
       if(var2.lt.var_2(n2).and.var2.gt.var_2(1))then
          res=(var_2(iy+1)-var2)/(var_2(iy+1)-var_2(iy)) * tab1(ix,iy) + (var2-var_2(iy)) &
               / (var_2(iy+1)-var_2(iy))* tab1(ix,iy+1)
       else if(var2.ge.var_2(n2)) then
          tmp1=tab1(n1,n2)
          res=tmp1
       else
          tmp1=tab1(1,1)
          res=tmp1
       end if
    end if

  end function interpol_2d

 end module m_plexrt_nca
