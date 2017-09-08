!-------------------------------------------------------------------------
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
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

!> @brief Module contains Raytracer to compute Transfer Coefficients for different 'stream' realizations
!> @author Fabian Jakub LMU/MIM

module m_boxmc
#if defined(__INTEL_COMPILER)
  use ifport
#endif
#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

  use m_helper_functions_dp, only : approx, mean, rmse, imp_reduce_sum, norm, deg2rad, compute_normal_3d, hit_plane, spherical_2_cartesian
  use m_helper_functions, only : CHKERR
  use iso_c_binding
  use m_mersenne
  use mpi
  use m_data_parameters, only: mpiint,iintegers,ireals,ireal_dp,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10, inil, pi_dp

  use m_optprop_parameters, only : delta_scale_truncate,stddev_atol,stddev_rtol,ldebug_optprop

  implicit none

  private
  public :: t_boxmc, t_boxmc_8_10, t_boxmc_1_2, t_boxmc_3_10, t_boxmc_wedge_5_5

  integer,parameter :: fg=1,bg=2,tot=3
  real(ireal_dp),parameter :: zero=0, one=1 ,nil=-9999

  integer(mpiint) :: mpierr,myid,numnodes

  logical :: lRNGseeded=.False.
  type(randomNumberSequence),save :: rndSeq

  ! ******************** TYPE DEFINITIONS ************************
  type,abstract :: t_boxmc
    integer(iintegers) :: dir_streams=inil,diff_streams=inil
    logical :: initialized=.False.
  contains
    procedure :: init
    procedure :: get_coeff
    procedure :: move_photon
    procedure(intersect_distance),deferred :: intersect_distance

    procedure(init_dir_photon),deferred  :: init_dir_photon
    procedure(init_diff_photon),deferred :: init_diff_photon
    procedure(update_dir_stream),deferred  :: update_dir_stream
    procedure(update_diff_stream),deferred :: update_diff_stream
  end type

  type,extends(t_boxmc) :: t_boxmc_1_2
  contains
    procedure :: intersect_distance => intersect_distance_1_2
    procedure :: init_dir_photon    => init_dir_photon_1_2
    procedure :: init_diff_photon   => init_diff_photon_1_2
    procedure :: update_dir_stream  => update_dir_stream_1_2
    procedure :: update_diff_stream => update_diff_stream_1_2
  end type t_boxmc_1_2

  type,extends(t_boxmc) :: t_boxmc_8_10
  contains
    procedure :: intersect_distance => intersect_distance_8_10
    procedure :: init_dir_photon    => init_dir_photon_8_10
    procedure :: init_diff_photon   => init_diff_photon_8_10
    procedure :: update_dir_stream  => update_dir_stream_8_10
    procedure :: update_diff_stream => update_diff_stream_8_10
  end type t_boxmc_8_10

  type,extends(t_boxmc) :: t_boxmc_3_10
  contains
    procedure :: intersect_distance => intersect_distance_3_10
    procedure :: init_dir_photon    => init_dir_photon_3_10
    procedure :: init_diff_photon   => init_diff_photon_3_10
    procedure :: update_dir_stream  => update_dir_stream_3_10
    procedure :: update_diff_stream => update_diff_stream_3_10
  end type t_boxmc_3_10

  type,extends(t_boxmc) :: t_boxmc_wedge_5_5
  contains
    procedure :: intersect_distance => intersect_distance_wedge_5_5
    procedure :: init_dir_photon    => init_dir_photon_wedge_5_5
    procedure :: init_diff_photon   => init_diff_photon_wedge_5_5
    procedure :: update_dir_stream  => update_dir_stream_wedge_5_5
    procedure :: update_diff_stream => update_diff_stream_wedge_5_5
  end type t_boxmc_wedge_5_5

  type photon
    real(ireal_dp) :: loc(3)=nil,dir(3)=nil,weight=nil,dx=nil,dy=nil,dz=nil
    logical :: alive=.True.,direct=.False.
    integer(iintegers) :: side=inil,src=inil,scattercnt=0
    real(ireal_dp) :: optprop(3) ! kabs,ksca,g
  end type

  type stddev
    real(ireal_dp),allocatable,dimension(:) :: inc,delta,mean,mean2,var,relvar
    logical :: converged=.False.
    real(ireal_dp) :: atol=zero, rtol=zero
  end type
  ! ******************** TYPE DEFINITIONS ************************


  ! ***************** INTERFACES ************
  abstract interface
    subroutine init_diff_photon(bmc,p,src,dx,dy,dz)
      import :: t_boxmc,photon,iintegers,ireal_dp
      class(t_boxmc) :: bmc
      type(photon),intent(inout) :: p
      real(ireal_dp),intent(in) :: dx,dy,dz
      integer(iintegers),intent(in) :: src
    end subroutine
  end interface

  abstract interface
    subroutine init_dir_photon(bmc,p,src,direct,initial_dir,dx,dy,dz)
      import :: t_boxmc,photon,iintegers,ireal_dp
      class(t_boxmc) :: bmc
      type(photon),intent(inout) :: p
      real(ireal_dp),intent(in) :: dx,dy,dz,initial_dir(3)
      integer(iintegers),intent(in) :: src
      logical,intent(in) :: direct
    end subroutine
  end interface

  abstract interface
    subroutine update_diff_stream(bmc,p,S)
      import :: t_boxmc,photon,iintegers,ireal_dp
      class(t_boxmc) :: bmc
      type(photon),intent(in) :: p
      real(ireal_dp),intent(inout) :: S(:)
    end subroutine
  end interface

  abstract interface
    subroutine update_dir_stream(bmc,p,T)
      import :: t_boxmc,photon,iintegers,ireal_dp
      class(t_boxmc) :: bmc
      type(photon),intent(in) :: p
      real(ireal_dp),intent(inout) :: T(:)
    end subroutine
  end interface

  abstract interface
    subroutine intersect_distance(bmc,p,max_dist)
      import :: t_boxmc,photon,ireal_dp
      class(t_boxmc) :: bmc
      type(photon),intent(inout) :: p
      real(ireal_dp),intent(out) :: max_dist
    end subroutine
  end interface
! ***************** INTERFACES ************

contains


  !> @brief Calculate Transfer Coefficients using MonteCarlo integration
  !> @details All MPI Nodes start photons from src stream and ray trace it including scattering events through the box until it leaves the box through one of the exit streams.\n
  !> Scattering Absorption is accounted for by carrying on a photon weight and succinctly lower it by lambert Beers Law \f$ \omega_{abso}^{'} = \omega_{abso} \cdot e^{- \rm{d}s \cdot {\rm k}_{sca}   }   \f$ \n
  !> New Photons are started until we reach a stdvariance which is lower than the given stddev in function call init_stddev. Once this precision is reached, we exit the photon loop and build the average with all the other MPI Nodes.
  subroutine get_coeff(bmc,comm,op_bg,src,ldir,phi0,theta0,dx,dy,dz, ret_S_out, ret_T_out, ret_S_tol,ret_T_tol, inp_atol, inp_rtol )
    class(t_boxmc)                :: bmc                          !< @param[in] bmc Raytracer Type - determines number of streams
    real(ireals),intent(in)       :: op_bg(3)                     !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
    real(ireals),intent(in)       :: phi0                         !< @param[in] phi0 solar azimuth angle
    real(ireals),intent(in)       :: theta0                       !< @param[in] theta0 solar zenith angle
    integer(iintegers),intent(in) :: src                          !< @param[in] src stream from which to start photons - see init_photon routines
    integer(mpiint),intent(in)    :: comm                         !< @param[in] comm MPI Communicator
    logical,intent(in)            :: ldir                         !< @param[in] ldir determines if photons should be started with a fixed incidence angle
    real(ireals),intent(in)       :: dx,dy,dz                     !< @param[in] dx,dy,dz box with dimensions in [m]
    real(ireals),intent(out)      :: ret_S_out(bmc%diff_streams)  !< @param[out] S_out diffuse streams transfer coefficients
    real(ireals),intent(out)      :: ret_T_out(bmc%dir_streams)   !< @param[out] T_out direct streams transfer coefficients
    real(ireals),intent(out)      :: ret_S_tol(bmc%diff_streams)  !< @param[out] absolute tolerances of results
    real(ireals),intent(out)      :: ret_T_tol(bmc%dir_streams)   !< @param[out] absolute tolerances of results
    real(ireals),intent(in),optional :: inp_atol                  !< @param[in] inp_atol if given, determines targeted absolute stddeviation
    real(ireals),intent(in),optional :: inp_rtol                  !< @param[in] inp_rtol if given, determines targeted relative stddeviation


    real(ireal_dp) :: S_out(bmc%diff_streams)
    real(ireal_dp) :: T_out(bmc%dir_streams)
    real(ireal_dp) :: S_tol(bmc%diff_streams)
    real(ireal_dp) :: T_tol(bmc%dir_streams)

    real(ireal_dp) :: atol,rtol, coeffnorm

    type(stddev) :: std_Sdir, std_Sdiff, std_abso

    integer(iintegers) :: Nphotons


    if(.not. bmc%initialized ) stop 'Box Monte Carlo Ray Tracer is not initialized! - This should not happen!'

    if(present(inp_atol)) then
      atol = inp_atol
    else
      atol = stddev_atol
    endif
    if(present(inp_rtol)) then
      rtol = inp_rtol
    else
      rtol = stddev_rtol
    endif

    call init_stddev( std_Sdir , bmc%dir_streams  ,atol, rtol )
    call init_stddev( std_Sdiff, bmc%diff_streams ,atol, rtol )
    call init_stddev( std_abso , i1               ,atol, rtol )

    if(.not.ldir) std_Sdir%converged=.True.

    if( (any(op_bg.lt.zero)) .or. (any(isnan(op_bg))) ) then
      print *,'corrupt optical properties: bg:: ',op_bg
      call exit
    endif

    if( any([phi0,theta0].lt.0) .or. (phi0.gt.360_ireals) .or. (theta0.gt.180_ireals) ) then
      print *,'corrupt sun angles :: ',phi0, theta0
      call exit
    endif

    if(dx.le.zero .or. dy.le.zero .or. dz.le.zero ) then
      print *,'ERROR: box dimensions have to be positive!',dx,dy,dz
      call exit()
    endif

    call run_photons(bmc,src,                   &
                     real(op_bg,kind=ireal_dp), &
                     real(dx, kind=ireal_dp),   &
                     real(dy, kind=ireal_dp),   &
                     real(dz, kind=ireal_dp),   &
                     ldir,                      &
                     real(phi0,   kind=ireal_dp), &
                     real(theta0, kind=ireal_dp), &
                     Nphotons, std_Sdir,std_Sdiff,std_abso)


    S_out = std_Sdiff%mean
    T_out = std_Sdir%mean

    ! tolerances that we achieved and report them back
    S_tol = std_Sdiff%var
    T_tol = std_Sdir%var

    if(numnodes.gt.1) then ! average reduce results from all ranks
      call reduce_output(Nphotons, comm, S_out, T_out, S_tol, T_tol)
    endif


    ! some debug output at the end...
    coeffnorm = sum(S_out)+sum(T_out)
    if( coeffnorm.gt.one ) then
      if(coeffnorm.ge.one+1e-5_ireal_dp) then
        print *,'ohoh something is wrong! - sum of streams is bigger 1, this cant be due to energy conservation',&
        sum(S_out),'+',sum(T_out),'=',sum(S_out)+sum(T_out),'.gt',one,':: op',op_bg,'eps',epsilon(one)
        call exit
      else
        S_out = S_out / (coeffnorm+epsilon(coeffnorm)*10)
        T_out = T_out / (coeffnorm+epsilon(coeffnorm)*10)
        if(ldebug_optprop) print *,'renormalizing coefficients :: ',coeffnorm,' => ',sum(S_out)+sum(T_out)
      endif
      if( (sum(S_out)+sum(T_out)).gt.one ) then
        print *,'norm still too big',sum(S_out)+sum(T_out)
        call exit
      endif
    endif
    if( (any(isnan(S_out) )) .or. (any(isnan(T_out)) ) ) then
      print *,'Found a NaN in output! this should not happen! dir',T_out,'diff',S_out
      print *,'Input:', op_bg, '::', phi0, theta0, src, ldir, '::', dx,dy,dz
      call exit()
    endif

    ret_S_out = real(S_out, kind=ireals)
    ret_T_out = real(T_out, kind=ireals)
    ret_S_tol = real(S_tol, kind=ireals)
    ret_T_tol = real(T_tol, kind=ireals)

  end subroutine

  subroutine run_photons(bmc,src,op,dx,dy,dz,ldir,phi0,theta0,Nphotons,std_Sdir,std_Sdiff,std_abso)
      class(t_boxmc),intent(inout) :: bmc
      integer(iintegers),intent(in) :: src
      real(ireal_dp),intent(in) :: op(3),dx,dy,dz,phi0,theta0
      logical,intent(in) :: ldir
      integer(iintegers) :: Nphotons
      type(stddev),intent(inout)   :: std_Sdir, std_Sdiff, std_abso

      type(photon)       :: p
      integer(iintegers) :: k,mycnt,mincnt
      real(ireal_dp)   :: initial_dir(3)
      real(ireal_dp)   :: time(2)


      call cpu_time(time(1))

      ! we turn the initial direction in x and y, against the convetion of sun angles...
      ! i.e. here we have azimuth phi = 0, beam going towards the north
      ! and phi = 90, beam going towards east
      initial_dir = spherical_2_cartesian(phi0, theta0) * [-one, -one, one]

      mincnt= max( 1000, int( 1e3 /numnodes ) )
      mycnt = int(1e8)/numnodes
      mycnt = min( max(mincnt, mycnt ), huge(k)-1 )
      do k=1, mycnt

          if(k.gt.mincnt .and. all([std_Sdir%converged, std_Sdiff%converged, std_abso%converged ]) ) exit

          if(ldir) then
              call bmc%init_dir_photon(p,src,ldir,initial_dir,real(dx,kind=ireal_dp),real(dy,kind=ireal_dp),real(dz,kind=ireal_dp))
          else
              call bmc%init_diff_photon(p,src,real(dx,kind=ireal_dp),real(dy,kind=ireal_dp),real(dz,kind=ireal_dp))
          endif
          p%optprop = op

          move: do
              call bmc%move_photon(p)
              call roulette(p)

              if(.not.p%alive) exit move
              call scatter_photon(p)
          enddo move

          if(ldir) call refill_direct_stream(p,initial_dir)

          std_abso%inc = one-p%weight
          std_Sdir%inc  = zero
          std_Sdiff%inc = zero

          if(p%direct) then
              call bmc%update_dir_stream(p,std_Sdir%inc)
          else
              call bmc%update_diff_stream(p,std_Sdiff%inc)
          endif

          if (ldir) call std_update( std_Sdir , k, i1*numnodes )
          call std_update( std_abso , k, i1*numnodes)
          call std_update( std_Sdiff, k, i1*numnodes )
      enddo ! k photons
      Nphotons = k

      call cpu_time(time(2))

      !if(rand().gt..99_ireal_dp) then
      !  write(*,FMT='("src ",I0," dz",I0," op ",3(ES12.3),"(delta",3(ES12.3),") sun(",I0,",",I0,") N_phot ",I0 ,"=>",ES12.3,"phot/sec/node took",ES12.3,"sec")') &
      !    src,int(dz),op,p%optprop,int(phi0),int(theta0),Nphotons, Nphotons/max(epsilon(time),time(2)-time(1))/numnodes,time(2)-time(1)
      !endif
  end subroutine

  !> @brief take weighted average over mpi processes
  subroutine reduce_output(Nlocal, comm, S_out, T_out, S_tol, T_tol)
    integer(iintegers),intent(in) :: Nlocal
    integer(mpiint),intent(in)    :: comm
    real(ireal_dp),intent(inout)      :: S_out(:)
    real(ireal_dp),intent(inout)      :: T_out(:)
    real(ireal_dp),intent(inout)      :: S_tol(:)
    real(ireal_dp),intent(inout)      :: T_tol(:)

    real(ireal_dp) :: Nglobal
    ! weight mean by calculated photons and compare it with results from other nodes
    Nglobal = Nlocal
    call imp_reduce_sum(comm, Nglobal, myid)

    call reduce_var(comm, Nlocal, Nglobal, S_out)
    call reduce_var(comm, Nlocal, Nglobal, T_out)
    !TODO: combining stddeviation is probably not just the arithmetic mean?
    call reduce_var(comm, Nlocal, Nglobal, S_tol)
    call reduce_var(comm, Nlocal, Nglobal, T_tol)

  contains
    subroutine reduce_var(comm, Nlocal, Nglobal, arr)
      integer(iintegers),intent(in) :: Nlocal
      real(ireal_dp),intent(in) :: Nglobal
      real(ireal_dp),intent(inout) :: arr(:)
      integer(mpiint),intent(in) :: comm
      integer(iintegers) :: k

      arr = arr*Nlocal
      do k=1,size(arr)
        call imp_reduce_sum(comm, arr(k), myid)
      enddo
      arr = arr/Nglobal
    end subroutine
  end subroutine

  !> @brief russian roulette helps to reduce computations with not much weight
  subroutine roulette(p)
    type(photon),intent(inout) :: p
    real(ireal_dp),parameter :: m=1e-3_ireal_dp,s=1e-6_ireal_dp*m

    if(p%weight.lt.s) then
      if(R().ge.p%weight/m) then
        p%weight=zero
        p%alive=.False.
      else
        p%weight=m
      endif
    endif
  end subroutine

  !> @brief in a ``postprocessing`` step put scattered direct radiation back into dir2dir streams
  subroutine refill_direct_stream(p,initial_dir)
    type(photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: initial_dir(3)

    real(ireal_dp) :: angle

    angle = dot_product(initial_dir, p%dir)

    if(angle.gt.delta_scale_truncate) then
      p%direct = .True.
      !          print *,'delta scaling photon initial', initial_dir,'dir',p%dir,'angle',angle,'cos', (acos(angle))*180/pi_dp
    endif
  end subroutine

  !> @brief return equally distributed random number in [-1,1]
  function s()
    real(ireal_dp) :: s
    s = -one + 2*R()
  end function
  !> @brief return cosine of deg(in degrees)
  function deg2mu(deg)
    real(ireal_dp),intent(in) :: deg
    real(ireal_dp) :: deg2mu
    deg2mu = cos(deg2rad(deg))
  end function
  !> @brief return uniform random number between a and b
  function interv_R(a,b)
    real(ireal_dp),intent(in) :: a,b
    real(ireal_dp) :: interv_R
    real(ireal_dp) :: lb,ub

    lb=min(a,b)
    ub=max(a,b)
    interv_R = lb + R()*(ub-lb)
  end function
  !> @brief return uniform random number between [0,v] with a certain cutoff at the edges to make sure that photons are started ``inside`` a box
  function L(v)
    real(ireal_dp) :: L
    real(ireal_dp),intent(in) ::v
    real(ireal_dp),parameter :: eps=epsilon(L)*1e3_ireal_dp
    L = ( eps + R()*(one-2*eps) ) *v
  end function

  !> @brief main function for a single photon
  !> @details this routine will incrementally move a photon until it is either out of the domain or it is time for a interaction with the medium
  subroutine move_photon(bmc,p)
    class(t_boxmc) :: bmc
    type(photon),intent(inout) :: p
    real(ireal_dp) :: dist,intersec_dist

    call bmc%intersect_distance(p,intersec_dist)

    dist = distance(tau(R() ), get_ksca(p) )

    if(intersec_dist .le. dist) then
      p%alive=.False.
      call update_photon_loc(p,intersec_dist)
    else
      call update_photon_loc(p,dist)
    endif

    if(p%scattercnt.gt.1e9) then
      print *,'Scattercnt:',p%scattercnt,' -- maybe this photon got stuck? -- I will move this one out of the box but keep in mind, that this is a dirty hack i.e. absorption will be wrong!'
      call print_photon(p)
      p%alive=.False.
      call update_photon_loc(p,intersec_dist)
      call print_photon(p)
    endif

  contains

    !> @brief update physical location of photon and consider absorption
    subroutine update_photon_loc(p,dist)
      type(photon),intent(inout) :: p
      real(ireal_dp),intent(in) :: dist
      call absorb_photon(p,dist)
      p%loc = p%loc + (dist*p%dir)
      if(any(isnan(p%loc))) then
        print *,'loc is now a NAN! ',p%loc,'dist',dist
        call print_photon(p)
        call exit
      endif
    end subroutine
    !> @brief compute physical distance according to travel_tau
    elemental function distance(tau,beta)
      real(ireal_dp),intent(in) :: tau,beta
      real(ireal_dp) :: distance
      if(approx(beta,zero) ) then
          distance = huge(distance)
      else
          distance = tau/beta
      endif
    end function

    !> @brief throw the dice for a random optical thickness -- after the corresponding dtau it is time to do some interaction
    elemental function tau(r)
      real(ireal_dp),intent(in) :: r
      real(ireal_dp) :: tau,arg
      arg = max( epsilon(arg), one-r )
      tau = -log(arg)
    end function
  end subroutine move_photon

  !> @brief cumulative sum of henyey greenstein phase function
  elemental function hengreen(r,g)
    real(ireal_dp),intent(in) :: r,g
    real(ireal_dp) :: hengreen
    real(ireal_dp),parameter :: two=2*one
    if( approx(g,zero) ) then
      hengreen = two*r-one
    else
      hengreen = one/(two*g) * (one+g**two - ( (one-g**two) / ( two*g*r + one-g) )**two )
    endif
    hengreen = min(max(hengreen,-one), one)
  end function

  !> @brief remove photon weight due to absorption
  pure subroutine absorb_photon(p,dist)
    type(photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: dist
    real(ireal_dp) :: new_weight,tau

    tau = get_kabs(p)*dist
    if(tau.gt.30) then
      p%weight = zero
    else
      new_weight = p%weight * exp(- tau)
      p%weight = new_weight
    endif
  end subroutine

  !> @brief compute new direction of photon after scattering event
  subroutine scatter_photon(p)
    type(photon),intent(inout) :: p
    real(ireal_dp) :: muxs,muys,muzs,muxd,muyd,muzd
    real(ireal_dp) :: mutheta,fi,costheta,sintheta,sinfi,cosfi,denom,muzcosfi

    mutheta = hengreen(R(),get_g(p))

    p%scattercnt = p%scattercnt+1
    p%direct=.False.

    muxs = p%dir(1)
    muys = p%dir(2)
    muzs = p%dir(3)

    fi = R()*pi_dp*2

    costheta = (mutheta)
    sintheta = sqrt(one-costheta**2)
    sinfi = sin(fi)
    cosfi = cos(fi)

    if( approx(muzs , one) ) then
      muxd = sintheta*cosfi
      muyd = sintheta*sinfi
      muzd = costheta
    else if ( approx( muzs ,-one) ) then
      muxd =  sintheta*cosfi
      muyd = -sintheta*sinfi
      muzd = -costheta
    else
      denom = sqrt(one-muzs**2)
      muzcosfi = muzs*cosfi

      muxd = sintheta*(muxs*muzcosfi-muys*sinfi)/denom + muxs*costheta
      muyd = sintheta*(muys*muzcosfi+muxs*sinfi)/denom + muys*costheta
      muzd = -denom*sintheta*cosfi + muzs*costheta
    endif

    !        if(isnan(muxd).or.isnan(muyd).or.isnan(muzd) ) print *,'new direction is NAN :( --',muxs,muys,muzs,mutheta,fi,'::',muxd,muyd,muzd

    p%dir(1) = muxd
    p%dir(2) = muyd
    p%dir(3) = muzd

  end subroutine

  pure function get_kabs(p)
    real(ireal_dp) :: get_kabs
    type(photon),intent(in) :: p
    get_kabs = p%optprop(1)
  end function
  pure function get_ksca(p)
    real(ireal_dp) :: get_ksca
    type(photon),intent(in) :: p
    get_ksca = p%optprop(2)
  end function
  pure function get_g(p)
    real(ireal_dp) :: get_g
    type(photon),intent(in) :: p
    get_g = p%optprop(3)
  end function

  subroutine print_photon(p)
    type(photon),intent(in) :: p
    print *,'---------------------------'
    print *,'Location  of Photon:',p%loc
    print *,'Direction of Photon:',p%dir
    print *,'weight',p%weight,'alive,direct',p%alive,p%direct,'scatter count',p%scattercnt
    print *,'side',p%side,'src',p%src
    print *,'kabs,ksca,g',get_kabs(p),get_ksca(p),get_g(p)
    print *,'dx,dy,dz', p%dx, p%dy, p%dz
  end subroutine

  function R()
    real(ireal_dp) :: R
    real :: rvec(1)
    ! R = getRandomDouble(rndSeq) ! mersenne twister from robert pinucs, see mersenne.f90
    ! call random_number(R)
    call RANLUX(rvec,1)  ! use Luxury Pseudorandom Numbers from M. Luscher
    R = real(rvec(1), kind=ireal_dp)
  end function

  subroutine init_random_seed(myid)
    integer,intent(in) :: myid
    INTEGER :: i, n, clock, s
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    real :: rn

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = myid*3*7*11*13*17 + clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)

    call random_number(rn)
    s = int(rn*1000)*(myid+1)

    !  s=myid
    !  print *,myid,'Seeding RNG with ',s
    rndSeq = new_RandomNumberSequence(seed=s) ! seed pincus's mersenne twister
    call RLUXGO(4, int(s), 0, 0) ! seed ranlux rng
    lRNGseeded=.True.
  end subroutine
  subroutine init_stddev( std, N, atol, rtol)
    type(stddev),intent(inout) :: std
    real(ireal_dp),intent(in) :: atol,rtol
    integer(iintegers) :: N
    if( allocated(std%inc     ) ) deallocate( std%inc  )
    if( allocated(std%delta   ) ) deallocate( std%delta)
    if( allocated(std%mean    ) ) deallocate( std%mean )
    if( allocated(std%mean2   ) ) deallocate( std%mean2)
    if( allocated(std%var     ) ) deallocate( std%var  )
    if( allocated(std%relvar  ) ) deallocate( std%relvar)
    allocate( std%inc   (N)) ; std%inc   = zero
    allocate( std%delta (N)) ; std%delta = zero
    allocate( std%mean  (N)) ; std%mean  = zero
    allocate( std%mean2 (N)) ; std%mean2 = zero
    allocate( std%var   (N)) ; std%var   = zero
    allocate( std%relvar(N)) ; std%relvar= zero
    std%atol = atol
    std%rtol = rtol
    std%converged = .False.
  end subroutine

  pure subroutine std_update(std, N, numnodes)
    type(stddev),intent(inout) :: std
    integer(iintegers),intent(in) :: N, numnodes
    real(ireal_dp),parameter :: relvar_limit=1e-6_ireal_dp

    std%delta = std%inc   - std%mean
    std%mean  = std%mean  + std%delta/N
    std%mean2 = std%mean2 + std%delta * ( std%inc - std%mean )
    std%var = sqrt( std%mean2/N ) / sqrt( one*N*numnodes )
    where(std%mean.gt.relvar_limit)
      std%relvar = std%var / std%mean
    elsewhere
      std%relvar = .1_ireal_dp/sqrt(one*N) ! consider adding a photon weight of .1 as worst case that could happen for the next update...
    end where

    if( all( (std%var .lt. std%atol) .and. (std%relvar .lt. std%rtol) ) ) then
      std%converged = .True.
    else
      std%converged = .False.
    endif
  end subroutine

  subroutine init(bmc, comm)
    class(t_boxmc) :: bmc
    integer(mpiint),intent(in) :: comm

    !        print *,'initializing boxmc'
    if(comm.eq.-1) then
      myid = -1
    else
      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    endif

    if(.not.lRNGseeded) call init_random_seed(myid+2)

    numnodes=1
    if(myid.ge.0) call mpi_comm_size(comm,numnodes,mpierr); call CHKERR(mpierr)

    select type (bmc)
    class is (t_boxmc_8_10)
    bmc%dir_streams  =  8
    bmc%diff_streams = 10
    class is (t_boxmc_3_10)
    bmc%dir_streams  =  3
    bmc%diff_streams = 10
    class is (t_boxmc_1_2)
    bmc%dir_streams  =  1
    bmc%diff_streams =  2
    class is (t_boxmc_wedge_5_5)
    bmc%dir_streams  =  5
    bmc%diff_streams =  5
    class default
    stop 'initialize: unexpected type for boxmc object!'
  end select

  bmc%initialized = .True.
end subroutine


include 'boxmc_8_10.inc'
include 'boxmc_3_10.inc'
include 'boxmc_1_2.inc'
include 'boxmc_wedge_5_5.inc'

end module
