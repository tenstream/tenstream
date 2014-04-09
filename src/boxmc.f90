!> @brief Module contains Raytracer to compute Transfer Coefficients for different 'stream' realizations
!> @author Fabian Jakub LMU/MIM

module boxmc
      use helper_functions, only : approx,mean,rmse,deg2rad,norm
      use iso_c_binding
      use mersenne
      use mpi
      use data_parameters, only: mpiint,iintegers,ireals,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10, zero,one,nil,inil,pi
      
      use boxmc_parameters, only : delta_scale_truncate

      implicit none

      private
      public :: t_boxmc,t_boxmc_8_10,t_boxmc_1_2

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

      type photon
              double precision :: loc(3)=nil,dir(3)=nil,weight=nil,dx=nil,dy=nil,dz=nil
              logical :: alive=.True.,direct=.False.
              integer(iintegers) :: side=inil,src=inil,scattercnt=0
              double precision :: optprop(3) ! kabs,ksca,g
      end type

      type stddev
        real(ireals),allocatable,dimension(:) :: inc,delta,mean,mean2,var
        logical :: converged=.False.
        real(ireals) :: atol=zero, rtol=zero
      end type

      ! ******************** TYPE DEFINITIONS ************************

      integer,parameter :: fg=1,bg=2,tot=3
      integer(iintegers),parameter :: Nphotons_max=int(1e9)

      type(stddev) :: std_Sdir, std_Sdiff, std_abso

      integer :: ierr,myid,numnodes

      logical :: lRNGseeded=.False.
      type(randomNumberSequence) :: rndSeq

      ! ***************** INTERFACES ************
      abstract interface
        subroutine init_diff_photon(bmc,p,src,dx,dy,dz)
          import :: t_boxmc,photon,iintegers
          class(t_boxmc) :: bmc
          type(photon),intent(inout) :: p
          double precision,intent(in) :: dx,dy,dz
          integer(iintegers),intent(in) :: src
        end subroutine
      end interface

      abstract interface
        subroutine init_dir_photon(bmc,p,src,direct,initial_dir,dx,dy,dz)
          import :: t_boxmc,photon,iintegers
          class(t_boxmc) :: bmc
          type(photon),intent(inout) :: p
          double precision,intent(in) :: dx,dy,dz,initial_dir(3)
          integer(iintegers),intent(in) :: src
          logical,intent(in) :: direct
        end subroutine
      end interface

      abstract interface
        subroutine update_diff_stream(bmc,p,S,N)
          import :: t_boxmc,photon,iintegers,ireals
          class(t_boxmc) :: bmc
          type(photon),intent(in) :: p
          real(ireals),intent(inout) :: S(:)
          integer(iintegers),intent(inout) :: N(:)
        end subroutine
      end interface

      abstract interface
        subroutine update_dir_stream(bmc,p,S,N)
          import :: t_boxmc,photon,iintegers,ireals
          class(t_boxmc) :: bmc
          type(photon),intent(in) :: p
          real(ireals),intent(inout) :: S(:)
          integer(iintegers),intent(inout) :: N(:)
        end subroutine
      end interface

      abstract interface
        subroutine intersect_distance(bmc,p,max_dist)
          import :: t_boxmc,photon
          class(t_boxmc) :: bmc
          type(photon),intent(inout) :: p
          double precision,intent(out) :: max_dist
        end subroutine
      end interface
      ! ***************** INTERFACES ************

contains


  !> @brief Calculate Transfer Coefficients using MonteCarlo integration
  !> @details All MPI Nodes start photons from src stream and ray trace it including scattering events through the box until it leaves the box through one of the exit streams.\n
  !> Scattering Absorption is accounted for by carrying on a photon weight and succinctly lower it by lambert Beers Law \f$ \omega_{abso}^{'} = \omega_{abso} \cdot e^{- \rm{d}s \cdot {\rm k}_{sca}   }   \f$ \n
  !> New Photons are started until we reach a stdvariance which is lower than the given stddev in function call init_stddev. Once this precision is reached, we exit the photon loop and build the average with all the other MPI Nodes.
  subroutine get_coeff(bmc,comm,op_bg,src,S_out,Sdir_out,ldir,ldeltascale,phi0,theta0,dx,dy,dz)
      class(t_boxmc)                :: bmc                       !< @param[in] bmc Raytracer Type - determines number of streams
      double precision,intent(in)   :: op_bg(3)                  !< @param[in] op_bg optical properties have to be given as [kabs,ksca,g]
      double precision,intent(in)   :: phi0                      !< @param[in] phi0 solar azimuth angle
      double precision,intent(in)   :: theta0                    !< @param[in] theta0 solar zenith angle
      logical,intent(in)            :: ldeltascale               !< @param[in] ldeltascale implemented as following: if a photon leaves the box with an angle, that is approximately the same as the incidence angle, it is counted as direct
      integer(iintegers),intent(in) :: src                       !< @param[in] src stream from which to start photons - see init_photon routines
      integer(mpiint),intent(in)    :: comm                      !< @param[in] comm MPI Communicator
      logical,intent(in)            :: ldir                      !< @param[in] ldir determines if photons should be started with a fixed incidence angle
      double precision,intent(in)   :: dx,dy,dz                  !< @param[in] dx,dy,dz box with dimensions in [m]
      double precision,intent(out)  :: S_out(bmc%diff_streams)   !< @param[out] S_out diffuse streams transfer coefficients
      double precision,intent(out)  :: Sdir_out(bmc%dir_streams) !< @param[out] Sdir_out direct streams transfer coefficients

      type(photon)       :: p
      integer(iintegers) :: k,mycnt
      double precision   :: time(2),total_photons,initial_dir(3)
      integer(iintegers) :: Ndir(bmc%dir_streams),Ndiff(bmc%diff_streams)

      if(.not. bmc%initialized ) stop 'Box Monte Carlo Ray Tracer is not initialized! - This should not happen!'

      Ndir=i0;Ndiff=i0

      call init_stddev( std_Sdir , bmc%dir_streams  ,1e-5_ireals, 5e-3_ireals )
      call init_stddev( std_Sdiff, bmc%diff_streams ,1e-5_ireals, 5e-3_ireals )
      call init_stddev( std_abso , i1               ,1e-5_ireals, 5e-4_ireals )

      if(.not.ldir) std_Sdir%converged=.True.

      initial_dir  = [ sin(deg2rad(theta0))*sin(deg2rad(phi0)) ,&
                        sin(deg2rad(theta0))*cos(deg2rad(phi0)) ,&
                                          - cos(deg2rad(theta0)) ]
      initial_dir = initial_dir/norm(initial_dir)

      if( (any(op_bg.lt.zero)) .or. (any(isnan(op_bg))) ) then
        print *,'corrupt optical properties: bg:: ',op_bg
        call exit
      endif

      if(dx.le.zero .or. dy.le.zero .or. dz.le.zero ) then
        print *,'ERROR: box dimensions have to be positive!',dx,dy,dz
        call exit()
      endif

      call cpu_time(time(1))

      mycnt = max(i1*numnodes,Nphotons_max/numnodes)
      do k=1,mycnt
        if(k*numnodes.gt.1e4 .and. all([std_Sdir%converged, std_Sdiff%converged, std_abso%converged ]) ) exit

        if(ldir) then
          call bmc%init_dir_photon(p,src,ldir,initial_dir,dx,dy,dz)
        else
          call bmc%init_diff_photon(p,src,dx,dy,dz)
        endif
        p%optprop = op_bg

        move: do
          call bmc%move_photon(p)
          call roulette(p)

          if(.not.p%alive) exit move
          call scatter_photon(p)
        enddo move

        if(ldir) call delta_scaling(p,ldeltascale,initial_dir)

        std_abso%inc = one-p%weight
        std_Sdir%inc  = zero
        std_Sdiff%inc = zero

        if(p%direct) then
          call bmc%update_dir_stream(p,std_Sdir%inc,Ndir)
        else
          call bmc%update_diff_stream(p,std_Sdiff%inc,Ndiff)
        endif

        if(ldir) call std_update( std_Sdir , k, i1*numnodes )
        call std_update( std_abso , k, i1*numnodes)
        call std_update( std_Sdiff, k, i1*numnodes )
      enddo

      total_photons=k
      S_out    = dble( std_Sdiff%mean*k )
      if(ldir) Sdir_out = dble( std_Sdir%mean *k )
      if(myid.ge.0) then
        call mpi_reduce_sum(total_photons,comm,myid)
        do k=1,ubound(S_out,1)
          call mpi_reduce_sum(S_out(k),comm,myid)
        enddo

        if(ldir) then
          do k=1,ubound(Sdir_out,1)
            call mpi_reduce_sum(Sdir_out(k),comm,myid)
          enddo
        endif
      endif
      S_out    = S_out / total_photons
      if(ldir) then
        Sdir_out = Sdir_out / total_photons
      else
        Sdir_out=zero
      endif

      if( (sum(S_out)+sum(Sdir_out)).gt.one+1e-6_ireals) then
        print *,'ohoh something is wrong! - sum of streams is bigger 1, this cant be due to energy conservation',sum(S_out),'+',sum(Sdir_out),'=',sum(S_out)+sum(Sdir_out),':: op',p%optprop
        call print_photon(p)
        call exit
      endif
      if( (any(isnan(S_out) )) .or. (any(isnan(Sdir_out)) ) ) then
        print *,'Found a NaN in output! this should not happen! dir',Sdir_out,'diff',S_out
        call print_photon(p)
        call exit()
      endif
      call cpu_time(time(2))

      if(myid.le.0.and.total_photons.ge.1e7) print *,src,dz,op_bg,'angles',phi0,theta0,'took',time(2)-time(1),'s',' phots*1e3 :: ',total_photons/1e3,' abso :',one-sum(Sdir_out)-sum(S_out),':',total_photons/(time(2)-time(1))/numnodes,'phots/sec/node'
  end subroutine

subroutine mpi_reduce_sum(v,comm,myid)
    double precision,intent(inout) :: v
    integer,intent(in) :: comm,myid
    integer :: ierr
    if(myid.eq.0) then
      call mpi_reduce(MPI_IN_PLACE, v, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
    else
      call mpi_reduce(v, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
    endif
end subroutine

subroutine roulette(p)
        type(photon),intent(inout) :: p
        double precision,parameter :: m=1e-1_ireals,s=1e-3_ireals*m

        if(p%weight.lt.s) then
                if(R().ge.p%weight/m) then
                        p%weight=zero
                        p%alive=.False.
                else
                        p%weight=m
                endif
        endif
end subroutine
subroutine delta_scaling(p,ldeltascale,initial_dir)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: initial_dir(3)
        logical,intent(in) :: ldeltascale

        double precision :: angle
        if((.not.ldeltascale).or.p%direct) return

        angle = dot_product(initial_dir, p%dir)

        if(angle.gt.delta_scale_truncate) then
          p%direct = .True.
!          print *,'delta scaling photon initial', initial_dir,'dir',p%dir,'angle',angle,'cos', (acos(angle))*180/pi
        endif
end subroutine

double precision function s()
        s = -one + 2*R()
end function
function L(v)
    double precision :: L
    double precision,intent(in) ::v
    double precision,parameter :: eps=1e-3
    L = min( max(R()*v,eps), v-eps)
end function

subroutine move_photon(bmc,p)
        class(t_boxmc) :: bmc
        type(photon),intent(inout) :: p
        double precision :: tau_travel,dist,intersec_dist

        tau_travel = tau(R())
        call bmc%intersect_distance(p,intersec_dist)

        dist = distance(tau(R() ), get_ksca(p) )

        if(intersec_dist .le. dist) then
          p%alive=.False.
          call update_photon_loc(p,intersec_dist)
        else
          call update_photon_loc(p,dist)
          return
        endif
end subroutine
subroutine update_photon_loc(p,dist)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dist
        call absorb_photon(p,dist)
        p%loc = p%loc + (dist*p%dir)
        if(any(isnan(p%loc))) then
          print *,'loc is now a NAN! ',p%loc,'dist',dist
          call print_photon(p)
          call exit
        endif
end subroutine
pure double precision function hit_plane(p,po,pn)
        type(photon),intent(in) :: p
        double precision,intent(in) :: po(3),pn(3)
        double precision :: discr
        discr = dot_product(p%dir,pn)
        if( approx(discr, zero) ) then
                hit_plane=huge(hit_plane)
        else        
                hit_plane = dot_product(po-p%loc, pn) / discr
        endif
end function

elemental function distance(tau,beta)
        double precision,intent(in) :: tau,beta
        double precision :: distance
        distance = tau/beta
        if(approx(beta,zero) ) distance=huge(distance)
end function

elemental function tau(r)
        double precision,intent(in) :: r
        double precision :: tau,arg
        arg = max( epsilon(arg), one-r )
        tau = -log(arg)
end function

elemental function hengreen(r,g)
        double precision,intent(in) :: r,g
        double precision :: hengreen
        double precision,parameter :: one=1.0,two=2.0
        if( approx(g,zero) ) then
          hengreen = two*r-one
        else
          hengreen = one/(two*g) * (one+g**two - ( (one-g**two) / ( two*g*r + one-g) )**two )
        endif
        hengreen = min(max(hengreen,-one), one)
end function

subroutine absorb_photon(p,dist)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dist
        double precision :: new_weight,tau

        tau = get_kabs(p)*dist
        if(tau.gt.20) then
          p%weight = zero
        else
          new_weight = p%weight * exp(- tau)
!          if( (new_weight.gt.one).or.(new_weight.lt.zero) ) then
!            print *,'something wrong with new weight after absorption',new_weight,'=',p%weight,'*',exp(-get_kabs(p)*dist),'(',get_kabs(p),dist,')'
!            call exit
!          endif
          p%weight = new_weight
        endif
end subroutine

subroutine scatter_photon(p)
        type(photon),intent(inout) :: p
        double precision :: muxs,muys,muzs,muxd,muyd,muzd
        double precision :: mutheta,fi,costheta,sintheta,sinfi,cosfi,denom,muzcosfi

        mutheta = hengreen(R(),get_g(p))

        p%scattercnt = p%scattercnt+1

        muxs = p%dir(1)  
        muys = p%dir(2)  
        muzs = p%dir(3)  
        if(p%direct) then
          p%direct=.False.
        endif
        fi      = R()*pi*2.

        costheta = (mutheta)
        sintheta = sqrt(one-costheta**2)
        sinfi = sin(fi)
        cosfi = cos(fi)

        if( muzs .ge. one-1e-8_ireals) then
                muxd = sintheta*cosfi
                muyd = sintheta*sinfi
                muzd = costheta
        else if ( muzs .le. -one+1e-8_ireals) then
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

        if(isnan(muxd).or.isnan(muyd).or.isnan(muzd) ) print *,'new direction is NAN :( --',muxs,muys,muzs,mutheta,fi,'::',muxd,muyd,muzd

        p%dir(1) = muxd
        p%dir(2) = muyd
        p%dir(3) = muzd

end subroutine

pure double precision function get_kabs(p)
        type(photon),intent(in) :: p
        get_kabs = p%optprop(1)
end function
pure double precision function get_ksca(p)
        type(photon),intent(in) :: p
        get_ksca = p%optprop(2)
end function
pure double precision function get_g(p)
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
end subroutine

function R()
    double precision :: R
    R = getRandomReal(rndSeq)
    call random_number(R)
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
  s = int(rn*1000)*myid

!  s=myid
!  print *,myid,'Seeding RNG with ',s
  rndSeq = new_RandomNumberSequence(seed=s)
  lRNGseeded=.True.
end subroutine
  subroutine init_stddev( std, N, rtol, atol)
      type(stddev),intent(inout) :: std
      real(ireals),intent(in) :: atol,rtol
      integer(iintegers) :: N
      if( allocated(std%inc     ) ) deallocate( std%inc  )
      if( allocated(std%delta   ) ) deallocate( std%delta)
      if( allocated(std%mean    ) ) deallocate( std%mean )
      if( allocated(std%mean2   ) ) deallocate( std%mean2)
      if( allocated(std%var     ) ) deallocate( std%var  )
      allocate( std%inc   (N)) ; std%inc   = zero
      allocate( std%delta (N)) ; std%delta = zero
      allocate( std%mean  (N)) ; std%mean  = zero
      allocate( std%mean2 (N)) ; std%mean2 = zero
      allocate( std%var   (N)) ; std%var   = zero
      std%atol = atol
      std%rtol = rtol
  end subroutine

  subroutine std_update(std, N, numnodes)
      type(stddev),intent(inout) :: std
      integer(iintegers) :: N, numnodes
      std%delta = std%inc   - std%mean
      std%mean  = std%mean  + std%delta/N
      std%mean2 = std%mean2 + std%delta * ( std%inc - std%mean )
      std%var = sqrt( std%mean2/max( i1,N ) ) / sqrt( one*N*numnodes )
!      print *,'atol',std%var,'rtol',std%var/max(std%mean,std%rtol)
      if( all( std%var .lt. std%atol .or. std%var/max(std%mean,epsilon(std%rtol)) .lt. std%rtol ) ) then
        std%converged = .True.
      else
        std%converged = .False.
      endif
  end subroutine

  subroutine init(bmc, comm)
        class(t_boxmc) :: bmc
        integer(mpiint),intent(in) :: comm

        if(comm.eq.-1) then
                myid = -1
        else
                call MPI_Comm_rank(comm, myid, ierr)
        endif

        if(.not.lRNGseeded) call init_random_seed(myid+2)

        numnodes=1
        if(myid.ge.0) call mpi_comm_size(comm,numnodes,ierr)

        select type (bmc)
          class is (t_boxmc_8_10)
                  bmc%dir_streams  =  8
                  bmc%diff_streams = 10
          class is (t_boxmc_1_2)
                  bmc%dir_streams  =  1
                  bmc%diff_streams =  2
          class default
              stop 'initialize: unexpected type for boxmc object!'
        end select

        bmc%initialized = .True.
  end subroutine

include 'boxmc_8_10.inc'
include 'boxmc_1_2.inc'

end module
