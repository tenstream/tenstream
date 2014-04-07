module boxmc
      use helper_functions, only : approx,mean,rmse,deg2rad,norm
      use iso_c_binding
      use mersenne
      use mpi
      use data_parameters, only: iintegers,ireals,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10, zero,one,nil,inil,pi
      
      use boxmc_parameters, only : delta_scale_truncate
      
      implicit none

      private
      public :: qreals,         & ! quad precision reals                                                                           
      fg,bg,tot,                & ! parameter integers for optical properties in cloudfraction boxes                               
      photon, stddev,           & ! struct definitions
      Nphotons_max,             & ! maximum number of photons that are globally started... if that does not reach convergence this gets too expensive...
      mpi_reduce_sum,           & ! mpi_reduce_sum
      roulette,                 & ! roulette kills photons by chance
      delta_scaling,            & ! resets direct or diffuse status depending on incident angle
      s,L,                      & ! small arbitrary helper functions
      update_photon_loc,        & ! shift photon location
      hit_plane,                & ! hit plane with given startlocation and normal vector
      distance,                 & ! dice roll distance, the photon has to travel
      tau,hengreen,             & ! random optical depth and henyey greenstein phase function
      absorb_photon,scatter_photon, & ! multiply photon weight with absorption and calculate scattering angle
      get_kabs,get_ksca,get_g,  & ! get'er functions for optical props
      print_photon,             & ! debug print statement
      init_random_seed,R,       & ! init mersenne RNG
      lRNGseeded,               & ! logical if init RNG
      init_stddev,std_update      ! standard deviation routines
      
      real(c_float128) :: quad_precision
      integer,parameter :: qreals=ireals

      integer,parameter :: fg=1,bg=2,tot=3

      type photon
              double precision :: loc(3)=nil,dir(3)=nil,weight=nil,dx=nil,dy=nil,dz=nil
              logical :: alive=.True.,direct=.False.
              integer(iintegers) :: side=inil,src=inil,scattercnt=0
              double precision :: optprop(3) ! kabs,ksca,g
      end type

      integer(iintegers),parameter :: Nphotons_max=int(1e9)

      type stddev
        real(qreals),allocatable,dimension(:) :: inc,delta,mean,mean2,var
        logical :: converged=.False.
        real(qreals) :: atol=zero, rtol=zero
      end type

      type(stddev) :: std_Sdir, std_Sdiff, std_abso

      integer :: ierr,myid,numnodes

      logical :: lRNGseeded=.False.
      type(randomNumberSequence) :: rndSeq


contains

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
subroutine delta_scaling(p,deltascale,initial_dir)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: initial_dir(3)
        logical,intent(in) :: deltascale

        double precision :: angle
        if((.not.deltascale).or.p%direct) return

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
      real(qreals),intent(in) :: atol,rtol
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


end module
