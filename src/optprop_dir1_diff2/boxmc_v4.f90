module boxmc_1_2
      use helper_functions, only : approx,mean,rmse,deg2rad
      use iso_c_binding
      use mersenne
      use mpi
      use data_parameters, only: iintegers,ireals,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10, nil,inil, one,zero,pi
      use boxmc_parameters_1_2, only: dir_streams,diff_streams, delta_scale_truncate
      implicit none

      private
      public :: bmc_get_coeff

      real(c_float128) :: quad_precision
!      integer,parameter :: qreals=kind(quad_precision)
      integer,parameter :: qreals=ireals

      integer,parameter :: fg=1,bg=2,tot=3,dp=kind(one)

      type photon
              double precision :: loc(3)=nil,dir(3)=nil,weight=nil,dx=nil,dy=nil,dz=nil
              logical :: alive=.True.,direct=.False.
              integer(iintegers) :: side=inil,src=inil,scattercnt=0
              double precision :: optprop(3) ! kabs,ksca,g
      end type

      integer(iintegers),parameter :: Nphotons=int(1e9)

      logical :: RNGseeded=.False.
      type(randomNumberSequence) :: rndSeq

      type stddev
        real(qreals),allocatable,dimension(:) :: inc,delta,mean,mean2,var
        logical :: converged=.False.
        real(qreals) :: atol=zero, rtol=zero
      end type

      type(stddev) :: std_Sdir, std_Sdiff, std_abso

      integer :: ierr,myid,numnodes
        
contains
subroutine bmc_get_coeff(comm,op_bg,src,S_out,Sdir_out,dir,deltascale,phi0,theta0,dx,dy,dz)
        double precision,intent(in) :: op_bg(3),phi0,theta0
        logical,intent(in) :: deltascale
        integer,intent(in) :: comm
        integer(iintegers),intent(in) :: src
        logical,intent(in) :: dir
        double precision,intent(in) :: dx,dy,dz
        double precision,intent(out):: S_out(diff_streams),Sdir_out(dir_streams)

        type(photon) :: p
        integer(iintegers) :: k,mycnt
        double precision :: time(2),total_photons,initial_dir(3)
        integer(iintegers) :: Ndir(dir_streams),Ndiff(diff_streams)

        Ndir=i0;Ndiff=i0
        
        call init_stddev( std_Sdir , dir_streams  ,1e-5_qreals, 5e-3_qreals )
        call init_stddev( std_Sdiff, diff_streams ,1e-5_qreals, 5e-3_qreals )
        call init_stddev( std_abso , i1           ,1e-5_qreals, 5e-4_qreals )

        initial_dir  = (/ sin(deg2rad(theta0))*sin(deg2rad(phi0)) ,&
                    sin(deg2rad(theta0))*cos(deg2rad(phi0)) ,&
                  - cos(deg2rad(theta0)) /)
        initial_dir = initial_dir/norm(initial_dir)

        if( (any(op_bg.lt.zero)) .or. (any(isnan(op_bg))) ) then
          print *,'corrupt optical properties: bg:: ',op_bg
          call exit
        endif

        if(comm.eq.-1) then
                myid = -1
        else
                call MPI_Comm_rank(comm, myid, ierr)
        endif

        if(.not.RNGseeded) call init_random_seed(myid+2)

        if(dx.le.zero .or. dy.le.zero .or. dz.le.zero ) print *,'ERROR: box dimensions have to be positive!',dx,dy,dz

        numnodes=1
        if(myid.ge.0) call mpi_comm_size(comm,numnodes,ierr)
!        print *,myid,'starting photon tracing',numnodes,'doing',cnt/numnodes,'in total'

        call cpu_time(time(1))

        mycnt = max(i1*numnodes,Nphotons/numnodes)
        do k=1,mycnt
                if(k*numnodes.gt.1e4 .and. all([std_Sdir%converged, std_Sdiff%converged, std_abso%converged ]) ) exit 

                if(dir) then
                        call init_dir_photon(p,src,dir,initial_dir,dx,dy,dz)
                else
                        call init_photon(p,src,dx,dy,dz)
                endif
                p%optprop = op_bg

                move: do 
                        call move_photon(p)
                        call roulette(p)
                
                        if(.not.p%alive) exit move
                        call scatter_photon(p)
                enddo move

                if(dir) call delta_scaling(p,deltascale,initial_dir)

                std_abso%inc = one-p%weight
                std_Sdir%inc  = zero
                std_Sdiff%inc = zero

                if(p%direct) then
                  call update_direct_stream(p,std_Sdir%inc,Ndir)
                else
                  call update_stream(p,std_Sdiff%inc,Ndiff)
                endif

                call std_update( std_abso , k, i1*numnodes)
                call std_update( std_Sdir , k, i1*numnodes )
                call std_update( std_Sdiff, k, i1*numnodes )
        enddo
        
        total_photons=k
        S_out    = dble( std_Sdiff%mean*k )
        Sdir_out = dble( std_Sdir%mean *k )
        if(myid.ge.0) then
                call mpi_reduce_sum(total_photons,comm,myid)
                do k=1,ubound(S_out,1)
                        call mpi_reduce_sum(S_out(k),comm,myid)
                enddo
                do k=1,ubound(Sdir_out,1)
                        call mpi_reduce_sum(Sdir_out(k),comm,myid)
                enddo
        endif
        S_out    = S_out / total_photons
        Sdir_out = Sdir_out / total_photons

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

        if(myid.le.0.and.total_photons.ge.1e7) then
          print *,src,dz,op_bg,'angles',phi0,theta0,'took',time(2)-time(1),'s',' phots*1e3 :: ',total_photons/1e3,' abso :',one-sum(Sdir_out)-sum(S_out),':',total_photons/(time(2)-time(1))/numnodes,'phots/sec/node'
        endif
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
subroutine update_direct_stream(p,S,N)
        type(photon),intent(in) :: p
        real(qreals),intent(inout) :: S(:)
        integer(iintegers),intent(inout) :: N(:)

    select case (p%side)
    case(2)
        S(1) = S(1)+p%weight
        N(1) = N(1)+i1
        return
    case default
      print *,'Dont know what to do with this p%side'
      call print_photon(p)
    end select
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
subroutine update_stream(p,S,N)
        type(photon),intent(in) :: p
        real(qreals),intent(inout) :: S(:)
        integer(iintegers),intent(inout) :: N(:)

!         _______1_______
!        |           10  |
!      5 |            8  | 6
!        |               |
!        |   9           |
!        |   7           |
!      3 |_______________| 4
!                2

        if(p%side.eq.1) then
                S(1) = S(1)+p%weight
                N(1) = N(1)+i1
                return
        else if(p%side.eq.2) then
                S(2) = S(2)+p%weight
                N(2) = N(2)+i1
                return
        else
                print *,'Couldnt find a stream on to which I can put the photon weight on?!'
                call print_photon(p)
                call exit
        endif
end subroutine

double precision function s()
        s = -one + 2*R()
end function
double precision function L(v)
        double precision,intent(in) ::v
        double precision,parameter :: eps=1e-3
        L = min( max(R()*v,eps), v-eps)
end function

subroutine init_dir_photon(p,src,direct,initial_dir,dx,dy,dz)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dx,dy,dz,initial_dir(3)
        integer(iintegers),intent(in) :: src
        logical,intent(in) :: direct

        p%alive = .False.

        if(src.eq.1) then
                p%loc = (/L(dx)     , L(dy)     ,    dz  /)
        else
                print *,'Dont know what to do with source spec:',src
                call exit
        endif

        p%weight=one
        p%dx   = dx
        p%dy   = dy
        p%dz   = dz
        p%alive = .True.
        p%direct= direct
        p%side = int(nil)
        p%src  = src
        p%dir = initial_dir
end subroutine

subroutine init_photon(p,src,dx,dy,dz)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dx,dy,dz
        integer(iintegers),intent(in) :: src
        double precision,parameter :: I=1.,O=0.

        if(src.eq.1) then
                p%loc = (/L(dx), L(dy),    dz  /)
                p%dir = (/s()   , s()   ,   -I  /)
        else if(src.eq.2) then
                p%loc = (/L(dx), L(dy),    O /)
                p%dir = (/s()   , s()   ,    I /)
        else
                print *,'Dont know what to do with source spec:',src
                call exit
        endif

        p%weight=one
        p%dx   = dx
        p%dy   = dy
        p%dz   = dz
        p%alive = .True.
        p%direct= .False.
        p%side = int(nil)
        p%src  = src
        p%dir = p%dir/norm(p%dir)

end subroutine

pure double precision function norm(v)
        double precision,intent(in) :: v(:)
        norm = sqrt(dot_product(v,v))
end function

subroutine move_photon(p)
        type(photon),intent(inout) :: p

        double precision :: tau_travel,dist,intersec_dist

        tau_travel = tau(R())
        call intersect_distance(p,intersec_dist)

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

subroutine intersect_distance(p,max_dist)
        type(photon),intent(inout) :: p
        double precision,intent(out) :: max_dist

        double precision :: x,y
!        double precision,parameter :: eps=1e-3

        !crossing with bottom and top plane:
          if(p%dir(3).ge.zero) then
            max_dist = hit_plane(p,[zero,zero,p%dz ],[zero,zero,one])
            p%side=1 
            x = p%loc(1)+p%dir(1)*max_dist
            y = p%loc(2)+p%dir(2)*max_dist
            return
          endif
          if(p%dir(3).le.zero) then
            max_dist = hit_plane(p,[zero,zero,zero ],[zero,zero,one])
            p%side=2
            x = p%loc(1)+p%dir(1)*max_dist
            y = p%loc(2)+p%dir(2)*max_dist
            return
          endif

          print *,'should actually not be here at the end of crossings in intersect distance! - however, please check if distance makes sense?:',max_dist
          call print_photon(p)

end subroutine

pure double precision function hit_plane(p,po,pn)
        type(photon),intent(in) :: p
        double precision,intent(in) :: po(3),pn(3)
        double precision :: discr
        discr = dot_product(p%dir,pn)
        if( approx(discr,zero) ) then
                hit_plane=huge(hit_plane)
        else        
                hit_plane = dot_product(po-p%loc, pn) / discr
        endif
end function

elemental function distance(tau,beta)
        double precision,intent(in) :: tau,beta
        double precision :: distance
        if( approx( beta,zero ) ) then
          distance=huge(distance)
        else
          distance = tau/beta
        endif

end function

elemental function tau(r)
        double precision,intent(in) :: r
        double precision :: tau
        tau = -log( max( epsilon(tau), one-r ) ) 
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

function R() 
      double precision :: R
!      call random_number(R)
      R = getRandomReal(rndSeq)
end function

subroutine absorb_photon(p,dist)!,tau_abs)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dist
        double precision :: new_weight,tau

        tau = get_kabs(p)*dist
        if(tau.gt.20) then
          p%weight = zero
        else
          new_weight = p%weight * exp(-tau)
!          if( (new_weight.gt.one).or.(new_weight.lt.zero) ) then
!            print *,'something wrong with new weight after absorption',new_weight,'=',p%weight,'*',exp(-tau),'(',get_kabs(p),dist,')'
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

        if( muzs .ge. one-epsilon(muzs)) then
                muxd = sintheta*cosfi
                muyd = sintheta*sinfi
                muzd = costheta
        else if ( muzs .le. -one+epsilon(muzs) ) then
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

subroutine init_random_seed(myid)
  integer myid
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
  rndseq = new_RandomNumberSequence(seed=s)
  RNGseeded=.True.
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
      if( all( std%var .lt. std%atol .or. std%var/std%mean .lt. std%rtol ) ) then
        std%converged = .True.
      else
        std%converged = .False.
      endif
  end subroutine


end module
