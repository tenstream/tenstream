module boxmc_8_10
      use helper_functions, only : approx,mean,rmse,deg2rad,norm
      use iso_c_binding
      use mersenne
      use mpi
      use data_parameters, only: iintegers,ireals,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10, zero,one,nil,inil,pi
      use boxmc_parameters_8_10,only : dir_streams,diff_streams
      use boxmc_parameters, only : delta_scale_truncate

      use boxmc, only :  qreals,& ! quad precision reals                                                                           
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


      implicit none

      private
      public :: bmc_get_coeff_8_10

      type(stddev) :: std_Sdir, std_Sdiff, std_abso

      integer :: ierr,myid,numnodes
        
contains
subroutine bmc_get_coeff_8_10(comm,op_bg,src,S_out,Sdir_out,dir,deltascale,phi0,theta0,dx,dy,dz)
        double precision,intent(in) :: op_bg(3),phi0,theta0
        logical,intent(in) :: deltascale
        integer(iintegers),intent(in) :: src
        integer :: comm
        logical,intent(in) :: dir
        double precision,intent(in) :: dx,dy,dz
        double precision,intent(out):: S_out(10),Sdir_out(8)

        type(photon) :: p
        integer(iintegers) :: k,mycnt
        double precision :: time(2),total_photons,initial_dir(3)
        integer(iintegers) :: Ndir(8),Ndiff(10)

        Ndir=i0;Ndiff=i0
        
        call init_stddev( std_Sdir , dir_streams  ,1e-5_qreals, 5e-3_qreals )
        call init_stddev( std_Sdiff, diff_streams ,1e-5_qreals, 5e-3_qreals )
        call init_stddev( std_abso , i1           ,1e-5_qreals, 5e-4_qreals )

        if(.not.dir) std_Sdir%converged=.True.

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

        if(.not.lRNGseeded) call init_random_seed(myid+2)

        if(dx.le.zero .or. dy.le.zero .or. dz.le.zero ) print *,'ERROR: box dimensions have to be positive!',dx,dy,dz

        numnodes=1
        if(myid.ge.0) call mpi_comm_size(comm,numnodes,ierr)
!        print *,myid,'starting photon tracing',numnodes,'doing',cnt/numnodes,'in total'

        call cpu_time(time(1))

        mycnt = max(i1*numnodes,Nphotons_max/numnodes)
        do k=1,mycnt
                if(k*numnodes.gt.1e4 .and. all([std_Sdir%converged, std_Sdiff%converged, std_abso%converged ]) ) exit 

                if(dir) then
                        call init_dir_photon_8_10(p,src,dir,initial_dir,dx,dy,dz)
                else
                        call init_photon_8_10(p,src,dx,dy,dz)
                endif
                p%optprop = op_bg

                move: do 
                        call move_photon_8_10(p)
                        call roulette(p)
                
                        if(.not.p%alive) exit move
                        call scatter_photon(p)
                enddo move

                if(dir) call delta_scaling(p,deltascale,initial_dir)

                std_abso%inc = one-p%weight
                std_Sdir%inc  = zero
                std_Sdiff%inc = zero

                if(p%direct) then
                  call update_direct_stream_8_10(p,std_Sdir%inc,Ndir)
                else
                  call update_stream_8_10(p,std_Sdiff%inc,Ndiff)
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

        if(myid.le.0.and.total_photons.ge.1e7) print *,src,dz,op_bg,'angles',phi0,theta0,'took',time(2)-time(1),'s',' phots*1e3 :: ',total_photons/1e3,' abso :',one-sum(Sdir_out)-sum(S_out),':',total_photons/(time(2)-time(1))/numnodes,'phots/sec/node'
!        if(myid.le.0) write(*, FMT='(" src ", I0, " optprop ", 3(f12.8), " angles ",2(f5.1), " total ", I0, " photons in ", (f5.2), "s (", I0, " p/s)" )' ) &
!            src,op_bg,[phi0,theta0],int(total_photons),time(2)-time(1),int(total_photons/(time(2)-time(1)))/numnodes
end subroutine

subroutine update_direct_stream_8_10(p,S,N)
        type(photon),intent(in) :: p
        real(qreals),intent(inout) :: S(:)
        integer(iintegers),intent(inout) :: N(:)

!            1       2
!               3       4
!         _______________
!        |               |
!      6 |           8   | 
!        |               |
!        |               |
!        |    7          |
!      5 |_______________| 
    select case (p%side)
    case(1:2)
      if(p%loc(1).le.p%dx/2.and.p%loc(2).le.p%dy/2) then
        S(1) = S(1)+p%weight
        N(1) = N(1)+i1
        return
      else if(p%loc(1).gt.p%dx/2.and.p%loc(2).le.p%dy/2) then
        S(2) = S(2)+p%weight
        N(2) = N(2)+i1
        return
      else if(p%loc(1).le.p%dx/2.and.p%loc(2).gt.p%dy/2) then
        S(3) = S(3)+p%weight
        N(3) = N(3)+i1
        return
      else if(p%loc(1).gt.p%dx/2.and.p%loc(2).gt.p%dy/2) then
        S(4) = S(4)+p%weight
        N(4) = N(4)+i1
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif

    case(3:4)

      if(p%loc(3).le.p%dz/2 ) then
        S(5) = S(5)+p%weight
        N(5) = N(5)+i1
        return
      else if(p%loc(3).gt.p%dz/2 ) then
        S(6) = S(6)+p%weight
        N(6) = N(6)+i1
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif
    case(5:6)

      if(p%loc(3).le.p%dz/2 ) then
        S(7) = S(7)+p%weight
        N(7) = N(7)+i1
        return
      else if(p%loc(3).gt.p%dz/2 ) then
        S(8) = S(8)+p%weight
        N(8) = N(8)+i1
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif
    case default
      print *,'Dont know what to do with this p%side'
      call print_photon(p)
    end select
end subroutine
subroutine update_stream_8_10(p,S,N)
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

        else if(p%side.eq.3 .and. p%dir(3).le.zero ) then
                S(3) = S(3)+p%weight
                N(3) = N(3)+i1
                return
        else if(p%side.eq.3 .and. p%dir(3).gt.zero ) then
                S(5) = S(5)+p%weight
                N(5) = N(5)+i1
                return

        else if(p%side.eq.4 .and. p%dir(3).le.zero ) then
                S(4) = S(4)+p%weight
                N(4) = N(4)+i1
                return
        else if(p%side.eq.4 .and. p%dir(3).gt.zero ) then
                S(6) = S(6)+p%weight
                N(6) = N(6)+i1
                return

        else if(p%side.eq.5 .and. p%dir(3).le.zero ) then
                S(7) = S(7)+p%weight
                N(7) = N(7)+i1
                return
        else if(p%side.eq.6 .and. p%dir(3).le.zero ) then
                S(8) = S(8)+p%weight
                N(8) = N(8)+i1
                return
        else if(p%side.eq.5 .and. p%dir(3).gt.zero ) then
                S(9) = S(9)+p%weight
                N(9) = N(9)+i1
                return
        else if(p%side.eq.6 .and. p%dir(3).gt.zero ) then
                S(10) =S(10)+p%weight
                N(10) =N(10)+i1
                return
        else
                print *,'Couldnt find a stream on to which I can put the photon weight on?!'
                call print_photon(p)
                call exit
        endif
end subroutine

subroutine init_dir_photon_8_10(p,src,direct,initial_dir,dx,dy,dz)
        type(photon),intent(inout) :: p
        double precision,intent(in) :: dx,dy,dz,initial_dir(3)
        integer(iintegers),intent(in) :: src
        logical,intent(in) :: direct

        p%alive = .False.

        if(src.eq.1) then
                p%loc = (/L(dx)/2     , L(dy)/2     ,    dz  /)
        else if(src.eq.2) then
                p%loc = (/dx/2+L(dx)/2, L(dy)/2     ,    dz  /)
        else if(src.eq.3) then
                p%loc = (/L(dx)/2     , dy/2+L(dy)/2,    dz  /)
        else if(src.eq.4) then
                p%loc = (/dx/2+L(dx)/2, dy/2+L(dy)/2,    dz  /)
        else if(src.eq.5) then
                p%loc = (/  zero      , L(dy)       , L(dz)/2/)
        else if(src.eq.6) then
                p%loc = (/  zero      , L(dy)       , dz/2+L(dz)/2/)
        else if(src.eq.7) then
                p%loc = (/L(dx)       ,   zero      , L(dz)/2/)
        else if(src.eq.8) then
                p%loc = (/L(dx)       ,   zero      , dz/2+L(dz)/2/)
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

subroutine init_photon_8_10(p,src,dx,dy,dz)
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
        else if(src.eq.3) then
                p%loc = (/  O   , L(dy), L(dz)/)
                p%dir = (/  I   , s()   , -R()   /)
        else if(src.eq.4) then
                p%loc = (/ dx   , L(dy), L(dz)/)
                p%dir = (/ -I   , s()   , -R()   /)
        else if(src.eq.5) then
                p%loc = (/  O   , L(dy), L(dz)/)
                p%dir = (/  I   , s()   ,R()   /)
        else if(src.eq.6) then
                p%loc = (/ dx   , L(dy), L(dz)/)
                p%dir = (/ -I   , s()   ,R()   /)
        else if(src.eq.7) then
                p%loc = (/L(dx),   O   , L(dz)/)
                p%dir = (/s()  ,   I   , -R()   /)
        else if(src.eq.8) then
                p%loc = (/L(dx),  dy   , L(dz)/)
                p%dir = (/s()  , -I    , -R()   /)
        else if(src.eq.9) then
                p%loc = (/L(dx),   O   , L(dz)/)
                p%dir = (/s()  ,   I   ,R()   /)
        else if(src.eq.10) then
                p%loc = (/L(dx),  dy   , L(dz)/)
                p%dir = (/s()  , -I    ,R()   /)
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

subroutine move_photon_8_10(p)
        type(photon),intent(inout) :: p
        double precision :: tau_travel,dist,intersec_dist

        tau_travel = tau(R())
        call intersect_distance_8_10(p,intersec_dist)

        dist = distance(tau(R() ), get_ksca(p) )

        if(intersec_dist .le. dist) then
          p%alive=.False.
          call update_photon_loc(p,intersec_dist)
        else
          call update_photon_loc(p,dist)
          return
        endif
end subroutine

subroutine intersect_distance_8_10(p,max_dist)
        type(photon),intent(inout) :: p
        double precision,intent(out) :: max_dist

        double precision :: x,y,z
        integer(iintegers) :: i,sides(3)

        double precision :: dist(3)

        !crossing with bottom and top plane:
          if(p%dir(3).ge.zero) then
            max_dist = hit_plane(p,[zero,zero,p%dz ],[zero,zero,one])
            p%side=1 
            x = p%loc(1)+p%dir(1)*max_dist
            y = p%loc(2)+p%dir(2)*max_dist
            if( ( x.ge.zero .and. x.le.p%dx) .and. ( y.ge.zero .and. y.le.p%dy) ) return
            dist(1) = max_dist; sides(1) = 1
          endif
          if(p%dir(3).le.zero) then
            max_dist = hit_plane(p,[zero,zero,zero ],[zero,zero,one])
            p%side=2
            x = p%loc(1)+p%dir(1)*max_dist
            y = p%loc(2)+p%dir(2)*max_dist
            if( ( x.ge.zero .and. x.le.p%dx) .and. ( y.ge.zero .and. y.le.p%dy) ) return
            dist(1) = max_dist; sides(1) = 2
          endif

          !crossing with left and right plane:
          if(p%dir(1).le.zero) then
            max_dist = hit_plane(p,[ zero ,zero,zero],[one,zero,zero])
            p%side=3
            y = p%loc(2)+p%dir(2)*max_dist
            z = p%loc(3)+p%dir(3)*max_dist
            if( ( y.ge.zero .and. y.le.p%dy) .and. ( z.ge.zero .and. z.le.p%dz) ) return
            dist(2) = max_dist; sides(2) = 3
          endif
          if(p%dir(1).ge.zero) then
            max_dist = hit_plane(p,[ p%dx ,zero,zero],[one,zero,zero])
            p%side=4
            y = p%loc(2)+p%dir(2)*max_dist
            z = p%loc(3)+p%dir(3)*max_dist
            if( ( y.ge.zero .and. y.le.p%dy) .and. ( z.ge.zero .and. z.le.p%dz) ) return
            dist(2) = max_dist; sides(2) = 4
          endif

          !crossing with back and forward plane:
          if(p%dir(2).le.zero) then
            max_dist = hit_plane(p,[zero, zero ,zero],[zero,one,zero])
            p%side=5
            x = p%loc(1)+p%dir(1)*max_dist
            z = p%loc(3)+p%dir(3)*max_dist
            if( ( x.ge.zero .and. x.le.p%dx) .and. ( z.ge.zero .and. z.le.p%dz) ) return
            dist(3) = max_dist; sides(3) = 5
          endif
          if(p%dir(2).ge.zero) then
            max_dist = hit_plane(p,[zero, p%dy ,zero],[zero,one,zero])
            p%side=6
            x = p%loc(1)+p%dir(1)*max_dist
            z = p%loc(3)+p%dir(3)*max_dist
            if( ( x.ge.zero .and. x.le.p%dx) .and. ( z.ge.zero .and. z.le.p%dz) ) return
            dist(3) = max_dist; sides(3) = 6
          endif

          !Ohhh there was a problem.. maybe with numerics, seems that it may happen that we dont find a solution if norm of p%dir is not equal to one....
          max_dist=huge(dist)
          do i=1,3
            if(p%dir(i).gt.1e-12_ireals.and.p%dir(i).lt.-1e-12_ireals) then
              if( dist(i).le.max_dist ) then
                p%side = sides(i)
                max_dist = dist(i)
              endif
            endif
          enddo

          print *,'should actually not be here at the end of crossings in intersect distance! - however, please check if distance makes sense?:',max_dist
          call print_photon(p)

end subroutine

end module
