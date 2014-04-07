
subroutine update_dir_stream_1_2(bmc,p,S,N)
        class(boxmc_1_2) :: bmc
        type(photon),intent(in) :: p
        real(ireals),intent(inout) :: S(:)
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
subroutine update_diff_stream_1_2(bmc,p,S,N)
        class(boxmc_1_2) :: bmc
        type(photon),intent(in) :: p
        real(ireals),intent(inout) :: S(:)
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

subroutine init_dir_photon_1_2(bmc,p,src,direct,initial_dir,dx,dy,dz)
        class(boxmc_1_2) :: bmc
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

subroutine init_diff_photon_1_2(bmc,p,src,dx,dy,dz)
        class(boxmc_1_2) :: bmc
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

subroutine intersect_distance_1_2(bmc,p,max_dist)
        class(boxmc_1_2) :: bmc
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
