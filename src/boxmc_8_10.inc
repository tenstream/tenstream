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


subroutine update_dir_stream_8_10(bmc,p,S,N)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(in) :: p
        real(ireals),intent(inout) :: S(:)
        integer(iintegers),intent(inout) :: N(:)
        real(ireals),parameter :: h=.5_ireals

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
      if(p%loc(1).le.p%dx*h.and.p%loc(2).le.p%dy*h) then
        S(1) = S(1)+p%weight
        N(1) = N(1)+i1
        return
      else if(p%loc(1).gt.p%dx*h.and.p%loc(2).le.p%dy*h) then
        S(2) = S(2)+p%weight
        N(2) = N(2)+i1
        return
      else if(p%loc(1).le.p%dx*h.and.p%loc(2).gt.p%dy*h) then
        S(3) = S(3)+p%weight
        N(3) = N(3)+i1
        return
      else if(p%loc(1).gt.p%dx*h.and.p%loc(2).gt.p%dy*h) then
        S(4) = S(4)+p%weight
        N(4) = N(4)+i1
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif

    case(3:4)

      if(p%loc(3).le.p%dz*h ) then
        S(5) = S(5)+p%weight
        N(5) = N(5)+i1
        return
      else if(p%loc(3).gt.p%dz*h ) then
        S(6) = S(6)+p%weight
        N(6) = N(6)+i1
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif
    case(5:6)

      if(p%loc(3).le.p%dz*h ) then
        S(7) = S(7)+p%weight
        N(7) = N(7)+i1
        return
      else if(p%loc(3).gt.p%dz*h ) then
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
        select type(bmc)
        end select

end subroutine
subroutine update_diff_stream_8_10(bmc,p,S,N)
        class(t_boxmc_8_10) :: bmc
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
select type(bmc)
end select

end subroutine

subroutine init_dir_photon_8_10(bmc,p,src,direct,initial_dir,dx,dy,dz)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(inout) :: p
        real(ireals),intent(in) :: dx,dy,dz,initial_dir(3)
        integer(iintegers),intent(in) :: src
        logical,intent(in) :: direct
        real(ireals),parameter :: h=.5_ireals
        real(ireals),parameter :: e=1e-8_ireals

        p%alive = .False.

        select case(src)
        case (1)
                p%loc = (/L(dx)*h     , L(dy)*h     ,    dz-e  /)
        case (2)
                p%loc = (/dx*h+L(dx)*h, L(dy)*h     ,    dz-e  /)
        case (3)
                p%loc = (/L(dx)*h     , dy*h+L(dy)*h,    dz-e  /)
        case (4)
                p%loc = (/dx*h+L(dx)*h, dy*h+L(dy)*h,    dz-e  /)
        case (5)
                p%loc = (/  zero+e      , L(dy)       , L(dz)*h/)
        case (6)
                p%loc = (/  zero+e      , L(dy)       , dz*h+L(dz)*h/)
        case (7)
                p%loc = (/L(dx)       ,   zero+e      , L(dz)*h/)
        case (8)
                p%loc = (/L(dx)       ,   zero+e      , dz*h+L(dz)*h/)
        case default
                print *,'Dont know what to do with source spec:',src
                call exit
        end select

        p%weight=one
        p%dx   = dx
        p%dy   = dy
        p%dz   = dz
        p%alive = .True.
        p%direct= direct
        p%side = int(nil)
        p%src  = src
        p%dir = initial_dir
        select type(bmc)
        end select

end subroutine

subroutine init_diff_photon_8_10(bmc,p,src,dx,dy,dz)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(inout) :: p
        real(ireals),intent(in) :: dx,dy,dz
        integer(iintegers),intent(in) :: src

        real(ireals) :: phi,mu,mu1,mu2

        mu1= sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
        mu2 = R() ! evenly distribution between 0 and 90 degree

        mu=mu1
        select case(src)
        case (-1) ! thermal emission
                p%loc = (/L(dx), L(dy), L(dz) /)
                mu = 2._ireals*mu2 - one
                phi= deg2rad( R()*360._ireals )
        case (1)
                p%loc = (/L(dx), L(dy),    dz  /)
                phi= deg2rad( R()*360._ireals )
        case (2)
                p%loc = (/L(dx), L(dy),    zero /)
                phi= deg2rad( R()*360._ireals )
        case (3)
                p%loc = (/  zero   , L(dy), L(dz)/)
                phi= deg2rad( interv_R(-90._ireals,90._ireals) )
        case (4)
                p%loc = (/ dx   , L(dy), L(dz)/)
                phi= deg2rad( interv_R(90._ireals,270._ireals) )
        case (5)
                p%loc = (/  zero   , L(dy), L(dz)/)
                phi= deg2rad( interv_R(90._ireals,270._ireals) )
        case (6)
                p%loc = (/ dx   , L(dy), L(dz)/)
                phi= deg2rad( interv_R(-90._ireals,90._ireals) )
        case (7)
                p%loc = (/L(dx),   zero   , L(dz)/)
                phi= deg2rad( interv_R(0._ireals,180._ireals) )
        case (8)
                p%loc = (/L(dx),  dy   , L(dz)/)
                phi= deg2rad( interv_R(180._ireals,360._ireals) )
        case (9)
                p%loc = (/L(dx),   zero   , L(dz)/)
                phi= deg2rad( interv_R(180._ireals,360._ireals) )
        case (10)
                p%loc = (/L(dx),  dy   , L(dz)/)
                phi= deg2rad( interv_R(0._ireals,180._ireals) )
        case default
                print *,'Dont know what to do with source spec:',src
                call exit
        end select
        select case(src)
          case(1:2)
            phi= deg2rad( R()*360._ireals )
          case(3:10)
            phi= deg2rad( interv_R(-90._ireals,90._ireals) )
        end select

        p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)

        select case(src)
          case(1)
            p%dir = rotate_y(p%dir,180)
          case(3)
            p%dir = rotate_y(p%dir,-90)
          case(4)
            p%dir = rotate_y(p%dir,-90)
            p%dir = rotate_z(p%dir,180)
          case(5)
            p%dir = rotate_y(p%dir,90)
            p%dir = rotate_z(p%dir,180)
          case(6)
            p%dir = rotate_y(p%dir,90)

          case(7)
            p%dir = rotate_y(p%dir,-90)
            p%dir = rotate_z(p%dir,-90)
          case(8)
            p%dir = rotate_y(p%dir,-90)
            p%dir = rotate_z(p%dir,180)
            p%dir = rotate_z(p%dir,-90)
          case(9)
            p%dir = rotate_y(p%dir,90)
            p%dir = rotate_z(p%dir,90)
          case(10)
            p%dir = rotate_y(p%dir,90)
            p%dir = rotate_z(p%dir,-90)

        end select
        p%weight=one
        p%dx   = dx
        p%dy   = dy
        p%dz   = dz
        p%alive = .True.
        p%direct= .False.
        p%side = int(nil)
        p%src  = src
        if(.not. approx( norm(p%dir) ,one ) ) then
          print *,'p%dir',p%dir,'norm',norm(p%dir)
          stop 'init_diff_photon_8_10 :: direction was not normed?!'
        endif
        select type(bmc)
        end select

        contains
          function rotate_x(v,angle)
              real(ireals) :: rotate_x(3)
              real(ireals),intent(in) :: v(3)
              integer :: angle
              real(ireals) :: M(3,3),s,c
              s=sin(deg2rad(one*angle))
              c=cos(deg2rad(one*angle))

              M(1,:)=[one ,zero ,zero]
              M(2,:)=[zero, c   , s  ]
              M(3,:)=[zero,-s   , c  ]

              rotate_x = matmul(M,v)
          end function
          function rotate_y(v,angle)
              real(ireals) :: rotate_y(3)
              real(ireals),intent(in) :: v(3)
              integer :: angle
              real(ireals) :: M(3,3),s,c
              s=sin(deg2rad(one*angle))
              c=cos(deg2rad(one*angle))

              M(1,:)=[ c  ,zero , -s ]
              M(2,:)=[zero, one ,zero]
              M(3,:)=[ s  , zero, c  ]

              rotate_y = matmul(M,v)
          end function
          function rotate_z(v,angle)
              real(ireals) :: rotate_z(3)
              real(ireals),intent(in) :: v(3)
              integer :: angle
              real(ireals) :: M(3,3),s,c
              s=sin(deg2rad(one*angle))
              c=cos(deg2rad(one*angle))

              M(1,:)=[ c  , s   ,zero]
              M(2,:)=[-s  , c   ,zero]
              M(3,:)=[zero, zero, one]

              rotate_z = matmul(M,v)
          end function

end subroutine

subroutine intersect_distance_8_10(bmc,p,max_dist)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(inout) :: p
        real(ireals),intent(out) :: max_dist

        real(ireal_dp) :: x,y,z
        integer(iintegers) :: i,sides(3)

        real(ireal_dp) :: dist(3)

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
          max_dist=huge(max_dist)
          do i=1,3
            if(.not. approx(p%dir(i),zero) ) then
              if( dist(i).le.max_dist ) then
                p%side = sides(i)
                max_dist = dist(i)
              endif
            endif
          enddo

          !print *,'should actually not be here at the end of crossings in intersect distance! - however, please check if distance makes sense?:',max_dist,'::',dist
          !call print_photon(p)

        select type(bmc)
        end select

end subroutine 
  ! ******************* BOXMC 8 10 ************************

