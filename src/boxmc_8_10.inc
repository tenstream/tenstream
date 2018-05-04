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


subroutine update_dir_stream_8_10(bmc,p,T)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(in) :: p
        real(ireal_dp),intent(inout) :: T(:)
        real(ireal_dp),parameter :: h=.5_ireal_dp

!          /     3       4 /
!         / 1       2     /
!        /_______________/
!        |               |
!      6 |           8   |
!        |               |
!        |               |
!        |    7          |
!      5 |_______________|
    select case (p%side)
    case(1:2)
      if(p%loc(1).le.p%dx*h.and.p%loc(2).le.p%dy*h) then
        T(1) = T(1)+p%weight
        return
      else if(p%loc(1).gt.p%dx*h.and.p%loc(2).le.p%dy*h) then
        T(2) = T(2)+p%weight
        return
      else if(p%loc(1).le.p%dx*h.and.p%loc(2).gt.p%dy*h) then
        T(3) = T(3)+p%weight
        return
      else if(p%loc(1).gt.p%dx*h.and.p%loc(2).gt.p%dy*h) then
        T(4) = T(4)+p%weight
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif

    case(3:4)

      if(p%loc(3).le.p%dz*h ) then
        T(5) = T(5)+p%weight
        return
      else if(p%loc(3).gt.p%dz*h ) then
        T(6) = T(6)+p%weight
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif
    case(5:6)

      if(p%loc(3).le.p%dz*h ) then
        T(7) = T(7)+p%weight
        return
      else if(p%loc(3).gt.p%dz*h ) then
        T(8) = T(8)+p%weight
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
subroutine update_diff_stream_8_10(bmc,p,S)
        class(t_boxmc_8_10) :: bmc
        type(photon),intent(in) :: p
        real(ireal_dp),intent(inout) :: S(:)

!         _______2_______
!        |            9  |
!      6 |            7  | 5
!        |               |
!        |   10          |
!        |   8           |
!      4 |_______________| 3
!                1

if(p%side.eq.1) then
  S(1) = S(1)+p%weight
  return

else if(p%side.eq.2) then
  S(2) = S(2)+p%weight
  return

else if(p%side.eq.3 .and. p%dir(3).le.zero ) then
  S(3) = S(3)+p%weight
  return
else if(p%side.eq.3 .and. p%dir(3).gt.zero ) then
  S(5) = S(5)+p%weight
  return

else if(p%side.eq.4 .and. p%dir(3).le.zero ) then
  S(4) = S(4)+p%weight
  return
else if(p%side.eq.4 .and. p%dir(3).gt.zero ) then
  S(6) = S(6)+p%weight
  return

else if(p%side.eq.5 .and. p%dir(3).le.zero ) then
  S(7) = S(7)+p%weight
  return
else if(p%side.eq.5 .and. p%dir(3).gt.zero ) then
  S(9) = S(7)+p%weight
  return

else if(p%side.eq.6 .and. p%dir(3).le.zero ) then
  S(8) = S(8)+p%weight
  return
else if(p%side.eq.6 .and. p%dir(3).gt.zero ) then
  S(10) =S(10)+p%weight
  return
else
  print *,'Couldnt find a stream on to which I can put the photon weight on?!'
  call print_photon(p)
  call exit
endif
select type(bmc)
end select

end subroutine

subroutine init_dir_photon_8_10(bmc,p,src,ldirect,initial_dir,vertices,ierr)
        class(t_boxmc_8_10)           :: bmc
        type(photon),intent(inout)    :: p
        real(ireal_dp),intent(in)     :: initial_dir(:), vertices(:)
        integer(iintegers),intent(in) :: src
        logical,intent(in)            :: ldirect
        integer(mpiint),intent(out)   :: ierr

        real(ireal_dp),parameter :: h=.5_ireal_dp
        real(ireal_dp),parameter :: e=1e-6_ireal_dp

        call setup_cube_coords_from_vertices(vertices, p%dx, p%dy, p%dz)

        p%alive = .False.

!          /     3       4 /
!         / 1       2     /
!        /_______________/
!        |               |
!      6 |           8   |
!        |               |
!        |               |
!        |    7          |
!      5 |_______________|
        select case(src)
        case (1)
                p%loc = [ L(p%dx)*h     , L(p%dy)*h     , p%dz*(one-e)   ]
        case (2)
                p%loc = [ p%dx*h+L(p%dx)*h, L(p%dy)*h     , p%dz*(one-e)   ]
        case (3)
                p%loc = [ L(p%dx)*h     , p%dy*h+L(p%dy)*h, p%dz*(one-e)   ]
        case (4)
                p%loc = [ p%dx*h+L(p%dx)*h, p%dy*h+L(p%dy)*h, p%dz*(one-e)   ]
        case (5)
                p%loc = [ zero+p%dx*e   , L(p%dy)       , L(p%dz)*h      ]
        case (6)
                p%loc = [ zero+p%dx*e   , L(p%dy)       , p%dz*h+L(p%dz)*h ]
        case (7)
                p%loc = [ L(p%dx)       , zero+p%dy*e   , L(p%dz)*h      ]
        case (8)
                p%loc = [ L(p%dx)       , zero+p%dy*e   , p%dz*h+L(p%dz)*h ]
        case default
                print *,'Dont know what to do with source spec:',src
                call exit
        end select

        p%weight=one
        p%alive = .True.
        p%direct= ldirect
        p%side = int(nil)
        p%src  = src
        p%scattercnt  = 0
        p%dir = initial_dir
        select type(bmc)
        end select

        ierr = 0
end subroutine

subroutine init_diff_photon_8_10(bmc,p,src,vertices,ierr)
        use m_helper_functions_dp, only : rotate_angle_x, rotate_angle_y, rotate_angle_z
        class(t_boxmc_8_10)           :: bmc
        type(photon),intent(inout)    :: p
        real(ireal_dp),intent(in)     :: vertices(:)
        integer(iintegers),intent(in) :: src
        integer(mpiint),intent(out)   :: ierr

        real(ireal_dp) :: phi,mu,mu1,mu2

        call setup_cube_coords_from_vertices(vertices, p%dx, p%dy, p%dz)

        mu1= sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
        mu2 = R() ! evenly distribution between 0 and 90 degree

        mu=mu1
        select case(src)
        case (-1) ! thermal emission
                p%loc = (/ L(p%dx) , L(p%dy) , L(p%dz) /)
                mu = 2*mu2 - one
                phi= deg2rad( R()*360 )
        case (1)
                p%loc = (/ L(p%dx) , L(p%dy) , zero  /)
                phi= deg2rad( R()*360 )
        case (2)
                p%loc = (/ L(p%dx) , L(p%dy) , p%dz    /)
                phi= deg2rad( R()*360 )
        case (3)
                p%loc = (/ p%dx    , L(p%dy) , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (4)
                p%loc = (/ zero  , L(p%dy) , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (5)
                p%loc = (/ p%dx    , L(p%dy) , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (6)
                p%loc = (/ zero  , L(p%dy) , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (7)
                p%loc = (/ L(p%dx) , p%dy    , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (8)
                p%loc = (/ L(p%dx) , zero  , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (9)
                p%loc = (/ L(p%dx) , p%dy    , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (10)
                p%loc = (/ L(p%dx) , zero  , L(p%dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case default
                print *,'Dont know what to do with source spec:',src
                call exit
        end select
        select case(src)
          case(1:2)
            phi= deg2rad( R()*360 )
          case(3:10)
            phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        end select

        p%dir = (/sin(phi)*sin(acos(mu)) , cos(phi)*sin(acos(mu)) , mu  /)
        ! phi = 0 is y-axes

        select case(src)
          case(2)
            p%dir = rotate_angle_y(p%dir,180*one)
          case(3)
            p%dir = rotate_angle_y(p%dir, 90*one)
            p%dir = rotate_angle_x(p%dir, 90*one)
          case(4)
            p%dir = rotate_angle_y(p%dir,-90*one)
            p%dir = rotate_angle_x(p%dir, 90*one)
          case(5)
            p%dir = rotate_angle_y(p%dir, 90*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(6)
            p%dir = rotate_angle_y(p%dir,-90*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(7)
            p%dir = rotate_angle_z(p%dir,180*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(8)
            p%dir = rotate_angle_x(p%dir, 90*one)

          case(9)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(10)
            p%dir = rotate_angle_z(p%dir,180*one)
            p%dir = rotate_angle_x(p%dir, 90*one)

        end select
        p%weight=one
        p%alive = .True.
        p%direct= .False.
        p%side = int(nil)
        p%src  = src
        p%scattercnt  = 0
        if(.not. approx( norm(p%dir) ,one ) ) then
          print *,'p%dir',p%dir,'norm',norm(p%dir)
          stop 'init_diff_photon_8_10 :: direction was not normed?!'
        endif
        select type(bmc)
        end select

        ierr = 0
end subroutine

subroutine intersect_distance_8_10(bmc,vertices,p,max_dist)
        class(t_boxmc_8_10) :: bmc
        real(ireal_dp),intent(in) :: vertices(:)
        type(photon),intent(inout) :: p
        real(ireal_dp),intent(out) :: max_dist

        call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
          p%side, max_dist)

        select type(bmc)
        end select

end subroutine
  ! ******************* BOXMC 8 10 ************************

