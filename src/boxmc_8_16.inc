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


subroutine update_dir_stream_8_16(bmc,vertices,p,T)
        class(t_boxmc_8_16) :: bmc
        type(t_photon),intent(in) :: p
        real(ireal_dp),intent(in)     :: vertices(:)
        real(ireal_dp),intent(inout) :: T(:)
        real(ireal_dp),parameter :: h=.5_ireal_dp
        real(ireal_dp) :: dx, dy, dz

!          /     3       4 /
!         / 1       2     /
!        /_______________/
!        |               |
!      6 |           8   |
!        |               |
!        |               |
!        |    7          |
!      5 |_______________|

    call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

    select case (p%side)
    case(1:2)
      if(p%loc(1).le.dx*h.and.p%loc(2).le.dy*h) then
        T(1) = T(1)+p%weight
        return
      else if(p%loc(1).gt.dx*h.and.p%loc(2).le.dy*h) then
        T(2) = T(2)+p%weight
        return
      else if(p%loc(1).le.dx*h.and.p%loc(2).gt.dy*h) then
        T(3) = T(3)+p%weight
        return
      else if(p%loc(1).gt.dx*h.and.p%loc(2).gt.dy*h) then
        T(4) = T(4)+p%weight
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif

    case(3:4)

      if(p%loc(3).le.dz*h ) then
        T(5) = T(5)+p%weight
        return
      else if(p%loc(3).gt.dz*h ) then
        T(6) = T(6)+p%weight
        return
      else
        print *,'Couldnt find a stream on to which I can put the photon weight on?!'
        call print_photon(p)
      endif
    case(5:6)

      if(p%loc(3).le.dz*h ) then
        T(7) = T(7)+p%weight
        return
      else if(p%loc(3).gt.dz*h ) then
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
subroutine update_diff_stream_8_16(bmc,vertices,p,S)
        class(t_boxmc_8_16) :: bmc
        real(ireal_dp),intent(in) :: vertices(:)
        type(t_photon),intent(in) :: p
        real(ireal_dp),intent(inout) :: S(:)
        integer(iintegers) :: imu
        real(ireal_dp) :: mu
        ! rad2deg(arccos(mu_smpl_pts)) =  [ 0, 41, 60, 75, 90 ]
        real(ireal_dp), parameter :: mu_smpl_pts(5)=[1._ireal_dp, .75_ireal_dp, .50_ireal_dp, .25_ireal_dp, 0._ireal_dp]
        integer(iintegers),parameter :: soff=(size(mu_smpl_pts)-1)*2 ! offset for sideward streams

!         _____1_3_5_____
!        |           14  |
!      9 |           12  | 10
!        |               |
!        |   13          |
!        |   11          |
!      7 |_______________| 8
!              2 4 6

if(p%side.eq.1) then
  mu = dot_product(p%dir, [zero, zero, one])
  if(mu.lt.0) then
    call print_photon(p)
    call CHKERR(1_mpiint, 'wrong p%side?')
  endif
  do imu=1,size(mu_smpl_pts)-1
    if(mu.le.mu_smpl_pts(imu) .and. mu.ge.mu_smpl_pts(imu+1)) then
      S(imu*2-1) = S(imu*2-1)+p%weight
      exit
    endif
  enddo
  return

else if(p%side.eq.2) then
  mu = dot_product(p%dir, [zero, zero, -one])
  if(mu.lt.0) then
    call print_photon(p)
    call CHKERR(1_mpiint, 'wrong p%side?')
  endif
  do imu=1,size(mu_smpl_pts)-1
    if(mu.le.mu_smpl_pts(imu) .and. mu.ge.mu_smpl_pts(imu+1)) then
      S(imu*2) = S(imu*2)+p%weight
      exit
    endif
  enddo
  return

else if(p%side.eq.3 .and. p%dir(3).le.zero ) then
  S(soff+1) = S(soff+1)+p%weight
  return
else if(p%side.eq.3 .and. p%dir(3).gt.zero ) then
  S(soff+3) = S(soff+3)+p%weight
  return

else if(p%side.eq.4 .and. p%dir(3).le.zero ) then
  S(soff+2) = S(soff+2)+p%weight
  return
else if(p%side.eq.4 .and. p%dir(3).gt.zero ) then
  S(soff+4) = S(soff+4)+p%weight
  return

else if(p%side.eq.5 .and. p%dir(3).le.zero ) then
  S(soff+5) = S(soff+5)+p%weight
  return
else if(p%side.eq.5 .and. p%dir(3).gt.zero ) then
  S(soff+7) = S(soff+7)+p%weight
  return

else if(p%side.eq.6 .and. p%dir(3).le.zero ) then
  S(soff+6) = S(soff+6)+p%weight
  return
else if(p%side.eq.6 .and. p%dir(3).gt.zero ) then
  S(soff+8) =S(soff+8)+p%weight
  return
else
  print *,'Couldnt find a stream on to which I can put the photon weight on?!'
  call print_photon(p)
  call exit
  print *,'remove unused compiler warning', vertices
endif
select type(bmc)
end select

end subroutine

subroutine init_dir_photon_8_16(bmc,p,src,ldirect,initial_dir,vertices,ierr)
        class(t_boxmc_8_16)           :: bmc
        type(t_photon),intent(inout)    :: p
        real(ireal_dp),intent(in)     :: initial_dir(:), vertices(:)
        integer(iintegers),intent(in) :: src
        logical,intent(in)            :: ldirect
        integer(mpiint),intent(out)   :: ierr

        real(ireal_dp),parameter :: h=.5_ireal_dp
        real(ireal_dp),parameter :: e=1e-6_ireal_dp
        real(ireal_dp) :: dx, dy, dz

        call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

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
                p%loc = [ L(dx)*h     , L(dy)*h     , dz*(one-e)   ]
        case (2)
                p%loc = [ dx*h+L(dx)*h, L(dy)*h     , dz*(one-e)   ]
        case (3)
                p%loc = [ L(dx)*h     , dy*h+L(dy)*h, dz*(one-e)   ]
        case (4)
                p%loc = [ dx*h+L(dx)*h, dy*h+L(dy)*h, dz*(one-e)   ]
        case (5)
                p%loc = [ zero+dx*e   , L(dy)       , L(dz)*h      ]
        case (6)
                p%loc = [ zero+dx*e   , L(dy)       , dz*h+L(dz)*h ]
        case (7)
                p%loc = [ L(dx)       , zero+dy*e   , L(dz)*h      ]
        case (8)
                p%loc = [ L(dx)       , zero+dy*e   , dz*h+L(dz)*h ]
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

subroutine init_diff_photon_8_16(bmc,p,src,vertices,ierr)
        use m_helper_functions_dp, only : rotate_angle_x, rotate_angle_y, rotate_angle_z
        class(t_boxmc_8_16)           :: bmc
        type(t_photon),intent(inout)    :: p
        real(ireal_dp),intent(in)     :: vertices(:)
        integer(iintegers),intent(in) :: src
        integer(mpiint),intent(out)   :: ierr

        real(ireal_dp) :: phi,mu
        real(ireal_dp) :: dx, dy, dz
        real(ireal_dp), parameter :: mu_smpl_pts(5)=[1._ireal_dp, .75_ireal_dp, .50_ireal_dp, .25_ireal_dp, 0._ireal_dp] ** 2
        integer(iintegers), parameter :: soff=(size(mu_smpl_pts)-1)*2

        call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

        mu= sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian

        select case(src)
        case (1)
                p%loc = (/ L(dx) , L(dy) , zero  /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(2), mu_smpl_pts(1)))
        case (2)
                p%loc = (/ L(dx) , L(dy) , dz    /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(2), mu_smpl_pts(1)))
        case (3)
                p%loc = (/ L(dx) , L(dy) , zero  /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(3), mu_smpl_pts(2)))
        case (4)
                p%loc = (/ L(dx) , L(dy) , dz    /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(3), mu_smpl_pts(2)))
        case (5)
                p%loc = (/ L(dx) , L(dy) , zero  /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(4), mu_smpl_pts(3)))
        case (6)
                p%loc = (/ L(dx) , L(dy) , dz    /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(4), mu_smpl_pts(3)))
        case (7)
                p%loc = (/ L(dx) , L(dy) , zero  /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(5), mu_smpl_pts(4)))
        case (8)
                p%loc = (/ L(dx) , L(dy) , dz    /)
                phi= deg2rad( R()*360 )
                mu = sqrt(interv_R(mu_smpl_pts(5), mu_smpl_pts(4)))
        case (soff+1)
                p%loc = (/ dx    , L(dy) , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+2)
                p%loc = (/ zero  , L(dy) , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+3)
                p%loc = (/ dx    , L(dy) , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+4)
                p%loc = (/ zero  , L(dy) , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+5)
                p%loc = (/ L(dx) , dy    , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+6)
                p%loc = (/ L(dx) , zero  , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+7)
                p%loc = (/ L(dx) , dy    , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case (soff+8)
                p%loc = (/ L(dx) , zero  , L(dz) /)
                phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        case default
                print *,'Dont know what to do with source spec:',src
                call exit
        end select
        select case(src)
          case(1:soff)
            phi= deg2rad( R()*360 )
          case(soff+1:soff+8)
            phi= deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
        end select

        p%dir = (/sin(phi)*sin(acos(mu)) , cos(phi)*sin(acos(mu)) , mu  /)

        select case(src)
          case(1,3,5,7)
            continue
          case(2,4,6,8)
            p%dir = rotate_angle_y(p%dir,180*one)
          case(soff+1)
            p%dir = rotate_angle_y(p%dir, 90*one)
            p%dir = rotate_angle_x(p%dir, 90*one)
          case(soff+2)
            p%dir = rotate_angle_y(p%dir,-90*one)
            p%dir = rotate_angle_x(p%dir, 90*one)
          case(soff+3)
            p%dir = rotate_angle_y(p%dir, 90*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(soff+4)
            p%dir = rotate_angle_y(p%dir,-90*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(soff+5)
            p%dir = rotate_angle_z(p%dir,180*one)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(soff+6)
            p%dir = rotate_angle_x(p%dir, 90*one)
          case(soff+7)
            p%dir = rotate_angle_x(p%dir,-90*one)
          case(soff+8)
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
          stop 'init_diff_photon_8_16 :: direction was not normed?!'
        endif
        select type(bmc)
        end select

        ierr = 0
end subroutine

subroutine intersect_distance_8_16(bmc,vertices,p,max_dist)
        class(t_boxmc_8_16) :: bmc
        real(ireal_dp),intent(in) :: vertices(:)
        type(t_photon),intent(inout) :: p
        real(ireal_dp),intent(out) :: max_dist

        call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
          p%side, max_dist)

        select type(bmc)
        end select

end subroutine
  ! ******************* BOXMC 8 10 ************************
