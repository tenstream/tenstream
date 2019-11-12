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

! Distribution Code for direct streams/targets:
!
!   C ------- D
!   |         |
!  2|   1     |
!   |         |
!   |         |
!   A ------- B
!        3
! Labeling for destinations is symmetric horizontally and vertically
!
! Face numbering is
! 1 -> top/bot face
! 2 -> side face between (A,C)/(A,D)
! 3 -> side face between (A,B)/(C,D)

! Distribution Code for diffuse streams/targets:
!         _______2_______
!        |            9  |
!      6 |            7  | 5
!        |               |
!        |   10          |
!        |   8           |
!      4 |_______________| 3
!                1
!
! Face numbering is
! 2 -> top
! 1 -> bottom face
! 3,4 -> side faces in x-plane downwards
! 5,6 -> side faces in x-plane upwards
! 7,8 -> side faces in y-plane downwards
! 9,10-> side faces in y-plane upwards

subroutine init_dir_photon_3_10(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_3_10)               :: bmc
  type(t_photon)    ,intent(inout)  :: p
  real(ireal_dp)    ,intent(in)     :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in)     :: src
  logical           ,intent(in)     :: ldirect
  integer(mpiint)   ,intent(out)    :: ierr

  real(ireal_dp), parameter  :: e=1e-6_ireal_dp
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  p%alive = .False.

  select case(src)
  case (1)
    p%loc = [R()*dx, R()*dy, dz*(one-e)]
  case (2)
    p%loc = [e, R()*dy, R()*dz]
  case (3)
    p%loc = [R()*dx, e, R()*dz]
  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  p%weight     = one
  p%alive      = .True.
  p%direct     = ldirect
  p%side       = int(nil)
  p%src        = src
  p%scattercnt = 0
  p%dir        = initial_dir
  select type(bmc)
  end select

  ierr = 0
end subroutine

subroutine update_dir_stream_3_10(bmc,vertices,p,T)
        class(t_boxmc_3_10)          :: bmc
        real(ireal_dp),intent(in)    :: vertices(:)
        type(t_photon),intent(in)    :: p
        real(irealbmc),intent(inout) :: T(:)
        select case(p%side)
        case(1,2)
          T(1)=T(1)+p%weight
        case(3,4)
          T(2)=T(2)+p%weight
        case(5,6)
          T(3)=T(3)+p%weight
        case default
          stop 'called update dirstream with incorrect source'
          print *,'remove unused compiler warning', vertices
        end select
        select type(bmc)
        end select
end subroutine

subroutine init_diff_photon_3_10(bmc, p, src, vertices, ierr)
        use m_helper_functions, only : rotate_angle_x, rotate_angle_y, rotate_angle_z
        class(t_boxmc_3_10)           :: bmc
        type(t_photon),intent(inout)    :: p
        real(ireal_dp),intent(in)     :: vertices(:)
        integer(iintegers),intent(in) :: src
        integer(mpiint),intent(out)   :: ierr

        real(ireal_dp) :: phi,mu,mu1,mu2
        real(ireal_dp) :: dx, dy, dz

        call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

        mu1= sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
        mu2 = R() ! evenly distribution between 0 and 90 degree

        mu=mu1
        select case(src)
        case (-1) ! thermal emission
                p%loc = (/ L(dx) , L(dy) , L(dz) /)
                mu = 2*mu2 - one
        case (1)
                p%loc = (/ L(dx) , L(dy) , zero  /)
        case (2)
                p%loc = (/ L(dx) , L(dy) , dz    /)
        case (3)
                p%loc = (/ dx    , L(dy) , L(dz) /)
        case (4)
                p%loc = (/ zero  , L(dy) , L(dz) /)
        case (5)
                p%loc = (/ dx    , L(dy) , L(dz) /)
        case (6)
                p%loc = (/ zero  , L(dy) , L(dz) /)
        case (7)
                p%loc = (/ L(dx) , dy    , L(dz) /)
        case (8)
                p%loc = (/ L(dx) , zero  , L(dz) /)
        case (9)
                p%loc = (/ L(dx) , dy    , L(dz) /)
        case (10)
                p%loc = (/ L(dx) , zero  , L(dz) /)
        case default
                call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
        end select

        select case(src)
          case(-1,1,2)
            phi = deg2rad( R()*360 )
          case(3:10)
            phi = deg2rad( interv_R(-90._ireal_dp,90._ireal_dp) )
          case default
            call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
            call exit
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

        if(.not. approx( norm2(p%dir) ,one ) ) then
          print *,'p%dir',p%dir,'norm',norm2(p%dir)
          stop 'init_diff_photon_3_10 :: direction was not normed?!'
        endif

        select type(bmc)
        end select

        ierr = 0
end subroutine

subroutine update_diff_stream_3_10(bmc,vertices,p,S)
  class(t_boxmc_3_10)          :: bmc
  real(ireal_dp),intent(in)    :: vertices(:)
  type(t_photon),intent(in)    :: p
  real(irealbmc),intent(inout) :: S(:)

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
    S(9) = S(9)+p%weight
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
    print *,'remove unused compiler warning', vertices
  endif
  select type(bmc)
  end select

end subroutine

subroutine intersect_distance_3_10(bmc,vertices,p,max_dist)
  class(t_boxmc_3_10) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(out) :: max_dist


  call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    p%side, max_dist)

  select type(bmc)
  end select
end subroutine
subroutine box_halfspaces_3_10(bmc, vertices, origins, normals)
    class(t_boxmc_3_10) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    real(ireal_dp),allocatable,intent(out) :: origins(:,:), normals(:,:)
    call box_halfspaces(vertices, origins, normals)
    select type(bmc)
    end select
end subroutine
  ! ******************* BOXMC 3 10 ************************

