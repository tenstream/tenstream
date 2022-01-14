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

subroutine update_dir_stream_1_2(bmc,vertices,p,T)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(in) :: p
    real(irealbmc),intent(inout) :: T(:)

    select case (p%side)
    case(1:2)
      T(1) = T(1)+p%weight
      return
    case(3:6)
      call CHKWARN(1_mpiint, 'photon ended up on side face,'//&
        & 'this can happen but should be a very rare event!'//&
        & 'if you encounter this more often, dont trust the results!')
      T(1) = T(1)+p%weight
      return
    case default
      call CHKERR(1_mpiint, 'Dont know what to do with this p%side')
      print *,'remove unused compiler warning', vertices
    end select
    select type(bmc)
    end select
end subroutine

subroutine update_diff_stream_1_2(bmc,vertices,p,S)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(in) :: p
    real(irealbmc),intent(inout) :: S(:)

    !         _______1_______
    !        |               |
    !        |               |
    !        |               |
    !        |               |
    !        |               |
    !        |_______________|
    !                2
    select case (p%side)
    case(1)
      S(1) = S(1)+p%weight
      return
    case(2)
      S(2) = S(2)+p%weight
      return
    case default
      call print_photon(p)
      call CHKERR(1_mpiint, 'Dont know where to put photon weight')
      print *,'remove unused compiler warning', vertices
    end select
    select type(bmc)
    end select
end subroutine

subroutine init_dir_photon_1_2(bmc, p, src, ldirect, initial_dir, vertices, ierr)
    class(t_boxmc_1_2) :: bmc
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: ldirect
    integer(mpiint),intent(out) :: ierr

    real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

    p%dir   = initial_dir
    p%alive = .False.

    call CHKERR(size(vertices)-24, 'wrong size of vertices, expected 24 but found '//toStr(size(vertices)))
    associate( &
        pA => vertices(1:3), &
        pB => vertices(4:6), &
        pC => vertices(7:9), &
        pD => vertices(10:12), &
        pE => vertices(13:15), &
        pF => vertices(16:18), &
        pG => vertices(19:21), &
        pH => vertices(22:24) &
        )
      select case(src)
      case (1)
        if(p%dir(3).le.0) then
          call rand_pnt_on_plane(pE,pG,pH,pF, p%loc, normal, U, V)
          p%src_side = 1
        else
          call rand_pnt_on_plane(pA,pB,pD,pC, p%loc, normal, U, V)
          p%src_side = 2
        endif
      case default
        call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
      end select

      ! always nudge the photon a little bit inside the volume
      cube_center = (pA+pB+pC+pD+pE+pF+pG+pH)/8._ireal_dp
      dir_p_to_cube_center = cube_center - p%loc

      p%loc = p%loc + sqrt(epsilon(p%loc)) * dir_p_to_cube_center
    end associate

    p%weight     = one
    p%alive      = .True.
    p%direct     = ldirect
    p%side       = int(nil)
    p%src        = src
    p%scattercnt = 0

    select type(bmc)
    end select

    if(dot_product(p%dir, normal).lt.zero) then
      ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
      call CHKWARN(1_mpiint, 'initial direction of the photon is pointing outside of the cube')
      print *,'UUPS, are we starting in a wrong direction?'// &
        & ' :: init_dir_1_2', rad2deg(angle_between_two_vec(p%dir, normal))
      print *,'src', src, 'initial_dir', initial_dir
      print *,'normal',normal,'direction',p%dir
      call print_photon(p)
      call CHKERR(ierr)
      return
    endif

    ierr = 0
    !call print_photon(p)
    !call CHKERR(1_mpiint, 'DEBUG')
end subroutine

subroutine init_diff_photon_1_2(bmc, p, src, vertices, ierr)
    class(t_boxmc_1_2) :: bmc
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: vertices(:)
    integer(iintegers),intent(in) :: src
    integer(mpiint),intent(out) :: ierr

    real(ireal_dp) :: phi,mu
    real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

    call CHKERR(size(vertices)-24, 'wrong size of vertices, expected 24 but found '//toStr(size(vertices)))
    associate( &
        pA => vertices(1:3), &
        pB => vertices(4:6), &
        pC => vertices(7:9), &
        pD => vertices(10:12), &
        pE => vertices(13:15), &
        pF => vertices(16:18), &
        pG => vertices(19:21), &
        pH => vertices(22:24) &
        )

      select case(src)
      case (1)
        call rand_pnt_on_plane(pA,pB,pD,pC, p%loc, normal, U, V)
        p%src_side = 2
      case (2)
        call rand_pnt_on_plane(pE,pG,pH,pF, p%loc, normal, U, V)
        p%src_side = 1
      case default
        call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
      end select

      select case(src)
      case (1,2)
        phi= deg2rad( interv_R(-90._ireal_dp,0._ireal_dp) )
      case default
        phi = -999
        call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
      end select

      select case(src)
      case(1:2)
        mu = sqrt(interv_R(zero, one))
      case default
        mu = -1
        call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
      end select

      p%dir = (/sin(phi)*sin(acos(mu)) , cos(phi)*sin(acos(mu)) , mu  /)
      ! phi = 0 is y-axes

      select case(src)
      case(1)
        continue
      case(2)
        p%dir(3) = -p%dir(3)
      case default
        call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
      end select

      ! always nudge the photon a little bit inside the volume
      cube_center = (pA+pB+pC+pD+pE+pF+pG+pH)/8._ireal_dp
      dir_p_to_cube_center = cube_center - p%loc

      p%loc = p%loc + sqrt(epsilon(p%loc)) * dir_p_to_cube_center
    end associate

    p%weight=one
    p%alive = .True.
    p%direct= .False.
    p%side = int(nil)
    p%src  = src
    p%scattercnt  = 0

    if(.not. approx( norm2(p%dir) ,one ) ) then
      print *,'p%dir',p%dir,'norm',norm2(p%dir)
      stop 'init_diff_photon_3_24 :: direction was not normed?!'
    endif

    select type(bmc)
    end select

    if(dot_product(p%dir, normal).lt.zero) then
      ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
      call CHKWARN(1_mpiint, 'initial direction of the photon is pointing outside of the cube')
      print *,'UUPS, are we starting in a wrong direction?'// &
        & ' :: init_diff_photon_bmx_3_24', rad2deg(angle_between_two_vec(p%dir, normal))
      print *,'src', src
      print *,'normal',normal,'direction',p%dir
      call print_photon(p)
      call CHKERR(ierr)
      return
    endif

    ierr = 0
    !call print_photon(p)
    !call CHKERR(1_mpiint, 'DEBUG')
  end subroutine

subroutine intersect_distance_1_2(bmc,vertices,p,max_dist)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(out) :: max_dist

    real(ireal_dp) :: hit(2)
    integer(iintegers) :: i

    associate( &
        pA => vertices(1:3), &
        pB => vertices(4:6), &
        pC => vertices(7:9), &
        pD => vertices(10:12), &
        pE => vertices(13:15), &
        pF => vertices(16:18), &
        pG => vertices(19:21), &
        pH => vertices(22:24) &
        )
      p%loc(1:2) = pA(1:2) ! horizontally move the photon back to the vicinity of the initial vertices (reduces numeric errors)

      hit(1) = hit_plane(p%loc, p%dir, pE, cross_3d(pF-pE, pG-pE))
      hit(2) = hit_plane(p%loc, p%dir, pA, cross_3d(pC-pA, pB-pA))
    end associate

    p%side = 0
    max_dist = huge(max_dist)
    do i = 1, 2
      if(hit(i).lt.zero) cycle
      if(p%scattercnt.eq.0 .and. i.eq.p%src_side) cycle
      if(hit(i).lt.max_dist) then
        max_dist = hit(i)
        p%side    = i
        p%subface = 1
      endif
    enddo
    if(p%side.eq.0) then
      call print_photon(p)
      call CHKERR(1_mpiint, 'did not hit anything')
    endif
    select type(bmc)
    end select
end subroutine

subroutine box_halfspaces_1_2(bmc, vertices, origins, normals)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    real(ireal_dp),allocatable,intent(out) :: origins(:,:), normals(:,:)
    call box_halfspaces(vertices, origins, normals)
    select type(bmc)
    end select
end subroutine
