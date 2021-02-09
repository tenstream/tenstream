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

module m_optprop
#include "petsc/finclude/petsc.h"
  use petsc

#ifdef _XLF
      use ieee_arithmetic
#define isnan ieee_is_nan
#endif

use m_optprop_parameters, only : ldebug_optprop, wedge_sphere_radius, param_eps
use m_helper_functions, only : rmse, CHKERR, CHKWARN, toStr, cstr, approx, deg2rad, rad2deg, swap, is_between, vec_proj_on_plane, &
  char_arr_to_str, triangle_area_by_vertices, compute_normal_3d, volume_hexahedron
use m_data_parameters, only: ireals,ireal_dp,irealLUT,ireal_params,iintegers,one,zero,i0,i1,inil,mpiint
use m_optprop_base, only: t_optprop_base, t_op_config, find_op_dim_by_name
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_3_6, t_optprop_LUT_3_10, &
  t_optprop_LUT_8_10, t_optprop_LUT_3_16, t_optprop_LUT_8_16, t_optprop_LUT_8_18, &
  t_optprop_LUT_wedge_5_8, t_optprop_LUT_rectilinear_wedge_5_8, t_optprop_LUT_wedge_18_8

use m_optprop_ANN, only : t_optprop_ANN, t_optprop_ANN_3_10
use m_boxmc_geometry, only : setup_default_unit_cube_geometry, setup_default_wedge_geometry
use m_eddington, only: eddington_coeff_zdun
use m_tenstream_options, only: twostr_ratio
use m_intersection, only: hit_plane, line_intersection_3d

use m_LUT_param_phi, only: theta_from_param_theta, iterative_phi_theta_from_param_phi_and_param_theta

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public ::                          &
  t_optprop,                       &
  t_optprop_cube,                  &
  t_optprop_wedge,                 &
  t_optprop_1_2,                   &
  t_optprop_3_6,                   &
  t_optprop_3_10,                  &
  t_optprop_3_10_ann,              &
  t_optprop_8_10,                  &
  t_optprop_3_16,                  &
  t_optprop_8_16,                  &
  t_optprop_8_18,                  &
  t_optprop_wedge_5_8,             &
  t_optprop_rectilinear_wedge_5_8, &
  t_optprop_wedge_18_8,            &
  OPP_1D_RETCODE,                  &
  OPP_TINYASPECT_RETCODE,          &
  dir2dir3_coeff_corr_zx,          &
  dir2dir3_coeff_corr_zy,          &
  dir2dir3_coeff_corr_xx,          &
  dir2dir3_coeff_corr_xy,          &
  dir2dir3_coeff_corr_yy,          &
  dir2dir3_coeff_corr

type,abstract :: t_optprop
  logical :: optprop_debug=ldebug_optprop
  class(t_optprop_LUT), allocatable :: LUT
  class(t_optprop_ANN), allocatable :: ANN
  class(t_optprop_base), pointer :: dev
  contains
    procedure :: init
    procedure :: get_coeff_bmc
    procedure :: destroy
end type

! Cube types
type,abstract,extends(t_optprop) :: t_optprop_cube
  contains
    procedure :: get_coeff => get_coeff_cube
    procedure :: dir2dir_coeff_symmetry => dir2dir_coeff_symmetry_none
    procedure :: dir2diff_coeff_symmetry => dir2diff_coeff_symmetry_none
    procedure :: diff2diff_coeff_symmetry => diff2diff_coeff_symmetry
end type

! we introduce one special cube type for 8 direct streams, this way, all of them can share dir2dir_coeff_symmetry
type,abstract,extends(t_optprop_cube) :: t_optprop_cube_dir8
  contains
    procedure :: dir2dir_coeff_symmetry => dir2dir8_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_1_2
end type

type,extends(t_optprop_cube) :: t_optprop_3_6
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff6_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_3_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff10_coeff_symmetry
end type

type,extends(t_optprop_cube) :: t_optprop_3_10_ann
end type

type,extends(t_optprop_cube) :: t_optprop_3_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff16_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff10_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff16_coeff_symmetry
end type

type,extends(t_optprop_cube_dir8) :: t_optprop_8_18
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff18_coeff_symmetry
end type

! Wedge types
type,abstract,extends(t_optprop) :: t_optprop_wedge
  contains
    procedure :: get_coeff => get_coeff_wedge
end type

type,extends(t_optprop_wedge) :: t_optprop_wedge_5_8
end type

type,extends(t_optprop_wedge) :: t_optprop_rectilinear_wedge_5_8
end type

type,extends(t_optprop_wedge) :: t_optprop_wedge_18_8
end type

integer(mpiint), parameter :: OPP_1D_RETCODE = -1_mpiint
integer(mpiint), parameter :: OPP_TINYASPECT_RETCODE = -2_mpiint

logical, parameter :: lDEBUG_geometric_coeff_correction = .true.

contains

  subroutine init(OPP, comm, skip_load_LUT)
      class(t_optprop), target, intent(inout) :: OPP
      integer(mpiint) ,intent(in) :: comm
      logical, intent(in), optional :: skip_load_LUT
      integer(mpiint) :: ierr

      select type(OPP)
      class is (t_optprop_1_2)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_1_2::OPP%LUT)

      class is (t_optprop_3_6)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_3_6::OPP%LUT)

      class is (t_optprop_3_10)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_3_10::OPP%LUT)

      class is (t_optprop_8_10)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_8_10::OPP%LUT)

      class is (t_optprop_3_16)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_3_16::OPP%LUT)

      class is (t_optprop_8_16)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_8_16::OPP%LUT)

      class is (t_optprop_8_18)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_8_18::OPP%LUT)

      class is (t_optprop_wedge_5_8)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_wedge_5_8::OPP%LUT)

      class is (t_optprop_rectilinear_wedge_5_8)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_rectilinear_wedge_5_8::OPP%LUT)

      class is (t_optprop_wedge_18_8)
        if(.not.allocated(OPP%LUT) ) allocate(t_optprop_LUT_wedge_18_8::OPP%LUT)

      class is (t_optprop_3_10_ann)
        if(.not.allocated(OPP%ANN) ) allocate(t_optprop_ANN_3_10::OPP%ANN)

      class default
        call CHKERR(1_mpiint, ' init optprop : unexpected type for optprop object!')
      end select

      if(allocated(OPP%LUT)) then
        call OPP%LUT%init(comm, skip_load_LUT)
        OPP%dev => OPP%LUT
      endif
      if(allocated(OPP%ANN)) then
        call OPP%ANN%init(comm, ierr); call CHKERR(ierr)
        OPP%dev => OPP%ANN
      endif

  end subroutine
  subroutine destroy(OPP)
      class(t_optprop) :: OPP
      integer(mpiint) :: ierr
      if(allocated(OPP%LUT)) then
          call OPP%LUT%destroy(ierr); call CHKERR(ierr)
          deallocate(OPP%LUT)
      endif
      if(allocated(OPP%ANN)) then
          call OPP%ANN%destroy(ierr); call CHKERR(ierr)
          deallocate(OPP%ANN)
      endif
  end subroutine

  subroutine get_coeff_wedge(OPP, tauz, w0, g, aspect_zx, ldir, C, ierr, wedge_coords, angles)
    class(t_optprop_wedge)              :: OPP
    logical,intent(in)                  :: ldir
    real(irealLUT),intent(in)           :: tauz, w0, g, aspect_zx
    real(irealLUT),intent(in)           :: wedge_coords(:) ! 2 coordinates of wedge C_point, only used for wedge OPP types
    real(irealLUT),intent(in),optional  :: angles(:)
    real(irealLUT),intent(out)          :: C(:)
    integer(mpiint), intent(out) :: ierr

    logical,parameter :: compute_coeff_online=.False.

    ierr = 0

    if(ldebug_optprop) then
      call check_inp(OPP, tauz, w0, g, aspect_zx, ldir, C, angles, wedge_coords)
    endif

    if(handle_aspect_zx_1D_case()) return

    if(compute_coeff_online) then
      call do_bmc_computation(C)
      return
    endif

    call do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)

    if(.False. .and. ldir) call print_coeff_diff()

    !if(ldir .and. present(angles)) then
    !  call handle_critical_azimuth()
    !endif
    contains
      subroutine do_bmc_computation(Cbmc)
        real(irealLUT), intent(out) :: Cbmc(:)
        real(ireals), allocatable :: vertices(:)
        real(ireal_params) :: phi, theta

        call setup_default_wedge_geometry(&
          real([0,0], ireals), &
          real([1,0], ireals), &
          real(wedge_coords, ireals), &
          real(aspect_zx, ireals), &
          vertices, &
          real(wedge_sphere_radius, ireals))

        if(present(angles)) then
          call iterative_phi_theta_from_param_phi_and_param_theta(&
            real(vertices, ireal_params), &
            real(angles(1), ireal_params), &
            real(angles(2), ireal_params), &
            phi, theta, ierr); call CHKERR(ierr)

          call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), ldir, Cbmc, &
            [rad2deg(real(phi, irealLUT)), rad2deg(real(theta, irealLUT))])
          print *,'Cbmc', tauz, w0, g, aspect_zx, wedge_coords, ':', angles, rad2deg(phi), rad2deg(theta), '=>', new_line(':'), Cbmc
        else
          call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), ldir, Cbmc)
          print *,'Cbmc', tauz, w0, g, aspect_zx, wedge_coords, '=>', new_line(':'), Cbmc
        endif

      end subroutine
      subroutine print_coeff_diff()
        real(irealLUT) :: Cbmc(size(C))
        real(ireals) :: err(2)
        integer(iintegers) :: isrc

        call do_bmc_computation(Cbmc)

        err = rmse(real(C, ireals), real(Cbmc, ireals))
        print *,'rmse', err
        do isrc=1,OPP%dev%dir_streams
          print *, 'lut src', isrc, ':', C(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams)
          print *, 'bmc src', isrc, ':', Cbmc(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams)
        enddo
        if(err(2).gt.one) then
          !call CHKERR(1_mpiint, 'DEBUG')
        endif
        !C = Cbmc
      end subroutine
      subroutine do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)
        real(irealLUT), intent(in) :: tauz, w0, aspect_zx
        logical,intent(in)       :: ldir
        real(irealLUT),intent(in),optional :: angles(:)
        real(irealLUT) :: save_param_phi, save_param_theta

        real(irealLUT), save, allocatable :: inp_arr_dir(:), inp_arr_diff(:)
        integer(iintegers), save :: dimidx_dir(8)
        integer(iintegers), save :: dimidx_diff(6)
        logical, save :: linit=.False.

        if(.not.linit) then
          dimidx_dir(1) = find_op_dim_by_name(OPP%dev%dirconfig, 'tau')
          dimidx_dir(2) = find_op_dim_by_name(OPP%dev%dirconfig, 'w0')
          dimidx_dir(3) = find_op_dim_by_name(OPP%dev%dirconfig, 'aspect_zx')
          dimidx_dir(4) = find_op_dim_by_name(OPP%dev%dirconfig, 'g')
          dimidx_dir(5) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cx')
          dimidx_dir(6) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cy')
          dimidx_dir(7) = find_op_dim_by_name(OPP%dev%dirconfig, 'param_phi')
          dimidx_dir(8) = find_op_dim_by_name(OPP%dev%dirconfig, 'param_theta')
          allocate(inp_arr_dir(count(dimidx_dir.gt.0)))

          dimidx_diff(1) = find_op_dim_by_name(OPP%dev%dirconfig, 'tau')
          dimidx_diff(2) = find_op_dim_by_name(OPP%dev%dirconfig, 'w0')
          dimidx_diff(3) = find_op_dim_by_name(OPP%dev%dirconfig, 'aspect_zx')
          dimidx_diff(4) = find_op_dim_by_name(OPP%dev%dirconfig, 'g')
          dimidx_diff(5) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cx')
          dimidx_diff(6) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cy')
          allocate(inp_arr_diff(count(dimidx_diff.gt.0)))
          linit = .True.
        endif

        if(present(angles)) then ! obviously we want the direct coefficients

          associate(&
              param_phi => angles(1), &
              param_theta => angles(2) )

            call handle_critical_param_phi(param_phi, save_param_phi)
            call handle_critical_param_theta(param_theta, save_param_theta)

            if(dimidx_dir(1).gt.0) inp_arr_dir(dimidx_dir(1)) = tauz
            if(dimidx_dir(2).gt.0) inp_arr_dir(dimidx_dir(2)) = w0
            if(dimidx_dir(3).gt.0) inp_arr_dir(dimidx_dir(3)) = aspect_zx
            if(dimidx_dir(4).gt.0) inp_arr_dir(dimidx_dir(4)) = g
            if(dimidx_dir(5).gt.0) inp_arr_dir(dimidx_dir(5)) = wedge_coords(1)
            if(dimidx_dir(6).gt.0) inp_arr_dir(dimidx_dir(6)) = wedge_coords(2)
            if(dimidx_dir(7).gt.0) inp_arr_dir(dimidx_dir(7)) = save_param_phi
            if(dimidx_dir(8).gt.0) inp_arr_dir(dimidx_dir(8)) = save_param_theta

            if(ldir) then ! dir2dir
              call OPP%dev%get_dir2dir(inp_arr_dir, C)
            else ! dir2diff
              call OPP%dev%get_dir2diff(inp_arr_dir, C)
            endif
          end associate
        else
          ! diff2diff
          if(dimidx_diff(1).gt.0) inp_arr_diff(dimidx_diff(1)) = tauz
          if(dimidx_diff(2).gt.0) inp_arr_diff(dimidx_diff(2)) = w0
          if(dimidx_diff(3).gt.0) inp_arr_diff(dimidx_diff(3)) = aspect_zx
          if(dimidx_diff(4).gt.0) inp_arr_diff(dimidx_diff(4)) = g
          if(dimidx_diff(5).gt.0) inp_arr_diff(dimidx_diff(5)) = wedge_coords(1)
          if(dimidx_diff(6).gt.0) inp_arr_diff(dimidx_diff(6)) = wedge_coords(2)
          call OPP%dev%get_diff2diff(inp_arr_diff, C)
        endif

      end subroutine

      subroutine handle_critical_param_phi(param_phi, save_param_phi)
        real(irealLUT), intent(in) :: param_phi
        real(irealLUT), intent(out) :: save_param_phi
        logical :: lsample_critical

        lsample_critical = .False.
        if(approx(abs(param_phi), 1._irealLUT, param_eps)) lsample_critical = .True.

        if(lsample_critical) then
          if(param_phi.lt.-1._irealLUT) then
            save_param_phi = -1._irealLUT-param_eps
          elseif(param_phi.ge.1._irealLUT) then !1.0001
            save_param_phi = 1._irealLUT+param_eps
          elseif(param_phi.lt.0._irealLUT) then !-.999
            save_param_phi = -1._irealLUT+param_eps
          else ! .999
            save_param_phi = 1._irealLUT-param_eps
          endif
        else
          save_param_phi = param_phi
        endif
      end subroutine
      subroutine handle_critical_param_theta(param_theta, save_param_theta)
        real(irealLUT), intent(in) :: param_theta
        real(irealLUT), intent(out) :: save_param_theta
        logical :: lsample_critical

        lsample_critical = .False.
        if(approx(abs(param_theta), 0._irealLUT, param_eps)) lsample_critical = .True.

        if(lsample_critical) then
          if(param_theta.le.0._irealLUT) then
            save_param_theta = -param_eps
          else ! .0001
            save_param_theta = param_eps
          endif
        else
          save_param_theta = param_theta
        endif
      end subroutine

      logical function handle_aspect_zx_1D_case()
        real(ireals) :: c11,c12,c13,c23,c33
        real(irealLUT) :: restricted_aspect_zx
        real(irealLUT) :: mu

        handle_aspect_zx_1D_case = .False.

        !TODO: here we may have incoming radiation at the sides and we just drop that
        ! this has to be fixed for anisotropic grids

        if(present(angles)) then

          if(aspect_zx.ge.twostr_ratio) then
            C = zero
            mu = real( cos(theta_from_param_theta(real(angles(2), ireal_params), 0._ireal_params)), irealLUT)

            call eddington_coeff_zdun(&
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              real(mu, ireals), &
              c11,c12,c13,c23,c33)

            if(ldir) then
              select type(OPP)
              class is (t_optprop_wedge_5_8)
                ! set the transport coeffs for src top to zero, leave the rest.
                C(5*4+1) = real(c33, irealLUT) ! from top to bot
                !C(22:24) = 1 ! from sides to bot
              class is (t_optprop_rectilinear_wedge_5_8)
                ! set the transport coeffs for src top to zero, leave the rest.
                C(5*4+1) = real(c33, irealLUT) ! from top to bot
                !C(22:24) = 1 ! from sides to bot
              class is (t_optprop_wedge_18_8)
                C(18*15+1) = real(c33, irealLUT) ! from top to bot
                C(18*16+2) = real(c33, irealLUT) ! from top to bot
                C(18*17+3) = real(c33, irealLUT) ! from top to bot
              class default
                call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
              end select

            else
              select type(OPP)
              class is (t_optprop_wedge_5_8)
                C(0*5+1) = real(c13, irealLUT) ! reflection
                C(7*5+1) = real(c23, irealLUT) ! transmission
              class is (t_optprop_rectilinear_wedge_5_8)
                C(0*5+1) = real(c13, irealLUT) ! reflection
                C(7*5+1) = real(c23, irealLUT) ! transmission
              class is (t_optprop_wedge_18_8)
                C(0*18+1) = real(c13, irealLUT)
                C(0*18+2) = real(c13, irealLUT)
                C(0*18+3) = real(c13, irealLUT)
                C(7*18+1) = real(c23, irealLUT)
                C(7*18+2) = real(c23, irealLUT)
                C(7*18+3) = real(c23, irealLUT)
              class default
                call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
              end select
            endif

            handle_aspect_zx_1D_case = .True.
            ierr = OPP_1D_RETCODE

          elseif(aspect_zx.lt.OPP%dev%dirconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%dev%dirconfig%dims(3)%vrange(1)), &
              OPP%dev%dirconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        else ! diffuse

          !if(aspect_zx.gt.OPP%dev%diffconfig%dims(3)%vrange(2)) then
          if(aspect_zx.ge.twostr_ratio) then
            C = zero

            call eddington_coeff_zdun(&
              real(tauz, ireals), &
              real(w0, ireals), &
              real(g, ireals), &
              one, &
              c11,c12,c13,c23,c33)

            ! transmission & reflection towards top plate
            C(1)     = real(c12, irealLUT)
            C(8)     = real(c11, irealLUT)
            ! and bot plate
            C(7*8+1) = real(c11, irealLUT)
            C(8*8)   = real(c12, irealLUT)

            handle_aspect_zx_1D_case = .True.
            ierr = OPP_1D_RETCODE

          elseif(aspect_zx.lt.OPP%dev%diffconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%dev%diffconfig%dims(3)%vrange(1)), &
              OPP%dev%diffconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        endif
      end function
  end subroutine

  subroutine get_coeff_cube(OPP, tauz, w0, g, aspect_zx, dir, C, ierr, angles, lswitch_east, lswitch_north, opt_vertices)
    class(t_optprop_cube)             :: OPP
    logical,intent(in)                :: dir
    real(irealLUT),intent(in)           :: tauz, w0, g, aspect_zx
    real(irealLUT),intent(in),optional  :: angles(:)
    logical,intent(in)                  :: lswitch_east, lswitch_north
    real(irealLUT),intent(out)          :: C(:)
    integer(mpiint), intent(out) :: ierr
    real(ireals), intent(in), optional :: opt_vertices(:)

    logical, save :: compute_coeff_online=.False., lset=.False.
    logical :: lflg
    real(ireals), allocatable :: vertices(:)
    real(irealLUT), allocatable :: Clut(:), Cbmc(:), Cbmc2(:)
    real(irealLUT) :: save_aspect_zx
    ierr = 0

    if(ldebug_optprop) call check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles)

    if(present(angles)) then ! obviously we want the direct coefficients
      if(aspect_zx.lt.OPP%dev%dirconfig%dims(3)%vrange(1)) then
        save_aspect_zx = OPP%dev%dirconfig%dims(3)%vrange(1)
        ierr = OPP_TINYASPECT_RETCODE
      else
        save_aspect_zx = aspect_zx
      endif
      if(dir) then ! dir2dir
        call OPP%dev%get_dir2dir([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
        call OPP%dir2dir_coeff_symmetry(C, lswitch_east, lswitch_north)
      else         ! dir2diff
        call OPP%dev%get_dir2diff([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
        call OPP%dir2diff_coeff_symmetry(C, lswitch_east, lswitch_north)
      endif
    else
      ! diff2diff
      if(aspect_zx.lt.OPP%dev%diffconfig%dims(3)%vrange(1)) then
        save_aspect_zx = OPP%dev%diffconfig%dims(3)%vrange(1)
        ierr = OPP_TINYASPECT_RETCODE
      else
        save_aspect_zx = aspect_zx
      endif
      call OPP%dev%get_diff2diff([tauz, w0, save_aspect_zx, g], C)
      call OPP%diff2diff_coeff_symmetry(C)
    endif

    if(.not.lset) then
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        & "-bmc_online", compute_coeff_online, lflg, ierr); call CHKERR(ierr)
      lset = .True.
    endif

    if(compute_coeff_online) then
      allocate(Clut(size(C)), Cbmc(size(C)))
      Clut = C
      call setup_default_unit_cube_geometry(one, one, real(aspect_zx, ireals), vertices)
      call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, Cbmc, angles)
      C = Cbmc

      if(present(opt_vertices)) then
        allocate(Cbmc2(size(C)))
        call get_coeff_bmc(OPP, opt_vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, Cbmc2, angles)
        C = Cbmc2
        print *,new_line(''),opt_vertices(3:24:3),':',angles,new_line('')//&
          cstr('LUT            '//toStr(Clut) , 'red')//new_line('')//&
          cstr('bmc (regular  )'//toStr(Cbmc) , 'blue' )//new_line('')//&
          cstr('bmc (distorted)'//toStr(Cbmc2), 'green')
      endif
    endif
  end subroutine

  subroutine get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: tauz, w0, g, vertices(:)
      logical,intent(in) :: dir
      real(irealLUT),intent(out):: C(:)
      real(irealLUT),intent(in),optional :: angles(2)

      real(irealLUT) :: S_diff(OPP%dev%diff_streams),T_dir(OPP%dev%dir_streams)
      real(irealLUT) :: S_tol (OPP%dev%diff_streams),T_tol(OPP%dev%dir_streams)
      integer(iintegers) :: isrc

      real(irealLUT), parameter :: atol=1e-2_irealLUT, rtol=5e-1_irealLUT

      if(present(angles)) then
          do isrc=1,OPP%dev%dir_streams
            call OPP%dev%bmc_wrapper(isrc, &
              real(vertices, ireal_dp), &
              real(tauz, irealLUT), &
              real(w0, irealLUT), &
              real(g, irealLUT), &
              .True.,   &
              real(angles(1), irealLUT), &
              real(angles(2), irealLUT), &
              mpi_comm_self, &
              S_diff, T_dir, S_tol, T_tol, &
              inp_atol=real(atol, irealLUT), inp_rtol=real(rtol, irealLUT))
            if(dir) then !dir2dir
              C(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams) = T_dir
            else ! dir2diff
              C(isrc:OPP%dev%dir_streams*OPP%dev%diff_streams:OPP%dev%dir_streams) = S_diff
            endif
            if(w0.ge.1) then
              if(any(S_tol.gt.0).or.any(T_tol.gt.0)) then
                print *,'SumT', sum(T_dir), 'SumS', sum(S_diff), &
                  'Divergence', (1- (sum(T_dir)+sum(S_diff))), any(S_tol.gt.0), any(T_tol.gt.0)
                if(abs(1- (sum(T_dir)+sum(S_diff))).ge.1e-6_irealLUT) then
                  call CHKWARN(1_mpiint, 'divergence '// &
                    toStr(1- (sum(T_dir)+sum(S_diff)))//' seems quite large for w0='//toStr(w0))
                endif
              endif
            endif
          enddo
      else
        ! diff2diff
        do isrc=1,OPP%dev%diff_streams
            call OPP%dev%bmc_wrapper(isrc, &
              real(vertices, ireal_dp), &
              real(tauz, irealLUT), &
              real(w0, irealLUT), &
              real(g, irealLUT), &
              .False.,   &
              0._irealLUT, 0._irealLUT, &
              mpi_comm_self, &
              S_diff, T_dir, S_tol, T_tol, &
              inp_atol=real(atol, irealLUT), inp_rtol=real(rtol, irealLUT))
          C(isrc:OPP%dev%diff_streams**2:OPP%dev%diff_streams) = S_diff
          if(w0.ge.1) then
            if(any(S_tol.gt.0).or.any(T_tol.gt.0)) then
              print *,'SumT', sum(T_dir), 'SumS', sum(S_diff), &
                'Divergence', (1- (sum(T_dir)+sum(S_diff))), any(S_tol.gt.0), any(T_tol.gt.0)
              if(abs(1- (sum(T_dir)+sum(S_diff))).ge.1e-6_irealLUT) then
                call CHKWARN(1_mpiint, 'divergence '// &
                  toStr(1- (sum(T_dir)+sum(S_diff)))//' seems quite large for w0='//toStr(w0))
              endif
            endif
          endif
        enddo
      endif ! angles_present

  end subroutine

  subroutine check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles, wedge_coords)
    class(t_optprop) :: OPP
    real(irealLUT),intent(in) :: tauz, w0, g, aspect_zx
    logical,intent(in) :: dir
    real(irealLUT),intent(in):: C(:)
    real(irealLUT),intent(in),optional  :: angles(:), wedge_coords(:)
    integer(iintegers), save :: dimidx_dir(8)=-2
    integer(iintegers), save :: dimidx_diff(6)=-2

    if( (any([tauz, w0, g, aspect_zx].lt.zero)) .or. (any(isnan([tauz, w0, g, aspect_zx]))) ) then
      call CHKERR(1_mpiint,'optprop_lookup_coeff :: '// &
        'corrupt optical properties: bg:: '//toStr([tauz, w0, g, aspect_zx]))
    endif

    if(dir) then
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,[      'tau'], tauz, dimidx_dir(1))
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,[       'w0'], w0  , dimidx_dir(2))
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,[        'g'], g   , dimidx_dir(3), default_val=0._irealLUT)
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,['aspect_zx'], aspect_zx, dimidx_dir(4))
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,['phi        ', 'param_phi  '], angles(1), dimidx_dir(5))
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,['theta      ', 'param_theta'], angles(2), dimidx_dir(6))
    else
      call check_LUT_dimension_limits(OPP%dev%diffconfig,[      'tau'], tauz, dimidx_diff(1))
      call check_LUT_dimension_limits(OPP%dev%diffconfig,[       'w0'], w0  , dimidx_diff(2))
      call check_LUT_dimension_limits(OPP%dev%diffconfig,[        'g'], g   , dimidx_diff(3), default_val=0._irealLUT)
      call check_LUT_dimension_limits(OPP%dev%diffconfig,['aspect_zx'], aspect_zx, dimidx_diff(4))
    endif

    if(present(angles)) then
      if(dir .and. size(C).ne. OPP%dev%dir_streams**2) then
        print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%dev%dir_streams**2
      endif
      if(.not.dir .and. size(C).ne. OPP%dev%diff_streams*OPP%dev%dir_streams) then
        print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C), &
          'should be',OPP%dev%diff_streams*OPP%dev%dir_streams
      endif
    else
      if(dir .and. size(C).ne. OPP%dev%diff_streams) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%dev%diff_streams
      endif
      if(.not.dir .and. size(C).ne. OPP%dev%diff_streams**2) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%dev%diff_streams**2
      endif
    endif

    if(present(wedge_coords)) then
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,['wedge_coord_Cx'], wedge_coords(1), dimidx_dir(7))
      call check_LUT_dimension_limits(OPP%dev%dirconfig ,['wedge_coord_Cy'], wedge_coords(2), dimidx_dir(8))
      call check_LUT_dimension_limits(OPP%dev%diffconfig ,['wedge_coord_Cx'], wedge_coords(1), dimidx_diff(5))
      call check_LUT_dimension_limits(OPP%dev%diffconfig ,['wedge_coord_Cy'], wedge_coords(2), dimidx_diff(6))
    endif
    contains
      subroutine check_LUT_dimension_limits(config, dimnames, val, dimindex, default_val)
        type(t_op_config), intent(in) :: config
        character(len=*), intent(in) :: dimnames(:)
        real(irealLUT), intent(in) :: val
        integer(iintegers), intent(inout) :: dimindex
        real(irealLUT), intent(in), optional :: default_val
        integer(iintegers) :: i, tmpidx
        if(dimindex.lt.-1) then
          tmpidx = dimindex
          do i = 1, size(dimnames)
            if(tmpidx.le.-1) dimindex = find_op_dim_by_name(config, trim(dimnames(i)))
            !print *,'Looking for dim: '//trim(dimnames(i)),' -> ', dimindex
          enddo
          dimindex = tmpidx
        endif
        if(dimindex.gt.-1) then
          if(.not.is_between(val, &
            config%dims(dimindex)%vrange(1)-tiny(val), &
            config%dims(dimindex)%vrange(2)+tiny(val))) then
              call CHKERR(1_mpiint, 'value ('//toStr(val)//') is not in the range of the LUT for dimension '// &
                trim(config%dims(dimindex)%dimname)//' ( '//toStr(config%dims(dimindex)%vrange)//' )')
          endif
        else if(dimindex.eq.-1) then
          if(.not.present(default_val)) then
            print *, 'Available LUT Dimensions are: '//char_arr_to_str(config%dims(:)%dimname, ', ')
            call CHKERR(1_mpiint, 'currently the LUT calls do not have a dimension for '// &
                char_arr_to_str(dimnames, ','))
          else
            if(.not.approx(val, default_val)) then
              print *, 'Available LUT Dimensions are: '//char_arr_to_str(config%dims(:)%dimname, ', ')
              call CHKERR(1_mpiint, 'currently the LUT calls do not have a dimension for '// &
                char_arr_to_str(dimnames, ',')// &
                ' and should probably be limited to the default val: '//toStr(default_val))
            endif
          endif
        endif
      end subroutine
  end subroutine

  subroutine diff2diff_coeff_symmetry(OPP, coeff)
    class(t_optprop_cube)        :: OPP
    real(irealLUT), target, intent(inout) :: coeff(:)
    real(irealLUT), pointer :: v(:,:) ! dim(src, dst)
    integer(iintegers) :: i
    real(irealLUT) :: norm(0:9)

    if(OPP%dev%diff_streams.eq.10) then
      v(0:9,0:9) => coeff(1:100)

      do i = 0,9
        norm(i) = sum(v(i,:))
      enddo

      v(0,0) = (v(0,0) + v(1,1)) * .5_irealLUT
      v(1,1) = v(0,0)

      v(0,1) = (v(0,1) + v(1,0)) * .5_irealLUT
      v(1,0) = v(0,1)

      v(1,0) = (v(0,1) + v(1,0)) * .5_irealLUT
      v(0,1) = v(1,0)


      ! v(0,2) = v(0,3) = v(0,6) = v(0,7) = v(1,4) = v(1,5) = v(1,8) = v(1,9) = np.mean(v(0,(2,3,6,7)) + v(1,(4,5,8,9)), axis=-1)/2
      v(0,2) = ( v(0,2) + v(0,3) + v(0,6) + v(0,7) + &
        & v(1,4) + v(1,5) + v(1,8) + v(1,9) ) / 8._irealLUT
      v(0,3) = v(0,2)
      v(0,6) = v(0,2)
      v(0,7) = v(0,2)
      v(1,4) = v(0,2)
      v(1,5) = v(0,2)
      v(1,8) = v(0,2)
      v(1,9) = v(0,2)

      !v(1,2) = v(1,3) = v(1,6) = v(1,7) = v(0,4) = v(0,5) = v(0,8) = v(0,9) = np.mean(v(1,(2,3,6,7)) + v(0,(4,5,8,9)), axis=-1)/2
      v(1,2) = ( v(0,4) + v(0,5) + v(0,8) + v(0,9) + &
        & v(1,2) + v(1,3) + v(1,6) + v(1,7) ) / 8._irealLUT
      v(1,3) = v(1,2)
      v(1,6) = v(1,2)
      v(1,7) = v(1,2)
      v(0,4) = v(1,2)
      v(0,5) = v(1,2)
      v(0,8) = v(1,2)
      v(0,9) = v(1,2)

      !v(2,0) = v(3,0) = v(4,1) = v(5,1) = v(6,0) = v(7,0) = v(8,1) = v(9,1) = np.mean(v((2,3,6,7),0) + v((4,5,8,9),1), axis=-1)/2
      v(2,0) = ( v(2,0) + v(3,0) + v(6,0) + v(7,0) + &
        & v(4,1) + v(5,1) + v(8,1) + v(9,1) ) / 8._irealLUT
      v(3,0) = v(2,0)
      v(6,0) = v(2,0)
      v(7,0) = v(2,0)
      v(4,1) = v(2,0)
      v(5,1) = v(2,0)
      v(8,1) = v(2,0)
      v(9,1) = v(2,0)


      !v(2,1) = v(3,1) = v(4,0) = v(5,0) = v(6,1) = v(7,1) = v(8,0) = v(9,0) = np.mean(v((2,3,6,7),1) + v((4,5,8,9),0), axis=-1)/2
      v(2,1) = ( v(4,0) + v(5,0) + v(8,0) + v(9,0) + &
        & v(2,1) + v(3,1) + v(6,1) + v(7,1) ) / 8._irealLUT
      v(4,0) = v(2,1)
      v(5,0) = v(2,1)
      v(8,0) = v(2,1)
      v(9,0) = v(2,1)
      v(3,1) = v(2,1)
      v(6,1) = v(2,1)
      v(7,1) = v(2,1)

      !v(2,2) = v(3,3) = v(4,4) = v(5,5) = v(6,6) = v(7,7) = v(8,8) = v(9,9) = np.mean(( v(_,_) for _ in range(2,10)), axis=0)
      v(2,2) = ( v(2,2) + v(3,3) + v(4,4) + v(5,5) + v(6,6) + v(7,7) + v(8,8) + v(9,9) ) / 8._irealLUT
      v(3,3) = v(2,2)
      v(4,4) = v(2,2)
      v(5,5) = v(2,2)
      v(6,6) = v(2,2)
      v(7,7) = v(2,2)
      v(8,8) = v(2,2)
      v(9,9) = v(2,2)

      !v(2,3) = v(3,2) = v(4,5) = v(5,4) = v(6,7) = v(7,6) = v(8,9) = v(9,8) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(3,2,5,4,7,6,9,8))), axis=0)
      v(2,3) = ( v(2,3) + v(3,2) + v(4,5) + v(5,4) + v(6,7) + v(7,6) + v(8,9) + v(9,8) ) / 8._irealLUT
      v(3,2) = v(2,3)
      v(4,5) = v(2,3)
      v(5,4) = v(2,3)
      v(6,7) = v(2,3)
      v(7,6) = v(2,3)
      v(8,9) = v(2,3)
      v(9,8) = v(2,3)

      !v(2,4) = v(3,5) = v(4,2) = v(5,3) = v(6,8) = v(7,9) = v(8,6) = v(9,7) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(4,5,2,3,8,9,6,7))), axis=0)
      v(2,4) = ( v(2,4) + v(3,5) + v(4,2) + v(5,3) + v(6,8) + v(7,9) + v(8,6) + v(9,7) ) / 8._irealLUT
      v(3,5) = v(2,4)
      v(4,2) = v(2,4)
      v(5,3) = v(2,4)
      v(6,8) = v(2,4)
      v(7,9) = v(2,4)
      v(8,6) = v(2,4)
      v(9,7) = v(2,4)

      !v(2,6) = v(2,7) = v(3,6) = v(3,7) = v(4,8) = v(4,9) = v(5,8) = v(5,9) = v(6,2) = v(6,3) = v(7,2) = v(7,3) = v(8,4) = v(8,5) = v(9,4) = v(9,5) = np.mean(( v(i1,i2) for i1,i2 in zip(sorted(list(range(2,10))*2),(6,7,6,7, 8,9,8,9, 2,3,2,3, 4,5,4,5,))), axis=0)
      v(2,6) = ( v(2,6) + v(2,7) + v(3,6) + v(3,7) + &
        & v(4,8) + v(4,9) + v(5,8) + v(5,9) + &
        & v(6,2) + v(6,3) + v(7,2) + v(7,3) + &
        & v(8,4) + v(8,5) + v(9,4) + v(9,5) ) / 16._irealLUT
      v(2,7) = v(2,6)
      v(3,6) = v(2,6)
      v(3,7) = v(2,6)
      v(4,8) = v(2,6)
      v(4,9) = v(2,6)
      v(5,8) = v(2,6)
      v(5,9) = v(2,6)
      v(6,2) = v(2,6)
      v(6,3) = v(2,6)
      v(7,2) = v(2,6)
      v(7,3) = v(2,6)
      v(8,4) = v(2,6)
      v(8,5) = v(2,6)
      v(9,4) = v(2,6)
      v(9,5) = v(2,6)

      !v(2,8) = v(2,9) = v(3,8) = v(3,9) = v(4,6) = v(4,7) = v(5,6) = v(5,7) = v(6,4) = v(6,5) = v(7,4) = v(7,5) = v(8,2) = v(8,3) = v(9,2) = v(9,3) = np.mean(( v(i1,i2) for i1,i2 in zip(sorted(list(range(2,10))*2),(8,9,8,9, 6,7,6,7, 4,5,4,5, 2,3,2,3,))), axis=0)
      v(2,8) = ( v(2,8) + v(2,9) + v(3,8) + v(3,9) + &
        & v(4,6) + v(4,7) + v(5,6) + v(5,7) + &
        & v(6,4) + v(6,5) + v(7,4) + v(7,5) + &
        & v(8,2) + v(8,3) + v(9,2) + v(9,3) ) / 16._irealLUT
      v(2,9) = v(2,8)
      v(3,8) = v(2,8)
      v(3,9) = v(2,8)
      v(4,6) = v(2,8)
      v(4,7) = v(2,8)
      v(5,6) = v(2,8)
      v(5,7) = v(2,8)
      v(6,4) = v(2,8)
      v(6,5) = v(2,8)
      v(7,4) = v(2,8)
      v(7,5) = v(2,8)
      v(8,2) = v(2,8)
      v(8,3) = v(2,8)
      v(9,2) = v(2,8)
      v(9,3) = v(2,8)

      !v(2,5) = v(3,4) = v(4,3) = v(5,2) = v(6,9) = v(7,8) = v(8,7) = v(9,6) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(5,4,3,2,9,8,7,6))), axis=0)
      v(2,5) = ( v(2,5) + v(3,4) + v(4,3) + v(5,2) + v(6,9) + v(7,8) + v(8,7) + v(9,6) ) / 8._irealLUT
      v(3,4) = v(2,5)
      v(4,3) = v(2,5)
      v(5,2) = v(2,5)
      v(6,9) = v(2,5)
      v(7,8) = v(2,5)
      v(8,7) = v(2,5)
      v(9,6) = v(2,5)

      do i = 0,9
        v(i,:) = v(i,:) / max(tiny(v), sum(v(i,:))) * norm(i)
      enddo
    endif

    if(.False.) then
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir2diff_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    return
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  !for solver_3_6 only the offset is changing in those sides which should be switched
  subroutine dir3_to_diff6_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_6)         :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 3
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      coeff(7:9)   = newcoeff([1, 2, 3] + dof*3)
      coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
    endif
    if(lswitch_north) then
      newcoeff = coeff
      coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
      coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  !for solver_3_10 the offset is chaning and the destination order
  subroutine dir3_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_10)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 3
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
       coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
       coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
       coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
       coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
    endif
    if(lswitch_north) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
       coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
       coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
       coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
       coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir3_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_16)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 3
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
       coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
       coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
       coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
       coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
       coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
       coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
       coeff(1+(11-1)*dof:11*dof) = newcoeff([1, 2, 3] + dof*(12-1) )
       coeff(1+(12-1)*dof:12*dof) = newcoeff([1, 2, 3] + dof*(11-1) )
      !coeff(1+(13-1)*dof:13*dof) = newcoeff([1, 2, 3] + dof*(13-1) )
      !coeff(1+(14-1)*dof:14*dof) = newcoeff([1, 2, 3] + dof*(14-1) )
      !coeff(1+(15-1)*dof:15*dof) = newcoeff([1, 2, 3] + dof*(15-1) )
      !coeff(1+(16-1)*dof:16*dof) = newcoeff([1, 2, 3] + dof*(16-1) )
    endif
    if(lswitch_north) then
      newcoeff = coeff
       coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
       coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
       coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
       coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
      !coeff(1+(11-1)*dof:11*dof) = newcoeff([1, 2, 3] + dof*(11-1) )
      !coeff(1+(12-1)*dof:12*dof) = newcoeff([1, 2, 3] + dof*(12-1) )
       coeff(1+(13-1)*dof:13*dof) = newcoeff([1, 2, 3] + dof*(14-1) )
       coeff(1+(14-1)*dof:14*dof) = newcoeff([1, 2, 3] + dof*(13-1) )
       coeff(1+(15-1)*dof:15*dof) = newcoeff([1, 2, 3] + dof*(16-1) )
       coeff(1+(16-1)*dof:16*dof) = newcoeff([1, 2, 3] + dof*(15-1) )
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  !for solver_8_10 the offset is chaning and the destination order
  subroutine dir8_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_10)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof = 8
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 2-1) )
       coeff(1+( 3-1)*dof: 3*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 4-1) )
       coeff(1+( 4-1)*dof: 4*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 3-1) )
       coeff(1+( 5-1)*dof: 5*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 6-1) )
       coeff(1+( 6-1)*dof: 6*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 5-1) )
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*(10-1) )
    endif
    if (lswitch_north) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 2-1) )
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 4-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 3-1) )
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 6-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 5-1) )
       coeff(1+( 7-1)*dof: 7*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 8-1) )
       coeff(1+( 8-1)*dof: 8*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 7-1) )
       coeff(1+( 9-1)*dof: 9*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*(10-1) )
       coeff(1+(10-1)*dof:10*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 9-1) )
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir8_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_16)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    real(irealLUT)               :: newcoeff(size(coeff)) ! dim(src,dst)
    integer(iintegers), parameter:: dof = 8
    if(lswitch_east) then
      newcoeff(:) = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 1-1)*dof)
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 2-1)*dof)
       coeff(1+( 3-1)*dof: 3*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 7-1)*dof)
       coeff(1+( 4-1)*dof: 4*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 8-1)*dof)
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 5-1)*dof)
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 6-1)*dof)
       coeff(1+( 7-1)*dof: 7*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 3-1)*dof)
       coeff(1+( 8-1)*dof: 8*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 4-1)*dof)
       coeff(1+( 9-1)*dof: 9*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(10-1)*dof)
       coeff(1+(10-1)*dof:10*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 9-1)*dof)
       coeff(1+(11-1)*dof:11*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(12-1)*dof)
       coeff(1+(12-1)*dof:12*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(11-1)*dof)
      !coeff(1+(13-1)*dof:13*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(13-1)*dof)
      !coeff(1+(14-1)*dof:14*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(14-1)*dof)
      !coeff(1+(15-1)*dof:15*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(15-1)*dof)
      !coeff(1+(16-1)*dof:16*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(16-1)*dof)
    endif
    if (lswitch_north) then
      newcoeff(:) = coeff
       coeff(1+( 1-1)*dof: 1*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 5-1)*dof)
       coeff(1+( 2-1)*dof: 2*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 6-1)*dof)
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 3-1)*dof)
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 4-1)*dof)
       coeff(1+( 5-1)*dof: 5*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 1-1)*dof)
       coeff(1+( 6-1)*dof: 6*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 2-1)*dof)
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 7-1)*dof)
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 8-1)*dof)
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 9-1)*dof)
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(10-1)*dof)
      !coeff(1+(11-1)*dof:11*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(11-1)*dof)
      !coeff(1+(12-1)*dof:12*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(12-1)*dof)
       coeff(1+(13-1)*dof:13*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(14-1)*dof)
       coeff(1+(14-1)*dof:14*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(13-1)*dof)
       coeff(1+(15-1)*dof:15*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(16-1)*dof)
       coeff(1+(16-1)*dof:16*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(15-1)*dof)
    endif
    if(lswitch_east) then
    endif
    if (lswitch_north) then
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir8_to_diff18_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_18)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'dir8_to_diff18_coeff_symmetry_lswitch_east_not yet implemented')
      coeff = coeff
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'dir8_to_diff18_coeff_symmetry_lswitch_north_not yet implemented')
      coeff = coeff
    endif
    select type(OPP)
    end select
  end subroutine

  subroutine dir2dir_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    return
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  subroutine dir2dir8_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube_dir8)   :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    integer(iintegers), parameter:: dof=8
    real(irealLUT)               :: newcoeff(size(coeff))

    if(lswitch_east) then
      newcoeff = coeff
      coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof  )
      coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]      )
      coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*3)
      coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*2)
      coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*4)
      coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*5)
      coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*6)
      coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*7)
    endif
    if (lswitch_north) then
      newcoeff = coeff
      coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*2)
      coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*3)
      coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]      )
      coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof  )
      coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*4)
      coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*5)
      coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*6)
      coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*7)
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
    endif
  end subroutine

  subroutine dir2dir3_coeff_corr_xx(coeffs, verts, verts_dtd, sundir)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:)
    real(ireals), intent(in) :: verts_dtd(:)
    real(ireals), intent(in) :: sundir(:)

    real(irealLUT) :: f, coeff_mod

    ! xx
    if (abs(sundir(1)) > epsilon(sundir(1)) * 10) then
      f = real(s(verts_dtd) / s(verts), irealLUT)
      coeff_mod = max(min((1._irealLUT - f) * coeffs(2), coeffs(2)), - coeffs(5))
      coeffs(2) = coeffs(2) - coeff_mod
      coeffs(5) = coeffs(5) + coeff_mod
    endif
    !! xy
    !if (abs(sundir(2)) > epsilon(sundir(2)) * 10) then
    !  s = max(min(abs(( &
    !    (verts(
    !endif

    contains
      function s(vertz)
        real(ireals), intent(in) :: vertz(:)
        real(ireals) :: s

        s = max(min(abs(( &
          (vertz(3) - vertz(6) + (vertz(4) - vertz(1)) * sundir(3) / sundir(1)) / (vertz(18) - vertz(6)) + &
          (vertz(9) - vertz(12) + (vertz(10) - vertz(7)) * sundir(3) / sundir(1)) / (vertz(24) - vertz(12)) &
          ) / 2._ireals), 1._ireals), 0._ireals)
      end function

  end subroutine

  subroutine dir2dir3_coeff_corr_xy(coeffs, verts, verts_dtd, sundir)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:)
    real(ireals), intent(in) :: verts_dtd(:)
    real(ireals), intent(in) :: sundir(:)

    real(irealLUT) :: f, coeff_mod

    if (abs(sundir(2)) > epsilon(sundir(2)) * 10) then
      f = real(s(verts) / s(verts_dtd), irealLUT)
      coeff_mod = max(min((1._irealLUT - f) * coeffs(2), coeffs(2)), - coeffs(8))
      print *, 'st', s(verts_dtd), 's', s(verts), 'f', f, 'coeff_mod', coeff_mod
      coeffs(2) = coeffs(2) + coeff_mod
      coeffs(8) = coeffs(8) - coeff_mod
    endif

    contains
      function s(vertz)
        real(ireals), intent(in) :: vertz(:)
        real(ireals) :: s

        s = max(min(abs(( &
          (vertz(3) - vertz(6) + (vertz(4) - vertz(1)) * sundir(3) / sundir(1)) / (vertz(18) - vertz(6)) + &
          (vertz(9) - vertz(12) + (vertz(10) - vertz(7)) * sundir(3) / sundir(1)) / (vertz(24) - vertz(12)) &
          ) / 2._ireals), 1._ireals), 0._ireals)
      end function
  end subroutine dir2dir3_coeff_corr_xy


  subroutine dir2dir3_coeff_corr_yy(coeffs, verts, verts_dtd, sundir)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:)
    real(ireals), intent(in) :: verts_dtd(:)
    real(ireals), intent(in) :: sundir(:)

    real(irealLUT) :: f, coeff_mod

    if (abs(sundir(2)) > epsilon(sundir(2)) * 10) then
      f = real(s(verts_dtd) / s(verts), irealLUT)
      coeff_mod = max(min((1._irealLUT - f) * coeffs(3), coeffs(3)), - coeffs(9))
      coeffs(3) = coeffs(3) - coeff_mod
      coeffs(9) = coeffs(9) + coeff_mod
    endif

    contains
      function s(vertz)
        real(ireals), intent(in) :: vertz(:)
        real(ireals) :: s

        s = max(min(abs(( &
          (vertz(9) - vertz(3) + (vertz(2) - vertz(8)) * sundir(3) / sundir(2)) / (vertz(15) - vertz(3)) + &
          (vertz(12) - vertz(6) + (vertz(5) - vertz(11)) * sundir(3) / sundir(2)) / (vertz(18) - vertz(6)) &
          ) / 2._ireals), 1._ireals), 0._ireals)
      end function
  end subroutine dir2dir3_coeff_corr_yy

  subroutine dir2dir3_coeff_corr_zx(coeffs, verts, verts_dtd, sundir)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:)
    real(ireals), intent(in) :: verts_dtd(:)
    real(ireals), intent(in) :: sundir(:)

    real(ireals) :: s, st
    real(irealLUT) :: f, coeff_mod

    if (abs(sundir(1)) > epsilon(sundir(1)) * 10) then
      s = max(min(abs(( &
        (verts(6) - verts(18)) / (verts(16) - verts(13)) * sundir(1) / sundir(3) + &
        (verts(12) - verts(24)) / (verts(22) - verts(19)) * sundir(1) / sundir(3) &
        ) / 2._ireals), 1._ireals), 0._ireals)
      st = max(min(abs(( &
        ((verts_dtd(18) - verts_dtd(6)) * sundir(1) / sundir(3)) /&
        (verts_dtd(13) - verts_dtd(16) + (verts_dtd(18) - verts_dtd(15)) * sundir(1) / sundir(3)) + &
        ((verts_dtd(24) - verts_dtd(12)) * sundir(1) / sundir(3)) /&
        (verts_dtd(19) - verts_dtd(22) + (verts_dtd(24) - verts_dtd(21)) * sundir(1) / sundir(3)) &
        ) / 2._ireals), 1._ireals), 0._ireals)
      f = real(abs(st) / abs(s), irealLUT)
      coeff_mod = max(min((1._irealLUT - f) * coeffs(4), coeffs(4)), - coeffs(1))
      coeffs(1) = coeffs(1) + coeff_mod
      coeffs(4) = coeffs(4) - coeff_mod
    endif

  end subroutine

  subroutine dir2dir3_coeff_corr_zy(coeffs, verts, verts_dtd, sundir)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:)
    real(ireals), intent(in) :: verts_dtd(:)
    real(ireals), intent(in) :: sundir(:)

    real(irealLUT) :: f, coeff_mod
    print *, 'called the functin', abs(sundir(2)), epsilon(sundir(2)) * 10
    if (abs(sundir(2)) > epsilon(sundir(2)) * 10) then
      f = real(s(verts_dtd) / s(verts), irealLUT)
      print *, 'st', s(verts_dtd), 's', s(verts), 'f', f
      coeff_mod = max(min((1._irealLUT - f) * coeffs(7), coeffs(7)), - coeffs(1))
      coeffs(1) = coeffs(1) + coeff_mod
      coeffs(7) = coeffs(7) - coeff_mod
    endif

    contains
      function s(vertz)
        real(ireals), intent(in) :: vertz(:)
        real(ireals) :: s

        s = max(min(abs(( &
          (vertz(24) - vertz(12)) / ((vertz(17) - vertz(23)) * sundir(3) / sundir(2) +  vertz(24) - vertz(18)) + &
          (vertz(21) - vertz(9)) / ((vertz(14) - vertz(20)) * sundir(3) / sundir(2) + vertz(21) - vertz(15)) &
          ) / 2._ireals), 1._ireals), 0._ireals)

      end function

  end subroutine dir2dir3_coeff_corr_zy

  subroutine dir2dir3_coeff_corr(verts, vertz, sundir, coeffs)
    real(irealLUT), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:), vertz(:), sundir(:)
    real(ireals) :: a(3), b(3), c(3), d(3), e(3), f(3), g(3), h(3)

    call reset_points()

    if (lDEBUG_geometric_coeff_correction) print *, 'sundir', sundir

    if(lDEBUG_geometric_coeff_correction) print *, 'src z'
    associate( &
        o1 => c, &
        o2 => d, &
        o3 => b, &
        o4 => a &
        )
      call project_points(sundir, h, compute_normal_3d(g, h, f), g, h, f, e, c, d, b, a)
      call rearange_projections(h, f, e, g, d, b, a, c)
      call correct_coeffs(h, f, e, g, d, b, a, c, o1, o2, o3, o4, [1,4,7], coeffs)
      call reset_points()
    end associate

    if (lDEBUG_geometric_coeff_correction)  print *, 'src x'
    associate( &
        o1 => h, &
        o2 => d, &
        o3 => b, &
        o4 => f &
        )
      call project_points(sundir, a, compute_normal_3d(c, a, e), g, c, a, e, h, d, b, f)
      call rearange_projections(g, c, a, e, h, d, b, f)
      call correct_coeffs(g, c, a, e, h, d, b, f, o1, o2, o3, o4, [5,8,2], coeffs)
      call reset_points()
    end associate

    if (lDEBUG_geometric_coeff_correction) print *, 'src y'
    associate( &
        o1 => g, &
        o2 => c, &
        o3 => d, &
        o4 => h &
        )
      call project_points(sundir, b, compute_normal_3d(a, b, f), e, a, b, f, g, c, d, h)
      call rearange_projections(f, b, a, e, h, d, c, g)
      call correct_coeffs(f, b, a, e, h, d, c, g, o1, o2, o3, o4, [9,6,3], coeffs)
      call reset_points()
    end associate

    if (.false.) then
      call respect_absorption(verts, vertz, sundir, coeffs)
    endif

    if (.false.) then
      call correct_by_gradient(e, f, g, h, coeffs)
    endif

  contains
    subroutine reset_points()
      a = verts(1:3)
      b = verts(4:6)
      c = verts(7:9)
      d = verts(10:12)
      e = verts(13:15)
      f = verts(16:18)
      g = verts(19:21)
      h = verts(22:24)
    end subroutine

    subroutine correct_by_gradient(e, f, g, h, coeffs)
      real(ireals), intent(in) :: e(3), f(3), g(3), h(3)
      real(irealLUT), intent(inout) :: coeffs(9)
      real(irealLUT) :: grad_x, grad_y

      print *, 'correct_by_gradient'
      ! src x
      grad_x = real(((f(3) - e(3)) / (f(1) - e(1)) + (h(3) - g(3)) / (h(1) - g(1))) / 2, irealLUT)
      if ( grad_x - real(sundir(1) / sundir(3), irealLUT) > epsilon(coeffs) ) then
        print *, 'correcting x'
        print *, 'beofre'
        print *, coeffs(2:9:3)
        !coeffs(5) = coeffs(5) + coeffs(5) / sum(coeffs([5,8])) * coeffs(2)
        !coeffs(8) = coeffs(8) + coeffs(8) / sum(coeffs([5,8])) * coeffs(2)
        !coeffs(2) = 0
      endif
      ! src y
      grad_y = real(((g(3) - e(3)) / (g(2) - e(2)) + (h(3) - f(3)) / (h(2) - f(2))) / 2, irealLUT)
      if ( grad_y - real(sundir(2) / sundir(3), irealLUT) > epsilon(coeffs) ) then
        print *, 'correcting y'
        print *, 'before'
        print *, coeffs(3:9:3)
        coeffs(3) = coeffs(3) + coeffs(3) / sum(coeffs([6,3])) * coeffs(9)
        coeffs(6) = coeffs(6) + coeffs(6) / sum(coeffs([6,3])) * coeffs(9)
        coeffs(9) = 0
      endif
    end subroutine

    subroutine project_points(sundir, origin, normal, f1, f2, f3, f4, v1, v2, v3, v4)
      real(ireals), intent(inout) :: f1(3), f2(3), f3(3), f4(3), v1(3), v2(3), v3(3), v4(3)
      real(ireals), intent(in) :: sundir(3), normal(3), origin(3)
      real(ireals), parameter :: eps = 1._ireals / sqrt(epsilon(eps)) ! try to delete this
      v1 = v1 + min(hit_plane(v1, sundir, origin, normal), eps) * sundir
      v2 = v2 + min(hit_plane(v2, sundir, origin, normal), eps) * sundir
      v3 = v3 + min(hit_plane(v3, sundir, origin, normal), eps) * sundir
      v4 = v4 + min(hit_plane(v4, sundir, origin, normal), eps) * sundir

      v1 = vec_proj_on_plane(v1, normal)
      v2 = vec_proj_on_plane(v2, normal)
      v3 = vec_proj_on_plane(v3, normal)
      v4 = vec_proj_on_plane(v4, normal)

      f1 = vec_proj_on_plane(f1, normal)
      f2 = vec_proj_on_plane(f2, normal)
      f3 = vec_proj_on_plane(f3, normal)
      f4 = vec_proj_on_plane(f4, normal)
    end subroutine

    subroutine rearange_projections(f1, f2, f3, f4, v1, v2, v3, v4)
      real(ireals), intent(in) :: f1(3), f2(3), f3(3), f4(3)
      real(ireals), intent(inout) :: v1(3), v2(3), v3(3), v4(3)

      if (lDEBUG_geometric_coeff_correction) then
        print *, 'before'
        print *, 'f1', f1
        print *, 'f2', f2
        print *, 'f3', f3
        print *, 'f4', f4
        print *, '______________________'
        print *, 'v1', v1
        print *, 'v2', v2
        print *, 'v3', v3
        print *, 'v4', v4
        print *, 'after'
      endif

      call rearange_projection(f1-v1, f3, f4-f3, f2, f3-f2, v1)
      call rearange_projection(f2-v2, f3, f4-f3, f4, f1-f4, v2)
      call rearange_projection(f3-v3, f2, f1-f2, f4, f1-f4, v3)
      call rearange_projection(f4-v4, f2, f1-f2, f2, f3-f2, v4)

      if (lDEBUG_geometric_coeff_correction) then
        print *, 'v1', v1
        print *, 'v2', v2
        print *, 'v3', v3
        print *, 'v4', v4
      endif
    end subroutine

    subroutine rearange_projection(direction1, origin2, direction2, origin3, direction3, origin1)
      real(ireals), intent(in) :: direction1(3), origin2(3), direction2(3), origin3(3), direction3(3)
      real(ireals), intent(inout) :: origin1(3)
      real(ireals) :: coeff21, coeff22, coeff31, coeff32
      integer(mpiint) :: ierr

      call line_intersection_3d(origin1, direction1, origin2, direction2, coeff21, coeff22, ierr)
      call line_intersection_3d(origin1, direction1, origin3, direction3, coeff31, coeff32, ierr)
      call rearange_point(origin1, direction1, min(max(coeff21, coeff31, zero), one), origin1)
    end subroutine

    subroutine rearange_point(origin, direction, coefficient, point)
      real(ireals), intent(inout) :: point(3)
      real(ireals), intent(in) :: origin(3), direction(3), coefficient

      point = origin + coefficient * direction
    end subroutine

    subroutine proj_var_to_edges(c1, c2, c3, c4, v, p1, p2, p3) ! c: corner
      real(ireals), intent(in) :: c1(3), c2(3), c3(3), c4(3), v(3)
      real(ireals), intent(out) :: p1(3), p2(3), p3(3)
      real(ireals) :: c, t
      integer(mpiint) :: ierr

      call line_intersection_3d(v, c4-c1, c3, c4-c3, c, t, ierr)
      call rearange_point(v, c4-c1, c, p1)
      call line_intersection_3d(v, c2-c1, c2, c3-c2, c, t, ierr)
      call rearange_point(v, c2-c1, c, p2)
      call line_intersection_3d(v, c1-c2, c4, c1-c4, c, t, ierr)
      call rearange_point(v, c1-c2, c, p3)
    end subroutine

    subroutine proj_vars_to_edges(f1,f2,f3,f4,v1,v2,v3,v4,p1l,p1b,p1t,p2l,p2t,p2b,p3r,p3t,p3b,p4r,p4b,p4t)
      real(ireals), intent(in) :: f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),v4(3)
      real(ireals), intent(out) :: p1l(3),p1b(3),p1t(3),p2l(3),p2t(3),p2b(3),p3r(3),p3t(3),p3b(3),p4r(3),p4b(3),p4t(3)

      call proj_var_to_edges(f1, f2, f3, f4, v1, p1l, p1b, p1t)
      call proj_var_to_edges(f2, f1, f4, f3, v2, p2l, p2t, p2b)
      call proj_var_to_edges(f4, f3, f2, f1, v3, p3r, p3b, p3t)
      call proj_var_to_edges(f4, f3, f2, f1, v4, p4r, p4b, p4t)
    end subroutine

    subroutine correct_coeffs(f1, f2, f3, f4, v1, v2, v3, v4, o1, o2, o3, o4, slice, coeffs)
      real(ireals), intent(in) :: f1(3), f2(3), f3(3), f4(3), v1(3), v2(3), v3(3), v4(3), o1(3), o2(3), o3(3), o4(3)
      integer(iintegers), intent(in) :: slice(3)
      real(irealLUT), intent(inout) :: coeffs(9)
      real(ireals) :: area1, area2, area3, area, areas(3), a2v1, a2v2, a2v3, a2v4, a3v1, a3v2, a3v3, a3v4
      real(ireals) :: v2v1, v2v2, v2v3, v2v4, v3v1, v3v2, v3v3, v3v4, volume2, volume3
      real(ireals) :: p1l(3), p1b(3), p2l(3), p2t(3), p3r(3), p3t(3), p4r(3), p4b(3), p1t(3), p2b(3), p3b(3), p4t(3)

      call proj_vars_to_edges(f1, f2, f3, f4, v1, v2, v3, v4, p1l, p1b, p1t, p2l, p2t, p2b, p3r, p3t, p3b, p4r, p4b, p4t)

      if (lDEBUG_geometric_coeff_correction) then
        print *, 'p1b', p1b, 'p1l', p1l, 'p1t', p1t
        print *, 'p2t', p2t, 'p2l', p2l, 'p2b', p2b
        print *, 'p3t', p3t, 'p3r', p3r, 'p3b', p3b
        print *, 'p4b', p4b, 'p4r', p4r, 'p4t', p4t
      endif

      area = triangle_area_by_vertices(f1, f2, f3) + triangle_area_by_vertices(f1, f3, f4)

      a2v1 = triangle_area_by_vertices(v1,f2, p1b) + triangle_area_by_vertices(v1, f1,f2)
      a2v2 = triangle_area_by_vertices(v2,p2t,f1 ) + triangle_area_by_vertices(p2t,f1,f2)
      a2v3 = triangle_area_by_vertices(v3,f4, p3t) + triangle_area_by_vertices(v3, f3,f4)
      a2v4 = triangle_area_by_vertices(v4,p4b,f3 ) + triangle_area_by_vertices(v4, f3,f4)

      area2 = max(a2v1, a2v2, a2v3, a2v4)

      a3v1 = triangle_area_by_vertices(v1,p1t,f1) + triangle_area_by_vertices(v1,p1l,f4 ) + triangle_area_by_vertices(v1,f4,p1t)
      a3v2 = triangle_area_by_vertices(v2,p2b,f3) + triangle_area_by_vertices(v2,f3, p2l) + triangle_area_by_vertices(v2,f2,p2b)
      a3v3 = triangle_area_by_vertices(v3,p3b,f3) + triangle_area_by_vertices(v3,p3b,f2 ) + triangle_area_by_vertices(v3,f2,p3b)
      a3v4 = triangle_area_by_vertices(v4,p4t,f1) + triangle_area_by_vertices(v4,f4, p4t) + triangle_area_by_vertices(v4,f1,p4r)

      area3 = max(a3v1, a3v2, a3v3, a3v4)

      area1 = area - area2 - area3

      areas = max([area1,area2,area3], zero)
      areas = areas / sum(areas)

      coeffs(slice) = real(areas, irealLUT) * sum(coeffs(slice))

      ! absorption correction
      v2v1 = volume_hexahedron(o1, o2, o3, o4, v1, f1, f2, p1b)
      v2v2 = volume_hexahedron(o1, o2, o3, o4, v2, p2t, f1, f2)
      v2v3 = volume_hexahedron(o1, o2, o3, o4, v3, f3, f4, p3t)
      v2v4 = volume_hexahedron(o1, o2, o3, o4, v4, f4, f1, p4r)

      volume2 = max(v2v1, v2v2, v2v3, v2v4)

      v3v1 = volume_hexahedron(o1, o2, o3, o4, v1, p1l, f4, f1)
      v3v2 = volume_hexahedron(o1, o2, o3, o4, v2, f2, f3, p2l)
      v3v3 = volume_hexahedron(o1, o2, o3, o4, v3, f3, f4, p3t)
      v3v4 = volume_hexahedron(o1, o2, o3, o4, v4, f4, f1, p4r)

      volume3 = max(v3v1, v3v2, v3v3, v3v4)
      ! absorption correction example
      !coeffs(1) = real(exp( - 0.5_ireals * (d_prime_dst - d_prime_reg)), irealLUT) * coeffs(1)

      ! THIS IS HOW IT WORKS
      !coeffs_total = sum(coeffs(3:9:3))
      !coeffs([9,6,3]) = real(areas, irealLUT) * coeffs_total

      if (lDEBUG_geometric_coeff_correction) print *, 'areas', areas
    end subroutine
  end subroutine dir2dir3_coeff_corr
end module
