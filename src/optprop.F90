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

#ifdef _XLF
      use ieee_arithmetic
#define isnan ieee_is_nan
#endif

use m_optprop_parameters, only : ldebug_optprop, coeff_mode, wedge_sphere_radius, param_eps
use m_helper_functions, only : rmse, CHKERR, CHKWARN, itoa, ftoa, approx, deg2rad, rad2deg, swap, is_between, char_arr_to_str
use m_data_parameters, only: ireals,ireal_dp,irealLUT,ireal_params,iintegers,one,zero,i0,i1,inil,mpiint
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_3_6, t_optprop_LUT_3_10, &
  t_optprop_LUT_8_10, t_optprop_LUT_3_16, t_optprop_LUT_8_16, t_optprop_LUT_8_18, &
  t_optprop_LUT_wedge_5_8, t_optprop_LUT_rectilinear_wedge_5_8, t_optprop_LUT_wedge_18_8, &
  find_lut_dim_by_name, t_lut_config
use m_optprop_ANN, only : ANN_init, ANN_get_dir2dir, ANN_get_dir2diff, ANN_get_diff2diff
use m_boxmc_geometry, only : setup_default_unit_cube_geometry, setup_default_wedge_geometry
use m_eddington, only: eddington_coeff_zdun
use m_tenstream_options, only: twostr_ratio

use m_LUT_param_phi, only: theta_from_param_theta, iterative_phi_theta_from_param_phi_and_param_theta

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public :: t_optprop, t_optprop_cube, t_optprop_wedge, &
  t_optprop_1_2, t_optprop_3_6, t_optprop_3_10, &
  t_optprop_wedge_5_8, t_optprop_rectilinear_wedge_5_8, t_optprop_wedge_18_8, &
  t_optprop_8_10, t_optprop_3_16, t_optprop_8_16, t_optprop_8_18, &
  OPP_1D_RETCODE, OPP_TINYASPECT_RETCODE

type,abstract :: t_optprop
  logical :: optprop_debug=ldebug_optprop
  class(t_optprop_LUT), allocatable :: OPP_LUT
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

contains

  subroutine init(OPP, comm, skip_load_LUT)
      class(t_optprop), intent(inout) :: OPP
      integer(mpiint) ,intent(in) :: comm
      logical, intent(in), optional :: skip_load_LUT
      integer(mpiint) :: ierr

      select case (coeff_mode)
          case(i0) ! LookUpTable Mode
            select type(OPP)
              class is (t_optprop_1_2)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_1_2::OPP%OPP_LUT)

              class is (t_optprop_3_6)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_6::OPP%OPP_LUT)

              class is (t_optprop_3_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_10::OPP%OPP_LUT)

              class is (t_optprop_8_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_10::OPP%OPP_LUT)

              class is (t_optprop_3_16)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_16::OPP%OPP_LUT)

              class is (t_optprop_8_16)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_16::OPP%OPP_LUT)

              class is (t_optprop_8_18)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_18::OPP%OPP_LUT)

              class is (t_optprop_wedge_5_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_wedge_5_8::OPP%OPP_LUT)

              class is (t_optprop_rectilinear_wedge_5_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_rectilinear_wedge_5_8::OPP%OPP_LUT)

              class is (t_optprop_wedge_18_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_wedge_18_8::OPP%OPP_LUT)
              class default
                call CHKERR(1_mpiint, ' init optprop : unexpected type for optprop object!')
            end select
            call OPP%OPP_LUT%init(comm, skip_load_LUT)

          case(i1) ! ANN
            call ANN_init(comm, ierr)
  !          stop 'ANN not yet implemented'
          case default
            call CHKERR(1_mpiint, 'coeff mode optprop initialization not defined')
        end select
  end subroutine
  subroutine destroy(OPP)
      class(t_optprop) :: OPP
      if(allocated(OPP%OPP_LUT)) then
          call OPP%OPP_LUT%destroy()
          deallocate(OPP%OPP_LUT)
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
        do isrc=1,OPP%OPP_LUT%dir_streams
          print *, 'lut src', isrc, ':', C(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams)
          print *, 'bmc src', isrc, ':', Cbmc(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams)
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
          dimidx_dir(1) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'tau')
          dimidx_dir(2) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'w0')
          dimidx_dir(3) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'aspect_zx')
          dimidx_dir(4) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'g')
          dimidx_dir(5) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cx')
          dimidx_dir(6) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cy')
          dimidx_dir(7) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'param_phi')
          dimidx_dir(8) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'param_theta')
          allocate(inp_arr_dir(count(dimidx_dir.gt.0)))

          dimidx_diff(1) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'tau')
          dimidx_diff(2) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'w0')
          dimidx_diff(3) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'aspect_zx')
          dimidx_diff(4) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'g')
          dimidx_diff(5) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cx')
          dimidx_diff(6) = find_lut_dim_by_name(OPP%OPP_LUT%dirconfig, 'wedge_coord_Cy')
          allocate(inp_arr_diff(count(dimidx_diff.gt.0)))
          linit = .True.
        endif

          select case (coeff_mode)
          case(i0) ! LookUpTable Mode
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
                  call OPP%OPP_LUT%LUT_get_dir2dir(inp_arr_dir, C)
                else ! dir2diff
                  call OPP%OPP_LUT%LUT_get_dir2diff(inp_arr_dir, C)
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
              call OPP%OPP_LUT%LUT_get_diff2diff(inp_arr_diff, C)
            endif

          case default
            call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
          end select
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

          elseif(aspect_zx.lt.OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)), &
              OPP%OPP_LUT%dirconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        else ! diffuse

          !if(aspect_zx.gt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(2)) then
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

          elseif(aspect_zx.lt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)) then
            restricted_aspect_zx = min(max(aspect_zx, OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)), &
              OPP%OPP_LUT%diffconfig%dims(3)%vrange(2))
            call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
            handle_aspect_zx_1D_case = .True.
            ierr = OPP_TINYASPECT_RETCODE
          endif

        endif
      end function
  end subroutine

  subroutine get_coeff_cube(OPP, tauz, w0, g, aspect_zx, dir, C, ierr, angles, lswitch_east, lswitch_north)
    class(t_optprop_cube)             :: OPP
    logical,intent(in)                :: dir
    real(irealLUT),intent(in)           :: tauz, w0, g, aspect_zx
    real(irealLUT),intent(in),optional  :: angles(:)
    logical,intent(in)                  :: lswitch_east, lswitch_north
    real(irealLUT),intent(out)          :: C(:)
    integer(mpiint), intent(out) :: ierr

    logical,parameter :: compute_coeff_online=.False.
    real(ireals), allocatable :: vertices(:)
    real(irealLUT) :: save_aspect_zx
    ierr = 0

    if(compute_coeff_online) then
      call setup_default_unit_cube_geometry(one, one, real(aspect_zx, ireals), vertices)
      call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, C, angles)
      return
    endif

    if(ldebug_optprop) call check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles)


    select case (coeff_mode)

    case(i0) ! LookUpTable Mode

      if(present(angles)) then ! obviously we want the direct coefficients
        if(aspect_zx.lt.OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)) then
          save_aspect_zx = OPP%OPP_LUT%dirconfig%dims(3)%vrange(1)
          ierr = OPP_TINYASPECT_RETCODE
        else
          save_aspect_zx = aspect_zx
        endif
        if(dir) then ! dir2dir
          call OPP%OPP_LUT%LUT_get_dir2dir([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
          call OPP%dir2dir_coeff_symmetry(C, lswitch_east, lswitch_north)
        else         ! dir2diff
          call OPP%OPP_LUT%LUT_get_dir2diff([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
          call OPP%dir2diff_coeff_symmetry(C, lswitch_east, lswitch_north)
        endif
      else
        ! diff2diff
        if(aspect_zx.lt.OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)) then
          save_aspect_zx = OPP%OPP_LUT%diffconfig%dims(3)%vrange(1)
          ierr = OPP_TINYASPECT_RETCODE
        else
          save_aspect_zx = aspect_zx
        endif
        call OPP%OPP_LUT%LUT_get_diff2diff([tauz, w0, save_aspect_zx, g], C)
      endif


    case(i1) ! ANN

      if(present(angles)) then ! obviously we want the direct coefficients
        if(dir) then ! specifically the dir2dir
          call ANN_get_dir2dir(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
        else ! dir2diff
          call ANN_get_dir2diff(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
        endif
      else
        ! diff2diff
        call ANN_get_diff2diff(tauz, w0, g, aspect_zx, C)
      endif

    case default
      call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
    end select

  end subroutine

  subroutine get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: tauz, w0, g, vertices(:)
      logical,intent(in) :: dir
      real(irealLUT),intent(out):: C(:)
      real(irealLUT),intent(in),optional :: angles(2)

      real(irealLUT) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
      real(irealLUT) :: S_tol (OPP%OPP_LUT%diff_streams),T_tol(OPP%OPP_LUT%dir_streams)
      integer(iintegers) :: isrc

      real(irealLUT), parameter :: atol=5e-3_irealLUT, rtol=5e-2_irealLUT

      if(present(angles)) then
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, &
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
              C(isrc:OPP%OPP_LUT%dir_streams**2:OPP%OPP_LUT%dir_streams) = T_dir
            else ! dir2diff
              C(isrc:OPP%OPP_LUT%dir_streams*OPP%OPP_LUT%diff_streams:OPP%OPP_LUT%dir_streams) = S_diff
            endif
            if(w0.ge.1) then
              if(any(S_tol.gt.0).or.any(T_tol.gt.0)) then
                print *,'SumT', sum(T_dir), 'SumS', sum(S_diff), &
                  'Divergence', (1- (sum(T_dir)+sum(S_diff))), any(S_tol.gt.0), any(T_tol.gt.0)
                if(abs(1- (sum(T_dir)+sum(S_diff))).ge.1e-6_irealLUT) then
                  call CHKWARN(1_mpiint, 'divergence '// &
                    ftoa(1- (sum(T_dir)+sum(S_diff)))//' seems quite large for w0='//ftoa(w0))
                endif
              endif
            endif
          enddo
      else
        ! diff2diff
        do isrc=1,OPP%OPP_LUT%diff_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, &
              real(vertices, ireal_dp), &
              real(tauz, irealLUT), &
              real(w0, irealLUT), &
              real(g, irealLUT), &
              .False.,   &
              0._irealLUT, 0._irealLUT, &
              mpi_comm_self, &
              S_diff, T_dir, S_tol, T_tol, &
              inp_atol=real(atol, irealLUT), inp_rtol=real(rtol, irealLUT))
          C(isrc:OPP%OPP_LUT%diff_streams**2:OPP%OPP_LUT%diff_streams) = S_diff
          if(w0.ge.1) then
            if(any(S_tol.gt.0).or.any(T_tol.gt.0)) then
              print *,'SumT', sum(T_dir), 'SumS', sum(S_diff), &
                'Divergence', (1- (sum(T_dir)+sum(S_diff))), any(S_tol.gt.0), any(T_tol.gt.0)
              if(abs(1- (sum(T_dir)+sum(S_diff))).ge.1e-6_irealLUT) then
                call CHKWARN(1_mpiint, 'divergence '// &
                  ftoa(1- (sum(T_dir)+sum(S_diff)))//' seems quite large for w0='//ftoa(w0))
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
        'corrupt optical properties: bg:: '//ftoa([tauz, w0, g, aspect_zx]))
    endif

    if(dir) then
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,[      'tau'], tauz, dimidx_dir(1))
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,[       'w0'], w0  , dimidx_dir(2))
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,[        'g'], g   , dimidx_dir(3), default_val=0._irealLUT)
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,['aspect_zx'], aspect_zx, dimidx_dir(4))
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,['phi        ', 'param_phi  '], angles(1), dimidx_dir(5))
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,['theta      ', 'param_theta'], angles(2), dimidx_dir(6))
    else
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig,[      'tau'], tauz, dimidx_diff(1))
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig,[       'w0'], w0  , dimidx_diff(2))
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig,[        'g'], g   , dimidx_diff(3), default_val=0._irealLUT)
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig,['aspect_zx'], aspect_zx, dimidx_diff(4))
    endif

    if(present(angles)) then
      if(dir .and. size(C).ne. OPP%OPP_LUT%dir_streams**2) then
        print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%dir_streams**2
      endif
      if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams) then
        print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C), &
          'should be',OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams
      endif
    else
      if(dir .and. size(C).ne. OPP%OPP_LUT%diff_streams) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams
      endif
      if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams**2) then
        print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams**2
      endif
    endif

    if(present(wedge_coords)) then
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,['wedge_coord_Cx'], wedge_coords(1), dimidx_dir(7))
      call check_LUT_dimension_limits(OPP%OPP_LUT%dirconfig ,['wedge_coord_Cy'], wedge_coords(2), dimidx_dir(8))
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig ,['wedge_coord_Cx'], wedge_coords(1), dimidx_diff(5))
      call check_LUT_dimension_limits(OPP%OPP_LUT%diffconfig ,['wedge_coord_Cy'], wedge_coords(2), dimidx_diff(6))
    endif
    contains
      subroutine check_LUT_dimension_limits(config, dimnames, val, dimindex, default_val)
        type(t_lut_config), intent(in) :: config
        character(len=*), intent(in) :: dimnames(:)
        real(irealLUT), intent(in) :: val
        integer(iintegers), intent(inout) :: dimindex
        real(irealLUT), intent(in), optional :: default_val
        integer(iintegers) :: i, tmpidx
        if(dimindex.lt.-1) then
          tmpidx = dimindex
          do i = 1, size(dimnames)
            if(tmpidx.le.-1) dimindex = find_lut_dim_by_name(config, trim(dimnames(i)))
            !print *,'Looking for dim: '//trim(dimnames(i)),' -> ', dimindex
          enddo
          dimindex = tmpidx
        endif
        if(dimindex.gt.-1) then
          if(.not.is_between(val, &
            config%dims(dimindex)%vrange(1)-tiny(val), &
            config%dims(dimindex)%vrange(2)+tiny(val))) then
              call CHKERR(1_mpiint, 'value ('//ftoa(val)//') is not in the range of the LUT for dimension '// &
                trim(config%dims(dimindex)%dimname)//' ( '//ftoa(config%dims(dimindex)%vrange)//' )')
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
                ' and should probably be limited to the default val: '//ftoa(default_val))
            endif
          endif
        endif
      end subroutine
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
      !coeff( 1: 3) = newcoeff([1, 2, 3]        )
      !coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
      coeff( 7: 9) = newcoeff([1, 2, 3] + dof*3)
      coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
      coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
      coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
      ! coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
      ! coeff(22:24) = newcoeff([1, 2, 3] + dof*7)
      ! coeff(25:27) = newcoeff([1, 2, 3] + dof*8)
      ! coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
    endif
    if(lswitch_north) then
      newcoeff = coeff
      !coeff( 1: 3) = newcoeff([1, 2, 3]        )
      !coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
      !coeff( 7: 9) = newcoeff([1, 2, 3] + dof*2)
      !coeff(10:12) = newcoeff([1, 2, 3] + dof*3)
      !coeff(13:15) = newcoeff([1, 2, 3] + dof*4)
      !coeff(16:18) = newcoeff([1, 2, 3] + dof*5)
      coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
      coeff(22:24) = newcoeff([1, 2, 3] + dof*8)
      coeff(25:27) = newcoeff([1, 2, 3] + dof*7)
      coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
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
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
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
    integer(iintegers), parameter :: dof = 8
    real(irealLUT)               :: newcoeff(size(coeff))
    if(lswitch_east) then
      newcoeff = coeff
      !coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]        )
      !coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*1 )
      coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*3 )
      coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*2 )
      coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*5 )
      coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*4 )
      !coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*6 )
      !coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*7 )
      !coeff(65:72) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*8 )
      !coeff(73:80) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*9 )
    endif
    if (lswitch_north) then
      newcoeff = coeff
      !coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]        )
      !coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*1 )
      !coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*2 )
      !coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*3 )
      !coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*4 )
      !coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*5 )
      coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*6 )
      coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*8 )
      coeff(65:72) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*7 )
      coeff(73:80) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*9 )
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
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
    endif
    if(.False.) then ! remove compiler unused warnings
      select type(OPP)
      end select
      if(lswitch_east .or. lswitch_north) coeff=coeff
    endif
  end subroutine

  subroutine dir8_to_diff18_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_18)        :: OPP
    logical, intent(in)          :: lswitch_east, lswitch_north
    real(irealLUT),intent(inout) :: coeff(:)
    if(lswitch_east) then
      call CHKERR(1_mpiint, 'not yet implemented')
      coeff = coeff
    endif
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'not yet implemented')
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
end module
