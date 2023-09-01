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

  use m_optprop_parameters, only: ldebug_optprop, wedge_sphere_radius, param_eps
  use m_helper_functions, only: &
    & approx, &
    & char_arr_to_str, &
    & CHKERR, CHKWARN, &
    & cstr, &
    & deg2rad, &
    & get_petsc_opt, &
    & is_inrange, &
    & rad2deg, &
    & rmse, &
    & swap, &
    & toStr
  use m_data_parameters, only: ireals, ireal_dp, irealLUT, ireal_params, iintegers, one, zero, i0, i1, inil, mpiint
  use m_optprop_base, only: t_optprop_base, t_op_config, find_op_dim_by_name
  use m_optprop_LUT, only: &
    & t_optprop_LUT, &
    & t_optprop_LUT_1_2, &
    & t_optprop_LUT_3_6, &
    & t_optprop_LUT_3_10, &
    & t_optprop_LUT_3_16, &
    & t_optprop_LUT_3_24, &
    & t_optprop_LUT_3_30, &
    & t_optprop_LUT_8_10, &
    & t_optprop_LUT_8_16, &
    & t_optprop_LUT_8_18, &
    & t_optprop_LUT_wedge_5_8, &
    & t_optprop_LUT_rectilinear_wedge_5_8, &
    & t_optprop_LUT_wedge_18_8

  use m_optprop_ANN, only: t_optprop_ANN, t_optprop_ANN_3_10
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry, setup_default_wedge_geometry
  use m_eddington, only: eddington_coeff_ec
  use m_tenstream_options, only: twostr_ratio

  use m_LUT_param_phi, only: theta_from_param_theta, iterative_phi_theta_from_param_phi_and_param_theta

  use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

  implicit none

  private
  public :: &
    t_optprop, &
    t_optprop_cube, &
    t_optprop_wedge, &
    t_optprop_1_2, &
    t_optprop_3_6, &
    t_optprop_3_10, &
    t_optprop_3_10_ann, &
    t_optprop_8_10, &
    t_optprop_3_16, &
    t_optprop_3_24, &
    t_optprop_3_30, &
    t_optprop_8_16, &
    t_optprop_8_18, &
    t_optprop_wedge_5_8, &
    t_optprop_rectilinear_wedge_5_8, &
    t_optprop_wedge_18_8, &
    OPP_1D_RETCODE, &
    OPP_TINYASPECT_RETCODE

  type, abstract :: t_optprop
    logical :: optprop_debug = ldebug_optprop
    class(t_optprop_LUT), allocatable :: LUT
    class(t_optprop_ANN), allocatable :: ANN
    class(t_optprop_base), pointer :: dev
  contains
    procedure :: init
    procedure :: get_coeff_bmc
    procedure :: destroy
  end type

! Cube types
  type, abstract, extends(t_optprop) :: t_optprop_cube
  contains
    procedure :: get_coeff => get_coeff_cube
    procedure :: dir2dir_coeff_symmetry => dir2dir_coeff_symmetry_none
    procedure :: dir2diff_coeff_symmetry => dir2diff_coeff_symmetry_none
    procedure :: diff2diff_coeff_symmetry => diff2diff_coeff_symmetry
  end type

! we introduce one special cube type for 8 direct streams, this way, all of them can share dir2dir_coeff_symmetry
  type, abstract, extends(t_optprop_cube) :: t_optprop_cube_dir8
  contains
    procedure :: dir2dir_coeff_symmetry => dir2dir8_coeff_symmetry
  end type

  type, extends(t_optprop_cube) :: t_optprop_1_2
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_6
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff6_coeff_symmetry
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff10_coeff_symmetry
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_10_ann
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff16_coeff_symmetry
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_24
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff24_coeff_symmetry
  end type

  type, extends(t_optprop_cube) :: t_optprop_3_30
  contains
    procedure :: dir2diff_coeff_symmetry => dir3_to_diff30_coeff_symmetry
  end type

  type, extends(t_optprop_cube_dir8) :: t_optprop_8_10
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff10_coeff_symmetry
  end type

  type, extends(t_optprop_cube_dir8) :: t_optprop_8_16
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff16_coeff_symmetry
  end type

  type, extends(t_optprop_cube_dir8) :: t_optprop_8_18
  contains
    procedure :: dir2diff_coeff_symmetry => dir8_to_diff18_coeff_symmetry
  end type

! Wedge types
  type, abstract, extends(t_optprop) :: t_optprop_wedge
  contains
    procedure :: get_coeff => get_coeff_wedge
  end type

  type, extends(t_optprop_wedge) :: t_optprop_wedge_5_8
  end type

  type, extends(t_optprop_wedge) :: t_optprop_rectilinear_wedge_5_8
  end type

  type, extends(t_optprop_wedge) :: t_optprop_wedge_18_8
  end type

  integer(mpiint), parameter :: OPP_1D_RETCODE = -1_mpiint
  integer(mpiint), parameter :: OPP_TINYASPECT_RETCODE = -2_mpiint

contains

  subroutine init(OPP, comm, skip_load_LUT)
    class(t_optprop), target, intent(inout) :: OPP
    integer(mpiint), intent(in) :: comm
    logical, intent(in), optional :: skip_load_LUT
    integer(mpiint) :: ierr

    select type (OPP)
    class is (t_optprop_1_2)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_1_2 :: OPP%LUT)

    class is (t_optprop_3_6)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_3_6 :: OPP%LUT)

    class is (t_optprop_3_10)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_3_10 :: OPP%LUT)

    class is (t_optprop_8_10)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_8_10 :: OPP%LUT)

    class is (t_optprop_3_16)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_3_16 :: OPP%LUT)

    class is (t_optprop_3_24)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_3_24 :: OPP%LUT)

    class is (t_optprop_3_30)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_3_30 :: OPP%LUT)

    class is (t_optprop_8_16)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_8_16 :: OPP%LUT)

    class is (t_optprop_8_18)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_8_18 :: OPP%LUT)

    class is (t_optprop_wedge_5_8)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_wedge_5_8 :: OPP%LUT)

    class is (t_optprop_rectilinear_wedge_5_8)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_rectilinear_wedge_5_8 :: OPP%LUT)

    class is (t_optprop_wedge_18_8)
      if (.not. allocated(OPP%LUT)) allocate (t_optprop_LUT_wedge_18_8 :: OPP%LUT)

    class is (t_optprop_3_10_ann)
      if (.not. allocated(OPP%ANN)) allocate (t_optprop_ANN_3_10 :: OPP%ANN)

    class default
      call CHKERR(1_mpiint, ' init optprop : unexpected type for optprop object!')
    end select

    if (allocated(OPP%LUT)) then
      call OPP%LUT%init(comm, skip_load_LUT)
      OPP%dev => OPP%LUT
    end if
    if (allocated(OPP%ANN)) then
      call OPP%ANN%init(comm, ierr); call CHKERR(ierr)
      OPP%dev => OPP%ANN
    end if

  end subroutine
  subroutine destroy(OPP, ierr)
    class(t_optprop) :: OPP
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if (allocated(OPP%LUT)) then
      call OPP%LUT%destroy(ierr); call CHKERR(ierr)
      deallocate (OPP%LUT)
    end if
    if (allocated(OPP%ANN)) then
      call OPP%ANN%destroy(ierr); call CHKERR(ierr)
      deallocate (OPP%ANN)
    end if
  end subroutine

  subroutine get_coeff_wedge(OPP, tauz, w0, g, aspect_zx, ldir, C, ierr, wedge_coords, angles)
    class(t_optprop_wedge) :: OPP
    logical, intent(in) :: ldir
    real(irealLUT), intent(in) :: tauz, w0, g, aspect_zx
    real(irealLUT), intent(in) :: wedge_coords(:) ! 2 coordinates of wedge C_point, only used for wedge OPP types
    real(irealLUT), intent(in), optional :: angles(:)
    real(irealLUT), intent(out) :: C(:)
    integer(mpiint), intent(out) :: ierr

    logical, parameter :: compute_coeff_online = .false.

    ierr = 0

    if (ldebug_optprop) then
      call check_inp(OPP, tauz, w0, g, aspect_zx, ldir, C, angles, wedge_coords)
    end if

    if (handle_aspect_zx_1D_case()) return

    if (compute_coeff_online) then
      call do_bmc_computation(C)
      return
    end if

    call do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)

    if (.false. .and. ldir) call print_coeff_diff()

    !if(ldir .and. present(angles)) then
    !  call handle_critical_azimuth()
    !endif
  contains
    subroutine do_bmc_computation(Cbmc)
      real(irealLUT), intent(out) :: Cbmc(:)
      real(ireals) :: vertices(18)
      real(ireal_params) :: phi, theta

      call setup_default_wedge_geometry( &
        real([0, 0], ireals), &
        real([1, 0], ireals), &
        real(wedge_coords, ireals), &
        real(aspect_zx, ireals), &
        vertices, &
        real(wedge_sphere_radius, ireals))

      if (present(angles)) then
        call iterative_phi_theta_from_param_phi_and_param_theta( &
          real(vertices, ireal_params), &
          real(angles(1), ireal_params), &
          real(angles(2), ireal_params), &
          phi, theta, ierr); call CHKERR(ierr)

        call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), ldir, Cbmc, &
                           [rad2deg(real(phi, irealLUT)), rad2deg(real(theta, irealLUT))])
        print *, 'Cbmc', tauz, w0, g, aspect_zx, wedge_coords, ':', angles, rad2deg(phi), rad2deg(theta), '=>', new_line(':'), Cbmc
      else
        call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), ldir, Cbmc)
        print *, 'Cbmc', tauz, w0, g, aspect_zx, wedge_coords, '=>', new_line(':'), Cbmc
      end if

    end subroutine
    subroutine print_coeff_diff()
      real(irealLUT) :: Cbmc(size(C))
      real(ireals) :: err(2)
      integer(iintegers) :: isrc

      call do_bmc_computation(Cbmc)

      err = rmse(real(C, ireals), real(Cbmc, ireals))
      print *, 'rmse', err
      do isrc = 1, OPP%dev%dir_streams
        print *, 'lut src', isrc, ':', C(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams)
        print *, 'bmc src', isrc, ':', Cbmc(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams)
      end do
      if (err(2) .gt. one) then
        !call CHKERR(1_mpiint, 'DEBUG')
      end if
      !C = Cbmc
    end subroutine
    subroutine do_wedge_lookup(tauz, w0, aspect_zx, ldir, angles)
      real(irealLUT), intent(in) :: tauz, w0, aspect_zx
      logical, intent(in) :: ldir
      real(irealLUT), intent(in), optional :: angles(:)
      real(irealLUT) :: save_param_phi, save_param_theta

      real(irealLUT), save, allocatable :: inp_arr_dir(:), inp_arr_diff(:)
      integer(iintegers), save :: dimidx_dir(8)
      integer(iintegers), save :: dimidx_diff(6)
      logical, save :: linit = .false.

      if (.not. linit) then
        dimidx_dir(1) = find_op_dim_by_name(OPP%dev%dirconfig, 'tau')
        dimidx_dir(2) = find_op_dim_by_name(OPP%dev%dirconfig, 'w0')
        dimidx_dir(3) = find_op_dim_by_name(OPP%dev%dirconfig, 'aspect_zx')
        dimidx_dir(4) = find_op_dim_by_name(OPP%dev%dirconfig, 'g')
        dimidx_dir(5) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cx')
        dimidx_dir(6) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cy')
        dimidx_dir(7) = find_op_dim_by_name(OPP%dev%dirconfig, 'param_phi')
        dimidx_dir(8) = find_op_dim_by_name(OPP%dev%dirconfig, 'param_theta')
        allocate (inp_arr_dir(count(dimidx_dir .gt. 0)))

        dimidx_diff(1) = find_op_dim_by_name(OPP%dev%dirconfig, 'tau')
        dimidx_diff(2) = find_op_dim_by_name(OPP%dev%dirconfig, 'w0')
        dimidx_diff(3) = find_op_dim_by_name(OPP%dev%dirconfig, 'aspect_zx')
        dimidx_diff(4) = find_op_dim_by_name(OPP%dev%dirconfig, 'g')
        dimidx_diff(5) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cx')
        dimidx_diff(6) = find_op_dim_by_name(OPP%dev%dirconfig, 'wedge_coord_Cy')
        allocate (inp_arr_diff(count(dimidx_diff .gt. 0)))
        linit = .true.
      end if

      if (present(angles)) then ! obviously we want the direct coefficients

        associate ( &
          param_phi => angles(1), &
          param_theta => angles(2))

          call handle_critical_param_phi(param_phi, save_param_phi)
          call handle_critical_param_theta(param_theta, save_param_theta)

          if (dimidx_dir(1) .gt. 0) inp_arr_dir(dimidx_dir(1)) = tauz
          if (dimidx_dir(2) .gt. 0) inp_arr_dir(dimidx_dir(2)) = w0
          if (dimidx_dir(3) .gt. 0) inp_arr_dir(dimidx_dir(3)) = aspect_zx
          if (dimidx_dir(4) .gt. 0) inp_arr_dir(dimidx_dir(4)) = g
          if (dimidx_dir(5) .gt. 0) inp_arr_dir(dimidx_dir(5)) = wedge_coords(1)
          if (dimidx_dir(6) .gt. 0) inp_arr_dir(dimidx_dir(6)) = wedge_coords(2)
          if (dimidx_dir(7) .gt. 0) inp_arr_dir(dimidx_dir(7)) = save_param_phi
          if (dimidx_dir(8) .gt. 0) inp_arr_dir(dimidx_dir(8)) = save_param_theta

          if (ldir) then ! dir2dir
            call OPP%dev%get_dir2dir(inp_arr_dir, C)
          else ! dir2diff
            call OPP%dev%get_dir2diff(inp_arr_dir, C)
          end if
        end associate
      else
        ! diff2diff
        if (dimidx_diff(1) .gt. 0) inp_arr_diff(dimidx_diff(1)) = tauz
        if (dimidx_diff(2) .gt. 0) inp_arr_diff(dimidx_diff(2)) = w0
        if (dimidx_diff(3) .gt. 0) inp_arr_diff(dimidx_diff(3)) = aspect_zx
        if (dimidx_diff(4) .gt. 0) inp_arr_diff(dimidx_diff(4)) = g
        if (dimidx_diff(5) .gt. 0) inp_arr_diff(dimidx_diff(5)) = wedge_coords(1)
        if (dimidx_diff(6) .gt. 0) inp_arr_diff(dimidx_diff(6)) = wedge_coords(2)
        call OPP%dev%get_diff2diff(inp_arr_diff, C)
      end if

    end subroutine

    subroutine handle_critical_param_phi(param_phi, save_param_phi)
      real(irealLUT), intent(in) :: param_phi
      real(irealLUT), intent(out) :: save_param_phi
      logical :: lsample_critical

      lsample_critical = .false.
      if (approx(abs(param_phi), 1._ireallut, param_eps)) lsample_critical = .true.

      if (lsample_critical) then
        if (param_phi .lt. -1._ireallut) then
          save_param_phi = -1._ireallut - param_eps
        elseif (param_phi .ge. 1._ireallut) then !1.0001
          save_param_phi = 1._ireallut + param_eps
        elseif (param_phi .lt. 0._ireallut) then !-.999
          save_param_phi = -1._ireallut + param_eps
        else ! .999
          save_param_phi = 1._ireallut - param_eps
        end if
      else
        save_param_phi = param_phi
      end if
    end subroutine
    subroutine handle_critical_param_theta(param_theta, save_param_theta)
      real(irealLUT), intent(in) :: param_theta
      real(irealLUT), intent(out) :: save_param_theta
      logical :: lsample_critical

      lsample_critical = .false.
      if (approx(abs(param_theta), 0._ireallut, param_eps)) lsample_critical = .true.

      if (lsample_critical) then
        if (param_theta .le. 0._ireallut) then
          save_param_theta = -param_eps
        else ! .0001
          save_param_theta = param_eps
        end if
      else
        save_param_theta = param_theta
      end if
    end subroutine

    logical function handle_aspect_zx_1D_case()
      real(ireals) :: c11, c12, c13, c23, c33
      real(irealLUT) :: restricted_aspect_zx
      real(irealLUT) :: mu

      handle_aspect_zx_1D_case = .false.

      !TODO: here we may have incoming radiation at the sides and we just drop that
      ! this has to be fixed for anisotropic grids

      if (present(angles)) then

        if (aspect_zx .ge. twostr_ratio) then
          C = zero
          mu = real(cos(theta_from_param_theta(real(angles(2), ireal_params), 0._ireal_params)), irealLUT)

          call eddington_coeff_ec( &
            real(tauz, ireals), &
            real(w0, ireals), &
            real(g, ireals), &
            real(mu, ireals), &
            c11, c12, c13, c23, c33)

          if (ldir) then
            select type (OPP)
            class is (t_optprop_wedge_5_8)
              ! set the transport coeffs for src top to zero, leave the rest.
              C(5 * 4 + 1) = real(c33, irealLUT) ! from top to bot
              !C(22:24) = 1 ! from sides to bot
            class is (t_optprop_rectilinear_wedge_5_8)
              ! set the transport coeffs for src top to zero, leave the rest.
              C(5 * 4 + 1) = real(c33, irealLUT) ! from top to bot
              !C(22:24) = 1 ! from sides to bot
            class is (t_optprop_wedge_18_8)
              C(18 * 15 + 1) = real(c33, irealLUT) ! from top to bot
              C(18 * 16 + 2) = real(c33, irealLUT) ! from top to bot
              C(18 * 17 + 3) = real(c33, irealLUT) ! from top to bot
            class default
              call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
            end select

          else
            select type (OPP)
            class is (t_optprop_wedge_5_8)
              C(0 * 5 + 1) = real(c13, irealLUT) ! reflection
              C(7 * 5 + 1) = real(c23, irealLUT) ! transmission
            class is (t_optprop_rectilinear_wedge_5_8)
              C(0 * 5 + 1) = real(c13, irealLUT) ! reflection
              C(7 * 5 + 1) = real(c23, irealLUT) ! transmission
            class is (t_optprop_wedge_18_8)
              C(0 * 18 + 1) = real(c13, irealLUT)
              C(0 * 18 + 2) = real(c13, irealLUT)
              C(0 * 18 + 3) = real(c13, irealLUT)
              C(7 * 18 + 1) = real(c23, irealLUT)
              C(7 * 18 + 2) = real(c23, irealLUT)
              C(7 * 18 + 3) = real(c23, irealLUT)
            class default
              call CHKERR(1_mpiint, 'wedge handle_aspect_zx_1D_case not implemented for this type of OPP')
            end select
          end if

          handle_aspect_zx_1D_case = .true.
          ierr = OPP_1D_RETCODE

        elseif (aspect_zx .lt. OPP%dev%dirconfig%dims(3)%vrange(1)) then
          restricted_aspect_zx = min(max(aspect_zx, OPP%dev%dirconfig%dims(3)%vrange(1)), &
                                     OPP%dev%dirconfig%dims(3)%vrange(2))
          call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
          handle_aspect_zx_1D_case = .true.
          ierr = OPP_TINYASPECT_RETCODE
        end if

      else ! diffuse

        !if(aspect_zx.gt.OPP%dev%diffconfig%dims(3)%vrange(2)) then
        if (aspect_zx .ge. twostr_ratio) then
          C = zero

          call eddington_coeff_ec( &
            real(tauz, ireals), &
            real(w0, ireals), &
            real(g, ireals), &
            one, &
            c11, c12, c13, c23, c33)

          ! transmission & reflection towards top plate
          C(1) = real(c12, irealLUT)
          C(8) = real(c11, irealLUT)
          ! and bot plate
          C(7 * 8 + 1) = real(c11, irealLUT)
          C(8 * 8) = real(c12, irealLUT)

          handle_aspect_zx_1D_case = .true.
          ierr = OPP_1D_RETCODE

        elseif (aspect_zx .lt. OPP%dev%diffconfig%dims(3)%vrange(1)) then
          restricted_aspect_zx = min(max(aspect_zx, OPP%dev%diffconfig%dims(3)%vrange(1)), &
                                     OPP%dev%diffconfig%dims(3)%vrange(2))
          call do_wedge_lookup(tauz, w0, restricted_aspect_zx, ldir, angles)
          handle_aspect_zx_1D_case = .true.
          ierr = OPP_TINYASPECT_RETCODE
        end if

      end if
    end function
  end subroutine

  subroutine get_coeff_cube(OPP, tauz, w0, g, aspect_zx, dir, C, ierr, angles, lswitch_east, lswitch_north, opt_vertices)
    class(t_optprop_cube) :: OPP
    logical, intent(in) :: dir
    real(irealLUT), intent(in) :: tauz, w0, g, aspect_zx
    real(irealLUT), intent(in), optional :: angles(:)
    logical, intent(in), optional :: lswitch_east, lswitch_north
    real(irealLUT), intent(out) :: C(:)
    integer(mpiint), intent(out) :: ierr
    real(ireals), intent(in), optional :: opt_vertices(:)

    logical, save :: compute_coeff_online = .false., lset = .false., bmc_default_unit_cube_reference = .false.
    logical :: lflg
    real(ireals) :: vertices(24)
    real(irealLUT), allocatable :: Clut(:), Cbmc(:), Cbmc2(:)
    real(irealLUT) :: save_aspect_zx
    ierr = 0

    if (ldebug_optprop) call check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles)

    if (present(angles)) then ! obviously we want the direct coefficients
      if (aspect_zx .lt. OPP%dev%dirconfig%dims(3)%vrange(1)) then
        save_aspect_zx = OPP%dev%dirconfig%dims(3)%vrange(1)
        ierr = OPP_TINYASPECT_RETCODE
      else
        save_aspect_zx = aspect_zx
      end if
      if (dir) then ! dir2dir
        call OPP%dev%get_dir2dir([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
        call OPP%dir2dir_coeff_symmetry(C, lswitch_east, lswitch_north)
      else         ! dir2diff
        call OPP%dev%get_dir2diff([tauz, w0, save_aspect_zx, g, angles(1), angles(2)], C)
        call OPP%dir2diff_coeff_symmetry(C, lswitch_east, lswitch_north)
      end if
    else
      ! diff2diff
      if (aspect_zx .lt. OPP%dev%diffconfig%dims(3)%vrange(1)) then
        save_aspect_zx = OPP%dev%diffconfig%dims(3)%vrange(1)
        ierr = OPP_TINYASPECT_RETCODE
      else
        save_aspect_zx = aspect_zx
      end if
      call OPP%dev%get_diff2diff([tauz, w0, save_aspect_zx, g], C)
      !call OPP%diff2diff_coeff_symmetry(C)
    end if

    if (.not. lset) then
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-bmc_online", compute_coeff_online, lflg, ierr); call CHKERR(ierr)
      if (present(opt_vertices)) then
        call get_petsc_opt(PETSC_NULL_CHARACTER, "-bmc_default_unit_cube_reference", bmc_default_unit_cube_reference, lflg, ierr)
        call CHKERR(ierr)
      end if
      lset = .true.
    end if

    if (compute_coeff_online) then
      allocate (Clut(size(C)), Cbmc(size(C)))
      Clut = C

      if (.not. present(opt_vertices) .or. bmc_default_unit_cube_reference) then
        call setup_default_unit_cube_geometry(one, one, real(aspect_zx, ireals), vertices)
        call get_coeff_bmc(OPP, vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, Cbmc, angles)
        C = Cbmc
      end if
      if (present(opt_vertices)) then
        allocate (Cbmc2(size(C)))
        call get_coeff_bmc(OPP, opt_vertices, real(tauz, ireals), real(w0, ireals), real(g, ireals), dir, Cbmc2, angles)
        C = Cbmc2
        if (bmc_default_unit_cube_reference) then
          if (present(angles)) then
            print *, new_line(''), opt_vertices(3:24:3), ':', angles, new_line('')// &
              cstr('LUT            '//toStr(Clut), 'black')//new_line('')// &
              cstr('bmc (regular  )'//toStr(Cbmc), 'blue')//new_line('')// &
              cstr('bmc (distorted)'//toStr(Cbmc2), 'green')
          else
            print *, new_line(''), opt_vertices(3:24:3), new_line('')// &
              cstr('LUT            '//toStr(Clut), 'black')//new_line('')// &
              cstr('bmc (regular  )'//toStr(Cbmc), 'blue')//new_line('')// &
              cstr('bmc (distorted)'//toStr(Cbmc2), 'green')
          end if
        else
          if (present(angles)) then
            print *, new_line(''), opt_vertices(3:24:3), ':', angles, new_line('')// &
              cstr('LUT            '//toStr(Clut), 'black')//new_line('')// &
              cstr('bmc (distorted)'//toStr(Cbmc2), 'green')
          else
            print *, new_line(''), opt_vertices(3:24:3), new_line('')// &
              cstr('LUT            '//toStr(Clut), 'black')//new_line('')// &
              cstr('bmc (distorted)'//toStr(Cbmc2), 'green')
          end if
        end if
      end if
    end if
  end subroutine

  subroutine get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, angles)
    class(t_optprop) :: OPP
    real(ireals), intent(in) :: tauz, w0, g, vertices(:)
    logical, intent(in) :: dir
    real(irealLUT), intent(out) :: C(:)
    real(irealLUT), intent(in), optional :: angles(2)

    real(irealLUT) :: S_diff(OPP%dev%diff_streams), T_dir(OPP%dev%dir_streams)
    real(irealLUT) :: S_tol(OPP%dev%diff_streams), T_tol(OPP%dev%dir_streams)
    integer(iintegers) :: isrc

    real(irealLUT), parameter :: atol = 1e-3_ireallut, rtol = 5e-1_ireallut

    if (present(angles)) then
      do isrc = 1, OPP%dev%dir_streams
        call OPP%dev%bmc_wrapper(isrc, &
                                 real(vertices, ireal_dp), &
                                 real(tauz, irealLUT), &
                                 real(w0, irealLUT), &
                                 real(g, irealLUT), &
                                 .true., &
                                 real(angles(1), irealLUT), &
                                 real(angles(2), irealLUT), &
                                 mpi_comm_self, &
                                 S_diff, T_dir, S_tol, T_tol, &
                                 inp_atol=real(atol, irealLUT), inp_rtol=real(rtol, irealLUT))
        if (dir) then !dir2dir
          C(isrc:OPP%dev%dir_streams**2:OPP%dev%dir_streams) = T_dir
        else ! dir2diff
          C(isrc:OPP%dev%dir_streams * OPP%dev%diff_streams:OPP%dev%dir_streams) = S_diff
        end if
        if (w0 .ge. 1) then
          if (any(S_tol .gt. 0) .or. any(T_tol .gt. 0)) then
            print *, 'SumT', sum(T_dir), 'SumS', sum(S_diff), &
              'Divergence', (1 - (sum(T_dir) + sum(S_diff))), any(S_tol .gt. 0), any(T_tol .gt. 0)
            if (abs(1 - (sum(T_dir) + sum(S_diff))) .ge. 1e-6_ireallut) then
              call CHKWARN(1_mpiint, 'divergence '// &
                           toStr(1 - (sum(T_dir) + sum(S_diff)))//' seems quite large for w0='//toStr(w0))
            end if
          end if
        end if
      end do
    else
      ! diff2diff
      do isrc = 1, OPP%dev%diff_streams
        call OPP%dev%bmc_wrapper(isrc, &
                                 real(vertices, ireal_dp), &
                                 real(tauz, irealLUT), &
                                 real(w0, irealLUT), &
                                 real(g, irealLUT), &
                                 .false., &
                                 0._ireallut, 0._ireallut, &
                                 mpi_comm_self, &
                                 S_diff, T_dir, S_tol, T_tol, &
                                 inp_atol=real(atol, irealLUT), inp_rtol=real(rtol, irealLUT))
        C(isrc:OPP%dev%diff_streams**2:OPP%dev%diff_streams) = S_diff
        if (w0 .ge. 1) then
          if (any(S_tol .gt. 0) .or. any(T_tol .gt. 0)) then
            print *, 'SumT', sum(T_dir), 'SumS', sum(S_diff), &
              'Divergence', (1 - (sum(T_dir) + sum(S_diff))), any(S_tol .gt. 0), any(T_tol .gt. 0)
            if (abs(1 - (sum(T_dir) + sum(S_diff))) .ge. 1e-6_ireallut) then
              call CHKWARN(1_mpiint, 'divergence '// &
                           toStr(1 - (sum(T_dir) + sum(S_diff)))//' seems quite large for w0='//toStr(w0))
            end if
          end if
        end if
      end do
    end if ! angles_present

  end subroutine

  subroutine check_inp(OPP, tauz, w0, g, aspect_zx, dir, C, angles, wedge_coords)
    class(t_optprop) :: OPP
    real(irealLUT), intent(in) :: tauz, w0, g, aspect_zx
    logical, intent(in) :: dir
    real(irealLUT), intent(in) :: C(:)
    real(irealLUT), intent(in), optional :: angles(:), wedge_coords(:)
    integer(iintegers), save :: dimidx_dir(8) = -2
    integer(iintegers), save :: dimidx_diff(6) = -2

    if ((any([tauz, w0, g, aspect_zx] .lt. zero)) .or. (any(isnan([tauz, w0, g, aspect_zx])))) then
      call CHKERR(1_mpiint, 'optprop_lookup_coeff :: '// &
                  'corrupt optical properties: bg:: '//toStr([tauz, w0, g, aspect_zx]))
    end if

    if (dir) then
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['tau'], tauz, dimidx_dir(1))
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['w0'], w0, dimidx_dir(2))
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['g'], g, dimidx_dir(3), default_val=0._ireallut)
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['aspect_zx'], aspect_zx, dimidx_dir(4))
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['phi        ', 'param_phi  '], angles(1), dimidx_dir(5))
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['theta      ', 'param_theta'], angles(2), dimidx_dir(6))
    else
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['tau'], tauz, dimidx_diff(1))
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['w0'], w0, dimidx_diff(2))
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['g'], g, dimidx_diff(3), default_val=0._ireallut)
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['aspect_zx'], aspect_zx, dimidx_diff(4))
    end if

    if (present(angles)) then
      if (dir .and. size(C) .ne. OPP%dev%dir_streams**2) then
        print *, 'direct called get_coeff with wrong shaped output array:', size(C), 'should be ', OPP%dev%dir_streams**2
      end if
      if (.not. dir .and. size(C) .ne. OPP%dev%diff_streams * OPP%dev%dir_streams) then
        print *, 'dir2diffuse called get_coeff with wrong shaped output array:', size(C), &
          'should be', OPP%dev%diff_streams * OPP%dev%dir_streams
      end if
    else
      if (dir .and. size(C) .ne. OPP%dev%diff_streams) then
        print *, 'diff2diff called get_coeff with wrong shaped output array:', size(C), 'should be ', OPP%dev%diff_streams
      end if
      if (.not. dir .and. size(C) .ne. OPP%dev%diff_streams**2) then
        print *, 'diff2diff called get_coeff with wrong shaped output array:', size(C), 'should be ', OPP%dev%diff_streams**2
      end if
    end if

    if (present(wedge_coords)) then
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['wedge_coord_Cx'], wedge_coords(1), dimidx_dir(7))
      call check_LUT_dimension_limits(OPP%dev%dirconfig, ['wedge_coord_Cy'], wedge_coords(2), dimidx_dir(8))
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['wedge_coord_Cx'], wedge_coords(1), dimidx_diff(5))
      call check_LUT_dimension_limits(OPP%dev%diffconfig, ['wedge_coord_Cy'], wedge_coords(2), dimidx_diff(6))
    end if
  contains
    subroutine check_LUT_dimension_limits(config, dimnames, val, dimindex, default_val)
      type(t_op_config), intent(in) :: config
      character(len=*), intent(in) :: dimnames(:)
      real(irealLUT), intent(in) :: val
      integer(iintegers), intent(inout) :: dimindex
      real(irealLUT), intent(in), optional :: default_val
      integer(iintegers) :: i
      if (dimindex .lt. -1) then
        do i = 1, size(dimnames)
          if (dimindex .le. -1) then
            dimindex = find_op_dim_by_name(config, trim(dimnames(i)))
          else
            exit
          end if
          !print *,'Looking for dim: '//trim(dimnames(i)),' -> ', dimindex
        end do
      end if
      if (dimindex .gt. -1) then
        if (.not. is_inrange(val, &
                             config%dims(dimindex)%vrange(1) - tiny(val), &
                             config%dims(dimindex)%vrange(2) + tiny(val))) then
          call CHKERR(1_mpiint, 'value ('//toStr(val)//') is not in the range of the LUT for dimension '// &
                      trim(config%dims(dimindex)%dimname)//' ( '//toStr(config%dims(dimindex)%vrange)//' )')
        end if
      else if (dimindex .eq. -1) then
        if (.not. present(default_val)) then
          print *, 'Available LUT Dimensions are: '//char_arr_to_str(config%dims(:)%dimname, ', ')
          call CHKERR(1_mpiint, 'currently the LUT calls do not have a dimension for '// &
                      char_arr_to_str(dimnames, ','))
        else
          if (.not. approx(val, default_val)) then
            print *, 'Available LUT Dimensions are: '//char_arr_to_str(config%dims(:)%dimname, ', ')
            call CHKERR(1_mpiint, 'currently the LUT calls do not have a dimension for '// &
                        char_arr_to_str(dimnames, ',')// &
                        ' and should probably be limited to the default val: '//toStr(default_val))
          end if
        end if
      end if
    end subroutine
  end subroutine

  subroutine diff2diff_coeff_symmetry(OPP, coeff)
    class(t_optprop_cube) :: OPP
    real(irealLUT), target, intent(inout) :: coeff(:)
    real(irealLUT), pointer :: v(:, :) ! dim(src, dst)
    integer(iintegers) :: i
    real(irealLUT) :: norm(0:9)

    if (OPP%dev%diff_streams .eq. 10) then
      v(0:9, 0:9) => coeff(1:100)

      do i = 0, 9
        norm(i) = sum(v(i, :))
      end do
      norm(0:1) = sum(norm(0:1)) / 2
      norm(2:5) = sum(norm(2:5)) / 4
      norm(6:9) = sum(norm(6:9)) / 4

      v(0, 0) = (v(0, 0) + v(1, 1))*.5_ireallut
      v(1, 1) = v(0, 0)

      v(0, 1) = (v(0, 1) + v(1, 0))*.5_ireallut
      v(1, 0) = v(0, 1)

      v(1, 0) = (v(0, 1) + v(1, 0))*.5_ireallut
      v(0, 1) = v(1, 0)

      ! v(0,2) = v(0,3) = v(0,6) = v(0,7) = v(1,4) = v(1,5) = v(1,8) = v(1,9) = np.mean(v(0,(2,3,6,7)) + v(1,(4,5,8,9)), axis=-1)/2
      v(0, 2) = (v(0, 2) + v(0, 3) + v(0, 6) + v(0, 7) + &
        & v(1, 4) + v(1, 5) + v(1, 8) + v(1, 9)) / 8._ireallut
      v(0, 3) = v(0, 2)
      v(0, 6) = v(0, 2)
      v(0, 7) = v(0, 2)
      v(1, 4) = v(0, 2)
      v(1, 5) = v(0, 2)
      v(1, 8) = v(0, 2)
      v(1, 9) = v(0, 2)

      !v(1,2) = v(1,3) = v(1,6) = v(1,7) = v(0,4) = v(0,5) = v(0,8) = v(0,9) = np.mean(v(1,(2,3,6,7)) + v(0,(4,5,8,9)), axis=-1)/2
      v(1, 2) = (v(0, 4) + v(0, 5) + v(0, 8) + v(0, 9) + &
        & v(1, 2) + v(1, 3) + v(1, 6) + v(1, 7)) / 8._ireallut
      v(1, 3) = v(1, 2)
      v(1, 6) = v(1, 2)
      v(1, 7) = v(1, 2)
      v(0, 4) = v(1, 2)
      v(0, 5) = v(1, 2)
      v(0, 8) = v(1, 2)
      v(0, 9) = v(1, 2)

      !v(2,0) = v(3,0) = v(4,1) = v(5,1) = v(6,0) = v(7,0) = v(8,1) = v(9,1) = np.mean(v((2,3,6,7),0) + v((4,5,8,9),1), axis=-1)/2
      v(2, 0) = (v(2, 0) + v(3, 0) + v(6, 0) + v(7, 0) + &
        & v(4, 1) + v(5, 1) + v(8, 1) + v(9, 1)) / 8._ireallut
      v(3, 0) = v(2, 0)
      v(6, 0) = v(2, 0)
      v(7, 0) = v(2, 0)
      v(4, 1) = v(2, 0)
      v(5, 1) = v(2, 0)
      v(8, 1) = v(2, 0)
      v(9, 1) = v(2, 0)

      !v(2,1) = v(3,1) = v(4,0) = v(5,0) = v(6,1) = v(7,1) = v(8,0) = v(9,0) = np.mean(v((2,3,6,7),1) + v((4,5,8,9),0), axis=-1)/2
      v(2, 1) = (v(4, 0) + v(5, 0) + v(8, 0) + v(9, 0) + &
        & v(2, 1) + v(3, 1) + v(6, 1) + v(7, 1)) / 8._ireallut
      v(4, 0) = v(2, 1)
      v(5, 0) = v(2, 1)
      v(8, 0) = v(2, 1)
      v(9, 0) = v(2, 1)
      v(3, 1) = v(2, 1)
      v(6, 1) = v(2, 1)
      v(7, 1) = v(2, 1)

      !v(2,2) = v(3,3) = v(4,4) = v(5,5) = v(6,6) = v(7,7) = v(8,8) = v(9,9) = np.mean(( v(_,_) for _ in range(2,10)), axis=0)
      v(2, 2) = (v(2, 2) + v(3, 3) + v(4, 4) + v(5, 5) + v(6, 6) + v(7, 7) + v(8, 8) + v(9, 9)) / 8._ireallut
      v(3, 3) = v(2, 2)
      v(4, 4) = v(2, 2)
      v(5, 5) = v(2, 2)
      v(6, 6) = v(2, 2)
      v(7, 7) = v(2, 2)
      v(8, 8) = v(2, 2)
      v(9, 9) = v(2, 2)

      !v(2,3) = v(3,2) = v(4,5) = v(5,4) = v(6,7) = v(7,6) = v(8,9) = v(9,8) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(3,2,5,4,7,6,9,8))), axis=0)
      v(2, 3) = (v(2, 3) + v(3, 2) + v(4, 5) + v(5, 4) + v(6, 7) + v(7, 6) + v(8, 9) + v(9, 8)) / 8._ireallut
      v(3, 2) = v(2, 3)
      v(4, 5) = v(2, 3)
      v(5, 4) = v(2, 3)
      v(6, 7) = v(2, 3)
      v(7, 6) = v(2, 3)
      v(8, 9) = v(2, 3)
      v(9, 8) = v(2, 3)

      !v(2,4) = v(3,5) = v(4,2) = v(5,3) = v(6,8) = v(7,9) = v(8,6) = v(9,7) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(4,5,2,3,8,9,6,7))), axis=0)
      v(2, 4) = (v(2, 4) + v(3, 5) + v(4, 2) + v(5, 3) + v(6, 8) + v(7, 9) + v(8, 6) + v(9, 7)) / 8._ireallut
      v(3, 5) = v(2, 4)
      v(4, 2) = v(2, 4)
      v(5, 3) = v(2, 4)
      v(6, 8) = v(2, 4)
      v(7, 9) = v(2, 4)
      v(8, 6) = v(2, 4)
      v(9, 7) = v(2, 4)

      !v(2,6) = v(2,7) = v(3,6) = v(3,7) = v(4,8) = v(4,9) = v(5,8) = v(5,9) = v(6,2) = v(6,3) = v(7,2) = v(7,3) = v(8,4) = v(8,5) = v(9,4) = v(9,5) = np.mean(( v(i1,i2) for i1,i2 in zip(sorted(list(range(2,10))*2),(6,7,6,7, 8,9,8,9, 2,3,2,3, 4,5,4,5,))), axis=0)
      v(2, 6) = (v(2, 6) + v(2, 7) + v(3, 6) + v(3, 7) + &
        & v(4, 8) + v(4, 9) + v(5, 8) + v(5, 9) + &
        & v(6, 2) + v(6, 3) + v(7, 2) + v(7, 3) + &
        & v(8, 4) + v(8, 5) + v(9, 4) + v(9, 5)) / 16._ireallut
      v(2, 7) = v(2, 6)
      v(3, 6) = v(2, 6)
      v(3, 7) = v(2, 6)
      v(4, 8) = v(2, 6)
      v(4, 9) = v(2, 6)
      v(5, 8) = v(2, 6)
      v(5, 9) = v(2, 6)
      v(6, 2) = v(2, 6)
      v(6, 3) = v(2, 6)
      v(7, 2) = v(2, 6)
      v(7, 3) = v(2, 6)
      v(8, 4) = v(2, 6)
      v(8, 5) = v(2, 6)
      v(9, 4) = v(2, 6)
      v(9, 5) = v(2, 6)

      !v(2,8) = v(2,9) = v(3,8) = v(3,9) = v(4,6) = v(4,7) = v(5,6) = v(5,7) = v(6,4) = v(6,5) = v(7,4) = v(7,5) = v(8,2) = v(8,3) = v(9,2) = v(9,3) = np.mean(( v(i1,i2) for i1,i2 in zip(sorted(list(range(2,10))*2),(8,9,8,9, 6,7,6,7, 4,5,4,5, 2,3,2,3,))), axis=0)
      v(2, 8) = (v(2, 8) + v(2, 9) + v(3, 8) + v(3, 9) + &
        & v(4, 6) + v(4, 7) + v(5, 6) + v(5, 7) + &
        & v(6, 4) + v(6, 5) + v(7, 4) + v(7, 5) + &
        & v(8, 2) + v(8, 3) + v(9, 2) + v(9, 3)) / 16._ireallut
      v(2, 9) = v(2, 8)
      v(3, 8) = v(2, 8)
      v(3, 9) = v(2, 8)
      v(4, 6) = v(2, 8)
      v(4, 7) = v(2, 8)
      v(5, 6) = v(2, 8)
      v(5, 7) = v(2, 8)
      v(6, 4) = v(2, 8)
      v(6, 5) = v(2, 8)
      v(7, 4) = v(2, 8)
      v(7, 5) = v(2, 8)
      v(8, 2) = v(2, 8)
      v(8, 3) = v(2, 8)
      v(9, 2) = v(2, 8)
      v(9, 3) = v(2, 8)

      !v(2,5) = v(3,4) = v(4,3) = v(5,2) = v(6,9) = v(7,8) = v(8,7) = v(9,6) = np.mean(( v(i1,i2) for i1,i2 in zip(range(2,10),(5,4,3,2,9,8,7,6))), axis=0)
      v(2, 5) = (v(2, 5) + v(3, 4) + v(4, 3) + v(5, 2) + v(6, 9) + v(7, 8) + v(8, 7) + v(9, 6)) / 8._ireallut
      v(3, 4) = v(2, 5)
      v(4, 3) = v(2, 5)
      v(5, 2) = v(2, 5)
      v(6, 9) = v(2, 5)
      v(7, 8) = v(2, 5)
      v(8, 7) = v(2, 5)
      v(9, 6) = v(2, 5)

      do i = 0, 9
        v(i, :) = v(i, :) / max(tiny(v), sum(v(i, :))) * norm(i)
      end do
    end if

    if (.false.) then
      select type (OPP)
      end select
    end if
  end subroutine

  subroutine dir2diff_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    return
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
      if (lswitch_east .or. lswitch_north) coeff = coeff
    end if
  end subroutine

  !for solver_3_6 only the offset is changing in those sides which should be switched
  subroutine dir3_to_diff6_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_6) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 3
    real(irealLUT) :: newcoeff(size(coeff))
    if (lswitch_east) then
      newcoeff = coeff
      coeff(7:9) = newcoeff([1, 2, 3] + dof * 3)
      coeff(10:12) = newcoeff([1, 2, 3] + dof * 2)
    end if
    if (lswitch_north) then
      newcoeff = coeff
      coeff(13:15) = newcoeff([1, 2, 3] + dof * 5)
      coeff(16:18) = newcoeff([1, 2, 3] + dof * 4)
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
    end if
  end subroutine

  !for solver_3_10 the offset is chaning and the destination order
  subroutine dir3_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_10) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 3
    real(irealLUT) :: newcoeff(size(coeff))
    if (lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
      coeff(1 + (3 - 1) * dof:3 * dof) = newcoeff([1, 2, 3] + dof * (4 - 1))
      coeff(1 + (4 - 1) * dof:4 * dof) = newcoeff([1, 2, 3] + dof * (3 - 1))
      coeff(1 + (5 - 1) * dof:5 * dof) = newcoeff([1, 2, 3] + dof * (6 - 1))
      coeff(1 + (6 - 1) * dof:6 * dof) = newcoeff([1, 2, 3] + dof * (5 - 1))
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
    end if
    if (lswitch_north) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
      coeff(1 + (7 - 1) * dof:7 * dof) = newcoeff([1, 2, 3] + dof * (8 - 1))
      coeff(1 + (8 - 1) * dof:8 * dof) = newcoeff([1, 2, 3] + dof * (7 - 1))
      coeff(1 + (9 - 1) * dof:9 * dof) = newcoeff([1, 2, 3] + dof * (10 - 1))
      coeff(1 + (10 - 1) * dof:10 * dof) = newcoeff([1, 2, 3] + dof * (9 - 1))
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
    end if
  end subroutine

  subroutine dir3_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_16) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 3
    real(irealLUT) :: newcoeff(size(coeff))
    if (lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([1, 2, 3] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([1, 2, 3] + dof*( 2-1) )
      coeff(1 + (3 - 1) * dof:3 * dof) = newcoeff([1, 2, 3] + dof * (7 - 1))
      coeff(1 + (4 - 1) * dof:4 * dof) = newcoeff([1, 2, 3] + dof * (8 - 1))
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([1, 2, 3] + dof*( 5-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([1, 2, 3] + dof*( 6-1) )
      coeff(1 + (7 - 1) * dof:7 * dof) = newcoeff([1, 2, 3] + dof * (3 - 1))
      coeff(1 + (8 - 1) * dof:8 * dof) = newcoeff([1, 2, 3] + dof * (4 - 1))
      coeff(1 + (9 - 1) * dof:9 * dof) = newcoeff([1, 2, 3] + dof * (10 - 1))
      coeff(1 + (10 - 1) * dof:10 * dof) = newcoeff([1, 2, 3] + dof * (9 - 1))
      coeff(1 + (11 - 1) * dof:11 * dof) = newcoeff([1, 2, 3] + dof * (12 - 1))
      coeff(1 + (12 - 1) * dof:12 * dof) = newcoeff([1, 2, 3] + dof * (11 - 1))
      !coeff(1+(13-1)*dof:13*dof) = newcoeff([1, 2, 3] + dof*(13-1) )
      !coeff(1+(14-1)*dof:14*dof) = newcoeff([1, 2, 3] + dof*(14-1) )
      !coeff(1+(15-1)*dof:15*dof) = newcoeff([1, 2, 3] + dof*(15-1) )
      !coeff(1+(16-1)*dof:16*dof) = newcoeff([1, 2, 3] + dof*(16-1) )
    end if
    if (lswitch_north) then
      newcoeff = coeff
      coeff(1 + (1 - 1) * dof:1 * dof) = newcoeff([1, 2, 3] + dof * (5 - 1))
      coeff(1 + (2 - 1) * dof:2 * dof) = newcoeff([1, 2, 3] + dof * (6 - 1))
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([1, 2, 3] + dof*( 3-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([1, 2, 3] + dof*( 4-1) )
      coeff(1 + (5 - 1) * dof:5 * dof) = newcoeff([1, 2, 3] + dof * (1 - 1))
      coeff(1 + (6 - 1) * dof:6 * dof) = newcoeff([1, 2, 3] + dof * (2 - 1))
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([1, 2, 3] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([1, 2, 3] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([1, 2, 3] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([1, 2, 3] + dof*(10-1) )
      !coeff(1+(11-1)*dof:11*dof) = newcoeff([1, 2, 3] + dof*(11-1) )
      !coeff(1+(12-1)*dof:12*dof) = newcoeff([1, 2, 3] + dof*(12-1) )
      coeff(1 + (13 - 1) * dof:13 * dof) = newcoeff([1, 2, 3] + dof * (14 - 1))
      coeff(1 + (14 - 1) * dof:14 * dof) = newcoeff([1, 2, 3] + dof * (13 - 1))
      coeff(1 + (15 - 1) * dof:15 * dof) = newcoeff([1, 2, 3] + dof * (16 - 1))
      coeff(1 + (16 - 1) * dof:16 * dof) = newcoeff([1, 2, 3] + dof * (15 - 1))
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
      if (lswitch_east .or. lswitch_north) coeff = coeff
    end if
  end subroutine

  subroutine dir3_to_diff24_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_24) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    !integer(iintegers), parameter:: dof = 3
    real(irealLUT), target :: bak(3, 24)
    real(irealLUT), pointer :: pC(:, :) ! dim(src,dst)
    if (lswitch_east) then
      pC(1:3, 1:24) => coeff
      bak(:, :) = pC(:, :)
      pC(:, [1, 3, 5, 7, 2, 4, 6, 8]) = bak(:, [3, 1, 7, 5, 4, 2, 8, 6])
      pC(:, 10:16:2) = bak(:, 9:15:2)
      pC(:, 9:15:2) = bak(:, 10:16:2)
    end if
    if (lswitch_north) then
      pC(1:3, 1:24) => coeff
      bak(:, :) = pC(:, :)
      pC(:, [1, 3, 5, 7, 2, 4, 6, 8]) = bak(:, [5, 7, 1, 3, 6, 8, 2, 4])
      pC(:, 18:24:2) = bak(:, 17:24:2)
      pC(:, 17:24:2) = bak(:, 18:24:2)
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
      if (lswitch_east .or. lswitch_north) coeff = coeff
    end if
  end subroutine

  subroutine dir3_to_diff30_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_3_30) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    !integer(iintegers), parameter:: dof = 3
    real(irealLUT) :: newcoeff(size(coeff))
    if (lswitch_east) then
      newcoeff = coeff
      call CHKERR(1_mpiint, 'not yet implemented')
    end if
    if (lswitch_north) then
      newcoeff = coeff
      call CHKERR(1_mpiint, 'not yet implemented')
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
      if (lswitch_east .or. lswitch_north) coeff = coeff
    end if
  end subroutine

  !for solver_8_10 the offset is chaning and the destination order
  subroutine dir8_to_diff10_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_10) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 8
    real(irealLUT) :: newcoeff(size(coeff))
    if (lswitch_east) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 2-1) )
      coeff(1 + (3 - 1) * dof:3 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * (4 - 1))
      coeff(1 + (4 - 1) * dof:4 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * (3 - 1))
      coeff(1 + (5 - 1) * dof:5 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * (6 - 1))
      coeff(1 + (6 - 1) * dof:6 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * (5 - 1))
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 7-1) )
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 8-1) )
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*( 9-1) )
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof*(10-1) )
    end if
    if (lswitch_north) then
      newcoeff = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 1-1) )
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 2-1) )
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 4-1) )
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 3-1) )
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 6-1) )
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof*( 5-1) )
      coeff(1 + (7 - 1) * dof:7 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * (8 - 1))
      coeff(1 + (8 - 1) * dof:8 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * (7 - 1))
      coeff(1 + (9 - 1) * dof:9 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * (10 - 1))
      coeff(1 + (10 - 1) * dof:10 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * (9 - 1))
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
    end if
  end subroutine

  subroutine dir8_to_diff16_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_16) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    real(irealLUT) :: newcoeff(size(coeff)) ! dim(src,dst)
    integer(iintegers), parameter :: dof = 8
    if (lswitch_east) then
      newcoeff(:) = coeff
      !coeff(1+( 1-1)*dof: 1*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 1-1)*dof)
      !coeff(1+( 2-1)*dof: 2*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 2-1)*dof)
      coeff(1 + (3 - 1) * dof:3 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (7 - 1) * dof)
      coeff(1 + (4 - 1) * dof:4 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (8 - 1) * dof)
      !coeff(1+( 5-1)*dof: 5*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 5-1)*dof)
      !coeff(1+( 6-1)*dof: 6*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+( 6-1)*dof)
      coeff(1 + (7 - 1) * dof:7 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (3 - 1) * dof)
      coeff(1 + (8 - 1) * dof:8 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (4 - 1) * dof)
      coeff(1 + (9 - 1) * dof:9 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (10 - 1) * dof)
      coeff(1 + (10 - 1) * dof:10 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (9 - 1) * dof)
      coeff(1 + (11 - 1) * dof:11 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (12 - 1) * dof)
      coeff(1 + (12 - 1) * dof:12 * dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + (11 - 1) * dof)
      !coeff(1+(13-1)*dof:13*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(13-1)*dof)
      !coeff(1+(14-1)*dof:14*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(14-1)*dof)
      !coeff(1+(15-1)*dof:15*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(15-1)*dof)
      !coeff(1+(16-1)*dof:16*dof) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+(16-1)*dof)
    end if
    if (lswitch_north) then
      newcoeff(:) = coeff
      coeff(1 + (1 - 1) * dof:1 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (5 - 1) * dof)
      coeff(1 + (2 - 1) * dof:2 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (6 - 1) * dof)
      !coeff(1+( 3-1)*dof: 3*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 3-1)*dof)
      !coeff(1+( 4-1)*dof: 4*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 4-1)*dof)
      coeff(1 + (5 - 1) * dof:5 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (1 - 1) * dof)
      coeff(1 + (6 - 1) * dof:6 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (2 - 1) * dof)
      !coeff(1+( 7-1)*dof: 7*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 7-1)*dof)
      !coeff(1+( 8-1)*dof: 8*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 8-1)*dof)
      !coeff(1+( 9-1)*dof: 9*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+( 9-1)*dof)
      !coeff(1+(10-1)*dof:10*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(10-1)*dof)
      !coeff(1+(11-1)*dof:11*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(11-1)*dof)
      !coeff(1+(12-1)*dof:12*dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+(12-1)*dof)
      coeff(1 + (13 - 1) * dof:13 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (14 - 1) * dof)
      coeff(1 + (14 - 1) * dof:14 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (13 - 1) * dof)
      coeff(1 + (15 - 1) * dof:15 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (16 - 1) * dof)
      coeff(1 + (16 - 1) * dof:16 * dof) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + (15 - 1) * dof)
    end if
    if (lswitch_east) then
    end if
    if (lswitch_north) then
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
    end if
  end subroutine

  subroutine dir8_to_diff18_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_8_18) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), target, intent(inout) :: coeff(:)
    if (lswitch_east) then
      call CHKERR(1_mpiint, 'dir8_to_diff18_coeff_symmetry_lswitch_east_not yet implemented')
      coeff = coeff
    end if
    if (lswitch_north) then
      call CHKERR(1_mpiint, 'dir8_to_diff18_coeff_symmetry_lswitch_north_not yet implemented')
      coeff = coeff
    end if
    select type (OPP)
    end select
  end subroutine

  subroutine dir2dir_coeff_symmetry_none(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), intent(inout) :: coeff(:)
    return
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
      if (lswitch_east .or. lswitch_north) coeff = coeff
    end if
  end subroutine

  subroutine dir2dir8_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop_cube_dir8) :: OPP
    logical, intent(in) :: lswitch_east, lswitch_north
    real(irealLUT), intent(inout) :: coeff(:)
    integer(iintegers), parameter :: dof = 8
    real(irealLUT) :: newcoeff(size(coeff))

    if (lswitch_east) then
      newcoeff = coeff
      coeff(1:8) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof)
      coeff(9:16) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8])
      coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 3)
      coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 2)
      coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 4)
      coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 5)
      coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 6)
      coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] + dof * 7)
    end if
    if (lswitch_north) then
      newcoeff = coeff
      coeff(1:8) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 2)
      coeff(9:16) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 3)
      coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8])
      coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof)
      coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 4)
      coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 5)
      coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 6)
      coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] + dof * 7)
    end if
    if (.false.) then ! remove compiler unused warnings
      select type (OPP)
      end select
    end if
  end subroutine
end module
