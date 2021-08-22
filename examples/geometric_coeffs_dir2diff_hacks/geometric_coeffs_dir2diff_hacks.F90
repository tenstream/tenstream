module m_example_geometric_coeffs_dir2diff_hacks

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, &
    irealLUT, ireal_dp
  use m_helper_functions, only : write_ascii_file_2d, toStr, deg2rad, spherical_2_cartesian
  use m_geometric_coeffs, only: dir2dir3_geometric_coeffs
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_boxmc, only : t_boxmc_3_10
  use m_optprop, only : t_optprop_3_10
  use m_tenstream_options, only: read_commandline_options

  implicit none
  type(t_boxmc_3_10) :: bmc_3_10
  type(t_optprop_3_10) :: OPP
  integer(mpiint) :: myid,mpierr,numnodes,comm

contains

  subroutine ex_geometric_coeffs_dir2diff_hacks(comm)
    integer(mpiint), intent(in) :: comm
    real(ireals) :: verts_dtd(24), bg(3), phi, theta, sundir(3), &
      c_scatter_gomtrc_reg(9), c_scatter_gomtrc_dst(9), c_gomtrc_dst(9)
    real(ireals), parameter :: dx=1, dy=dx, dz=dx
    real(ireals), allocatable :: verts(:)
    integer(iintegers), parameter :: imax=5
    integer(iintegers) :: i, src, iscat, iabso
    real(ireals) :: S(10), T(3)
    real(irealLUT) :: S_LUT(30), T_LUT(9)
    real(ireals) :: S_tol(10), T_tol(3)
    real(ireal_dp), parameter :: atol=1e-3, rtol=1e-2
    real(ireals) :: S_bmc(30), S_bmc_reg_cscat_cabso(imax-1,30), S_LUT_reg_cscat_cabso(imax-1,30), &
      T_gomtrc_reg_cscat(imax-1,9), S_bmc_dst_cscat_cabso(imax-1,30), T_gomtrc_dst_cscat(imax-1,9), &
      T_gomtrc_dst(imax-1,9), T_LUT_reg_cscat_cabso(imax-1,9) ! alpha_x, alpha_y, alpha_xy
    real(irealLUT) :: aspect_zx, tauz, w0, g
    integer(mpiint) :: ierr
    logical :: xinc, yinc

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
    call bmc_3_10%init(comm)

    call OPP%init(comm)

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    phi = 0
    theta = 30

    if (sin(deg2rad(phi)) .gt. zero) then
      xinc = .False.
    else
      xinc = .True.
    endif

    if (cos(deg2rad(phi)) .lt. zero) then
      yinc = .False.
    else
      yinc = .True.
    endif

    do iscat = 0, imax - 2
      do iabso = 0, imax - 2
        sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]
        !bg = [real(ireals) :: 0._ireals, 1._ireals, 0.85]
        bg = [real(ireals) :: exp(real(iabso, ireals)) - 1, exp(real(iscat, ireals)), 0.85]

        call dir2dir3_geometric_coeffs(verts, sundir, bg(2), c_scatter_gomtrc_reg)

        tauz = real(bg(1) + bg(2), irealLUT)
        w0 = real(min(bg(2) / max(bg(1) + bg(2), epsilon(bg)), 0.999989986_ireals), irealLUT)
        g = real(bg(3), irealLUT)
        aspect_zx = real(dz / dx, irealLUT)

        call OPP%get_coeff(&
          tauz, & !tauz
          w0, & !w0
          g, & ! g
          aspect_zx, & !aspect_zx
          .False., & !ldir
          S_LUT, & !coeff
          ierr, & !ierr
          [real(phi, irealLUT), real(theta, irealLUT)], & !angles
          xinc, & !lswitch_east (180-360 False)
          yinc & !lswitch_north (90-270 False)
          )
        call OPP%get_coeff(&
          tauz, & !tauz
          w0, & !w0
          g, & ! g
          aspect_zx, & !aspect_zx
          .True., & !ldir
          T_LUT, & !coeff
          ierr, & !ierr
          [real(phi, irealLUT), real(theta, irealLUT)], & !angles
          xinc, & !lswitch_east (180-360 False)
          yinc & !lswitch_north (90-270 False)
          )

        do src = 1,3
          call bmc_3_10%get_coeff(&
            comm,real(bg, ireal_dp),src,.True.,real(phi, ireal_dp),real(theta, ireal_dp),&
            real(verts, ireal_dp),S,T,S_tol,T_tol,inp_atol=atol,inp_rtol=rtol)
          S_bmc(src:size(S_bmc):3) = S
        enddo

        do i = 1, imax - 1
          T_gomtrc_reg_cscat(i,:) = c_scatter_gomtrc_reg
          S_LUT_reg_cscat_cabso(i,:) = S_LUT
          T_LUT_reg_cscat_cabso(i,:) = T_LUT
          S_bmc_reg_cscat_cabso(i,:) = S_bmc
        enddo

        do i = 2, imax
          verts_dtd = verts
          verts_dtd( 3) = verts( 3) + 2 * dz / real(i, ireals)
          verts_dtd( 6) = verts( 6) + dz / real(i, ireals)
          verts_dtd( 9) = verts( 9) + dz / (2 * real(i, ireals))
          verts_dtd(12) = verts( 6) + verts_dtd(9) - verts_dtd(3)
          verts_dtd(15:24:3) = verts_dtd(3:12:3) + dz

          call dir2dir3_geometric_coeffs(&
            verts_dtd, sundir, bg(2), c_scatter_gomtrc_dst)

          call dir2dir3_geometric_coeffs(&
            verts_dtd, sundir, bg(1) + bg(2), c_gomtrc_dst)

          do src = 1,3
            call bmc_3_10%get_coeff(&
              comm,real(bg, ireal_dp),src,.True.,real(phi, ireal_dp),real(theta, ireal_dp),real(verts_dtd, ireal_dp), &
              S,T,S_tol,T_tol,inp_atol=atol,inp_rtol=rtol)
            S_bmc(src:size(S_bmc):3) = S
          enddo

          T_gomtrc_dst_cscat(i-1,:) = c_scatter_gomtrc_dst
          T_gomtrc_dst(i-1,:) = c_gomtrc_dst
          S_bmc_dst_cscat_cabso(i-1,:) = S_bmc
        enddo

        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'S_LUT_reg_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          S_LUT_reg_cscat_cabso, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))

        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'T_LUT_reg_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          T_LUT_reg_cscat_cabso, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))

        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'S_bmc_reg_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          S_bmc_reg_cscat_cabso, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))
        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'S_bmc_dst_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          S_bmc_dst_cscat_cabso, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))
        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'gomtrc_reg_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          T_gomtrc_reg_cscat, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))
        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'gomtrc_dst_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          T_gomtrc_dst_cscat, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))
        call write_ascii_file_2d( &
          '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/'//&
          'tenstream_geometric_coeffs/data_ex_geometric_coeffs_dir2diff_hacks/'//&
          'T_gomtrc_dst_phi_'//sStr(phi)//'_theta_'//sStr(theta)//'_cscatter_'//sStr(bg(2))//'_cabso_'//sStr(bg(1))//'.out', &
          T_gomtrc_dst, ierr, header='# src_z src_x src_y ; cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2)))
      enddo
    enddo

    contains
      function sStr(num)
        real(ireals), intent(in) :: num
        character(:), allocatable :: sstr
        sstr = trim(adjustl(toStr(num)))
      end function
  end subroutine
end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : mpiint
  use m_helper_functions, only : CHKERR
  use m_example_geometric_coeffs_dir2diff_hacks, only : ex_geometric_coeffs_dir2diff_hacks

  implicit none

  integer(mpiint) :: ierr, myid

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER ,ierr); call CHKERR(ierr)

  call ex_geometric_coeffs_dir2diff_hacks(mpi_comm_world)

end program
