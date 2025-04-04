! Example computes geometric_coeffs coefficients using differen optional num_intervals (intervals for numerical iteration of
! triangle areas with extinction. Output is written into asci files to evaluate in python

module m_example_geometric_coeffs_num_intervals

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, one, ireal_dp

  use m_helper_functions, only: CHKERR, spherical_2_cartesian, cstr, write_ascii_file_2d, toStr

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_geometric_coeffs, only: dir2dir3_geometric_coeffs
  use m_boxmc_geometry, only: setup_default_unit_cube_geometry
  use m_boxmc, only: t_boxmc_3_10

  implicit none
  type(t_boxmc_3_10) :: bmc_3_10
  integer(mpiint) :: myid, mpierr, numnodes, comm

contains

  subroutine ex_geometric_coeffs_num_intervals(comm)
    integer(mpiint), intent(in) :: comm
    real(ireals) :: bg(3), phi, theta, sundir(3)
    real(ireals), parameter :: dx = 1, dy = dx, dz = dx
    real(ireals), allocatable :: verts(:)
    integer(iintegers), parameter :: imax = 100
    integer(iintegers) :: i, src
    real(ireals) :: S(10), T(3)
    real(ireals) :: S_tol(10), T_tol(3)
    real(ireal_dp), parameter :: atol = 1e-4, rtol = 1e-3
    real(ireals) :: c_gomtrc_reg_benchmark(9), c_gomtrc_reg(9, imax + 1), c_ext
    integer(mpiint) :: ierr

    call init_mpi_data_parameters(comm)
    call bmc_3_10%init(comm)

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    phi = 210
    theta = 40
    sundir = spherical_2_cartesian(phi, theta)*[-one, -one, one]
    bg = [real(ireals) :: 1._ireals, 0._ireals, 0.85]
    c_ext = bg(1) + bg(2)

    call dir2dir3_geometric_coeffs(verts, sundir, c_ext, c_gomtrc_reg_benchmark, num_intervals=10001_iintegers)
    do i = 1, imax
      call dir2dir3_geometric_coeffs(verts, sundir, c_ext, c_gomtrc_reg(:, i), num_intervals=i)
      if (i .eq. imax) print *, 'gomtrc', c_gomtrc_reg(:, i)
      c_gomtrc_reg(:, i) = (c_gomtrc_reg_benchmark - c_gomtrc_reg(:, i)) / &
                           (c_gomtrc_reg_benchmark + epsilon(c_gomtrc_reg))
    end do

    c_gomtrc_reg(:, imax + 1) = [real(ireals) :: &
                                 -2.2629924908925581e-003, -5.3639276943289577e-004, 1.6083701599837531e-003, &
                                 -2.3252373907683343e-003, 0.0000000000000000, -8.2549137394604885e-003, &
                                 -9.2645650761771661e-004, 1.4179496081441498e-003, 0.0000000000000000 &
                                 ]

    do src = 1, 3
      if (.false.) then
        call bmc_3_10%get_coeff( &
          comm, &
          real(bg, ireal_dp), &
          src, &
          .true., &
          real(phi, ireal_dp), &
          real(theta, ireal_dp), &
          real(verts, ireal_dp), &
          S, T, &
          S_tol, T_tol, &
          inp_atol=atol, &
          inp_rtol=rtol &
          )
        c_gomtrc_reg(src:9:3, imax + 1) = real(T, ireals)
      else
        print *, cstr('Montecarlo simulation not started, it is commented and the hard coded values are used.'//&
          &'If you changed the box geometry or bg, please make sure you uncomment the montecarlo call.', 'red')
      end if
    end do

    print *, 'bmc', c_gomtrc_reg(:, imax + 1)
    c_gomtrc_reg(:, imax + 1) = (c_gomtrc_reg_benchmark - c_gomtrc_reg(:, imax + 1)) / &
                                (c_gomtrc_reg_benchmark + epsilon(c_gomtrc_reg))

    call write_ascii_file_2d( &
      '/project/meteo/work/Hermann.Boettcher/hermannboettchermasterthesis/tenstream_geometric_coeffs/'// &
      'data_ex_geometric_coeffs_num_intervals/gomtrc_benchmark_reg_cabso_'//sStr(bg(1))//'.out', &
      c_gomtrc_reg, ierr, &
      header='# cabso='//sStr(bg(1))//' ; cscatter='//sStr(bg(2))//' ; c_bmc_reg as last entry')
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
  use m_data_parameters, only: mpiint
  use m_helper_functions, only: CHKERR
  use m_example_geometric_coeffs_num_intervals, only: ex_geometric_coeffs_num_intervals

  implicit none

  integer(mpiint) :: ierr, myid

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)

  call ex_geometric_coeffs_num_intervals(mpi_comm_world)

end program
