module test_plexrt_fish

  use m_examples_plex_ex_fish, only: plex_ex_fish
!#include "petsc/finclude/petsc.h"
!use petsc
!
use m_tenstream_options, only : read_commandline_options

use m_data_parameters, only: &
  & init_mpi_data_parameters, &
  & finalize_mpi, &
  & iintegers, ireals, mpiint, &
  & zero, one

use m_helper_functions, only: &
  & CHKERR, &
  & normalize_vec, &
  & spherical_2_cartesian

!use m_helper_functions, only: chkerr, cstr, itoa, &
!  & approx, meanval, imp_allreduce_mean, meanvec, reverse

!use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline
!use m_plex_grid, only: t_plexgrid, setup_plexgrid
!use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D
!
!use m_plex_rt, only: &
!  init_plex_rt_solver
!
!use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg
!
!use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, print_tenstr_atm, &
!  & hydrostat_plev

use pfunit_mod

implicit none

contains

  @before
  subroutine setup(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm     = this%getMpiCommunicator()
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm     = this%getMpiCommunicator()
    call finalize_mpi(comm, lfinalize_mpi=.False., lfinalize_petsc=.True.)
  end subroutine teardown

  @test(npes =[1,2,4])
  subroutine plexrt_regular_sw_sza0(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=4, Ny=5, Nz=3
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh=.True.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0.,0.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz,:)+edn(Nz,:))*Ag, eup(Nz,:), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes =[1,2,4])
  subroutine plexrt_regular_sw_sza20(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=6, Ny=3, Nz=3
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=2e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh=.True.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90.,20.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz,:)+edn(Nz,:))*Ag, eup(Nz,:), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes =[1,2,4])
  subroutine plexrt_regular_sw_sza40_no_scatter_no_Ag(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=6, Ny=3, Nz=4
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=5e-2

    logical, parameter :: lregular_mesh=.True.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90.,40.)
    w0 = 0
    Ag = 0

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(zero, eup, eps, 'no scatter, no albedo, eup should be 0')
    @assertEqual(zero, edn, eps, 'no scatter, no albedo, edn should be 0')
  end subroutine

  @test(npes =[1,2,4])
  subroutine plexrt_fish_sw_sza0(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=4, Ny=5, Nz=3
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=1e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh=.False.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(0.,0.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz,:)+edn(Nz,:))*Ag, eup(Nz,:), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes =[1,2,4])
  subroutine plexrt_fish_sw_sza20(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=6, Ny=3, Nz=3
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=2e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh=.False.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90.,20.)
    w0 = .5
    Ag = .3

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual((edir(Nz,:)+edn(Nz,:))*Ag, eup(Nz,:), eps, 'eup should be (edir+edn)*Ag')
  end subroutine

  @test(npes =[1,2,4])
  subroutine plexrt_fish_sw_sza40_no_scatter_no_Ag(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    integer(iintegers), parameter :: Nx=6, Ny=3, Nz=3
    logical, parameter :: lthermal=.False., lsolar=.True., lverbose=.True.
    real(ireals), parameter :: dz=1, dtau=1, Bplck = -1
    real(ireals), parameter :: eps=5e-2 ! better than 1 percent

    logical, parameter :: lregular_mesh=.False.
    real(ireals) :: w0, Ag
    real(ireals) :: sundir(3)
    real(ireals), allocatable, dimension(:,:) :: edir, edn, eup, abso

    comm     = this%getMpiCommunicator()

    sundir(:) = spherical_2_cartesian(90.,40.)
    w0 = 0
    Ag = 0

    call plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    @assertEqual(-sundir(3), edir(1,:), eps, 'TOA direct radiation should only depend on sundir')
    @assertEqual(-sundir(3)*exp(-one/abs(sundir(3))), edir(Nz,:), eps, 'Srfc direct radiation should be TOA direct incoming with LambertBeer')
    @assertEqual(zero, eup, eps, 'no scatter, no albedo, eup should be 0')
    @assertEqual(zero, edn, eps, 'no scatter, no albedo, edn should be 0')
  end subroutine
end module
