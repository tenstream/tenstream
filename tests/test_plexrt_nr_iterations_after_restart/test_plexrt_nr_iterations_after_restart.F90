module test_plexrt_nr_iterations_after_restart

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, one, zero, i0

  use m_helper_functions, only: CHKERR, reverse, itoa

  use m_adaptive_spectral_integration, only: need_new_solution

  use m_plex_grid, only: t_plexgrid, setup_plexgrid
  use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline
  use m_plex_rt, only: &
    init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
    plexrt_get_result, destroy_plexrt_solver

  use pfunit_mod

  implicit none

  class(t_plex_solver), allocatable :: solver

contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
  end subroutine teardown

  @test(npes =[2,1])
  subroutine error_growth_tracking(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: iter, k, vStart, vEnd
    integer(mpiint) :: comm, myid, numnodes, ierr

    integer(iintegers),parameter :: Nx=20,Ny=25,Nz=10
    real(ireals),parameter :: dx=100
    real(ireals),parameter :: sundir(3) = [zero, zero, -one]
    real(ireals),parameter :: albedo=0, dz=dx

    real(ireals),allocatable,dimension(:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: time
    logical :: lneed

    type(tDM) :: dm_serial, dm2d, dm3d
    real(ireals) :: hhl(Nz)
    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)

    call create_2d_fish_plex(comm, Nx, Ny, dm_serial, dm2d, opt_dx=dx)
    call DMDestroy(dm_serial, ierr); call CHKERR(ierr)

    hhl(1) = zero
    do k=2,Nz
      hhl(k) = hhl(k-1) + dz
    enddo
    hhl = reverse(hhl)

    call DMPlexGetDepthStratum(dm2d, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! 2D vertices

    call dmplex_2D_to_3D(dm2d, Nz, hhl, dm3d, zindex, lpolar_coords=.False.)

    call setup_plexgrid(dm2d,dm3d, Nz-1, zindex, plex, hhl)
    deallocate(zindex)

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)
    deallocate(plex)

    call set_plex_rt_optprop(solver, vert_integrated_kabs=one, vert_integrated_ksca=.5_ireals)

    allocate(solver%albedo)
    call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
    call VecSet(solver%albedo, albedo, ierr); call CHKERR(ierr)

    ! Solar test
    k = 1
    call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., sundir=sundir, &
      opt_solution_uid=k)
    print *,'First KSP took '//itoa(solver%solutions(k)%Niter_dir)//' solar iterations'// &
      ' ... and '//itoa(solver%solutions(k)%Niter_diff)//' solar diffuse iterations'

    do iter=1,3
        call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., sundir=sundir, &
          opt_solution_uid=k)

        print *,'KSP took '//itoa(solver%solutions(k)%Niter_dir)//' solar iterations'// &
                ' ... and '//itoa(solver%solutions(k)%Niter_diff)//' solar diffuse iterations on call nr'//itoa(iter)
        call assertTrue(solver%solutions(k)%Niter_dir .eq.0)
        call assertTrue(solver%solutions(k)%Niter_diff.eq.0)
    enddo


    ! Thermal test
    if(.not.allocated(solver%plck)) then
      allocate(solver%plck)
      call DMCreateGlobalVector(solver%plex%horizface1_dm, solver%plck, ierr); call CHKERR(ierr)
    endif
    call VecSet(solver%plck, 100._ireals, ierr); call CHKERR(ierr)

    k = 2
    call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=sundir, &
      opt_solution_uid=k)
    print *,'First KSP took '//itoa(solver%solutions(k)%Niter_dir)//' thermal iterations'// &
      ' ... and '//itoa(solver%solutions(k)%Niter_diff)//' thermal diffuse iterations'

    do iter=1,3
        call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=sundir, &
          opt_solution_uid=k)

        print *,'KSP took '//itoa(solver%solutions(k)%Niter_dir)//' thermal iterations'// &
                ' ... and '//itoa(solver%solutions(k)%Niter_diff)//' thermal diffuse iterations on call nr'//itoa(iter)
        call assertTrue(solver%solutions(k)%Niter_diff.eq.0)
    enddo


    call destroy_plexrt_solver(solver, lfinalizepetsc=.True.)
  end subroutine

end module
