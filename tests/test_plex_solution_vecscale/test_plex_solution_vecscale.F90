module test_plex_solution_vecscale

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only: CHKERR, approx
  use m_data_parameters, only: ireals, iintegers, mpiint, &
                               i0, i1, i2, i3, i4, i5, &
                               zero, one, &
                               init_mpi_data_parameters

  use m_icon_grid, only: t_icongrid, read_icon_grid_file, &
                         bcast_icongrid, distribute_icon_grid

  use m_plex_grid, only: t_plexgrid, &
                         setup_edir_dmplex, setup_abso_dmplex, compute_face_geometry, &
                         ncvar2d_to_globalvec, setup_plexgrid, &
                         gen_test_mat, get_normal_of_first_toa_face

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: compute_face_geometry, &
                       init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
                       plexrt_get_result, destroy_plexrt_solver, scale_flx

  use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D

  use m_pprts_base, only: t_state_container, prepare_solution, destroy_solution

  use pfunit_mod
  implicit none

contains

  @test(npes=[1])
  subroutine plex_ex3(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: Nx = 2, Ny = 3, Nz = 2
    real(ireals), parameter :: dz = 1.5_ireals

    type(tDM) :: dm2d_serial, dm2d, dm3d
    real(ireals) :: hhl(Nz)

    integer(mpiint) :: myid, numnodes, comm, ierr
    integer(iintegers) :: k
    type(tVec), allocatable :: dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2
    type(tVec), allocatable :: diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2
    type(t_state_container) :: solution

    class(t_plex_solver), allocatable :: solver
    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)

    call create_2d_fish_plex(comm, Nx, Ny, dm2d_serial, dm2d)
    call DMDestroy(dm2d_serial, ierr); call CHKERR(ierr)

    hhl(1) = zero
    do k = 2, Nz
      hhl(k) = hhl(k - 1) - dz
    end do

    call dmplex_2D_to_3D(dm2d, Nz, hhl, [zero, zero, -huge(one)], dm3d, zindex)

    call setup_plexgrid(dm2d, dm3d, Nz - 1, zindex, plex, hhl)
    deallocate (zindex)

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)

    call prepare_solution(solver%plex%edir_dm, solver%plex%ediff_dm, solver%plex%abso_dm, &
      & lsolar=.true., lthermal=.false., solution=solution)

    print *, 'Testing Scalevec Direct'
    call init_and_scalevecs(solution, one, solution%edir, solution%lWm2_dir)
    print *, 'Testing Scalevec Diffuse'
    call init_and_scalevecs(solution, one, solution%ediff, solution%lWm2_diff)

    call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-show_solution_edir_vec1', ierr); call CHKERR(ierr)

  contains
    subroutine init_and_scalevecs(solution, initialvar, solution_vec, solution_lWm2)
      type(t_state_container), intent(inout) :: solution
      real(ireals), intent(in) :: initialvar
      type(tVec), intent(inout) :: solution_vec
      logical, intent(inout) :: solution_lWm2
      type(tVec) :: tmp_vec
      real(ireals), pointer :: xa(:), xb(:)
      integer(iintegers) :: i
      real(ireals), parameter :: eps = 10 * epsilon(eps)
      real(ireals), parameter :: dx = 1, dy = 1, h = sqrt(dy**2 - (dx / 2)**2)

      solution_lWm2 = .true.
      call VecDuplicate(solution_vec, tmp_vec, ierr); call CHKERR(ierr)
      call VecSet(solution_vec, initialvar, ierr); call CHKERR(ierr)
      call VecSet(tmp_vec, initialvar, ierr); call CHKERR(ierr)

      call scale_flx(solver, solver%plex, &
                     dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, &
                     diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, &
                     solution, lWm2=.true.)

      call VecGetArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      do i = 1, size(xa)
        @assertEqual(xa(i), xb(i), eps, 'Vec should not have changed after scaling from Wm-2 to Wm-2')
      end do
      call VecRestoreArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)

      call scale_flx(solver, solver%plex, &
                     dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, &
                     diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, &
                     solution, lWm2=.false.)

      call VecGetArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      do i = 1, size(xa)
        @assertFalse(approx(xa(i), xb(i), eps), 'Vec should have changed after scaling from Wm-2 to W')
      end do
      ! Some hardcoded values:
      if (size(xa) .eq. 52) then ! diffuse fish
        do i = 1, 16
          @assertEqual(xb(i) * dx * h / 2, xa(i), eps, 'Should match the top/bot area of size 52 fish')
        end do
        do i = 17, 52
          @assertEqual(xb(i) * dz * sqrt((dx / 2)**2 + h**2), xa(i), eps, 'Should match the side area')
        end do
      elseif (size(xa) .eq. 17) then
        do i = 1, 8
          @assertEqual(xb(i) * dx * h / 2, xa(i), eps, 'Should match the top/bot area of size 17 fish')
        end do
        do i = 9, 17
          @assertEqual(xb(i) * dz * sqrt((dx / 2)**2 + h**2), xa(i), eps, 'Should match the side area')
        end do
      end if

      call VecRestoreArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)

      call scale_flx(solver, solver%plex, &
                     dir_scalevec_Wm2_to_W, dir_scalevec_W_to_Wm2, &
                     diff_scalevec_Wm2_to_W, diff_scalevec_W_to_Wm2, &
                     solution, lWm2=.true.)

      call VecGetArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)
      call VecGetArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      do i = 1, size(xa)
        @assertEqual(xa(i), xb(i), eps, 'Vec should not have changed after scaling back and forth from Wm-2 to Wm-2')
      end do
      call VecRestoreArrayReadF90(tmp_vec, xb, ierr); call CHKERR(ierr)
      call VecRestoreArrayReadF90(solution_vec, xa, ierr); call CHKERR(ierr)
    end subroutine
  end subroutine

end module
