module test_plexrt_nca
#include "petsc/finclude/petsc.h"
use petsc

use m_tenstream_options, only : read_commandline_options
use m_data_parameters, only: init_mpi_data_parameters, &
  iintegers, ireals, mpiint, default_str_len, &
  i0
use m_helper_functions, only: triangle_area_by_edgelengths, chkerr, meanval, itoa

use m_plexrt_nca, only: plexrt_nca_init, plexrt_nca

use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline
use m_plex_grid, only: t_plexgrid, setup_plexgrid
use m_icon_plex_utils, only: create_2d_fish_plex, dmplex_2D_to_3D, dump_ownership

use m_plex_rt, only: &
  init_plex_rt_solver, &
  set_plex_rt_optprop, &
  run_plex_rt_solver, &
  plexrt_get_result, &
  destroy_plexrt_solver

use pfunit_mod

implicit none

contains
  @test(npes =[1])
  subroutine test_single_box_nca(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    real(ireals), parameter :: dx1=1e3, dx2=dx1, dx3=dx1, dz=1e3
    real(ireals) :: atop, abot, a1, a2, a3, vol, hr
    real(ireals) :: base_info(7), side_info(3*5)
    real(ireals), parameter :: kabs=1e-5

    comm     = this%getMpiCommunicator()

    atop = triangle_area_by_edgelengths(dx1, dx2, dx3)
    abot = atop
    a1   = dx1 * dz
    a2   = dx2 * dz
    a3   = dx3 * dz

    vol = atop * dz

    base_info = [        &
      kabs,        & !kabs
      kabs,        & !kabs_top
      0._ireals,   & !Ldn_top
      1._ireals,   & !Btop
      kabs,        & !kabs_bot
      0._ireals,   & !Lup_bot
      1._ireals    & !Bbot
      ]
    side_info(1:5) = [   &
      kabs,        & !kabs
      0._ireals,   & !Ldn_top
      0._ireals,   & !Lup_top
      0._ireals,   & !Ldn_bot
      0._ireals    & !Lup_bot
      ]

    side_info( 6:10) = side_info(1:5)
    side_info(11:15) = side_info(1:5)

    call plexrt_nca_init(comm)

    call plexrt_nca (dx1, dx2, dx3, &
      dz, atop, abot, a1, a2, a3, vol, &
      base_info, side_info, hr)

    @assertEqual(-9.9e-5_ireals, hr, 1e-6_ireals, 'NCA Heating Test 1 false')
  end subroutine

  @test(npes =[1])
  subroutine test_horizontally_homogeneous_layer(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    real(ireals), parameter :: dx1=1e3, dx2=dx1, dx3=dx1, dz=1e3
    real(ireals) :: atop, abot, a1, a2, a3, vol, hr
    real(ireals) :: base_info(7), side_info(3*5)
    real(ireals), parameter :: kabs=10

    comm     = this%getMpiCommunicator()

    atop = triangle_area_by_edgelengths(dx1, dx2, dx3)
    abot = atop
    a1   = dx1 * dz
    a2   = dx2 * dz
    a3   = dx3 * dz

    vol = atop * dz

    base_info = [        &
      kabs,        & !kabs
      0._ireals,   & !kabs_top
      0._ireals,   & !Ldn_top
      1._ireals,   & !Btop
      0._ireals,   & !kabs_bot
      0._ireals,   & !Lup_bot
      1._ireals    & !Bbot
      ]
    side_info(1:5) = [   &
      kabs,        & !kabs
      0._ireals,   & !Ldn_top
     10._ireals,   & !Lup_top
     10._ireals,   & !Ldn_bot
      0._ireals    & !Lup_bot
      ]

    side_info( 6:10) = side_info(1:5)
    side_info(11:15) = side_info(1:5)

    call plexrt_nca_init(comm)

    call plexrt_nca (dx1, dx2, dx3, &
      dz, atop, abot, a1, a2, a3, vol, &
      base_info, side_info, hr)

    @assertEqual(0._ireals, hr, 1e-6_ireals)
  end subroutine

  @test(npes =[1])
  subroutine test_nca_dmplex_interface(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: myid, numnodes, comm, ierr

    integer(iintegers), parameter :: Nx=2, Ny=3, Nz=4
    real(ireals), parameter :: dz=1_ireals, Ag=1
    real(ireals), parameter :: sundir(3)=[0,0,1]

    type(tDM) :: dm2d_serial, dm2d, dm3d
    real(ireals) :: hhl(Nz)

    class(t_plex_solver), allocatable :: solver
    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    real(ireals), allocatable, dimension(:,:) :: edn, eup, abso
    real(ireals), pointer :: xv(:), xxv(:,:)
    integer(iintegers) :: k, i, fStart, fEnd, Ncol

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)

    call create_2d_fish_plex(comm, Nx, Ny, dm2d_serial, dm2d)
    call DMDestroy(dm2d_serial, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm2d, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart

    hhl(1) = 0
    do k=2,Nz
    hhl(k) = hhl(k-1) - dz
    enddo

    call dmplex_2D_to_3D(dm2d, Nz, hhl, dm3d, zindex, lpolar_coords=.False.)
    call dump_ownership(dm3d, '-dump_ownership', '-show_plex')

    call setup_plexgrid(dm3d, Nz-1, zindex, plex)
    deallocate(zindex)

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)
    call set_plex_rt_optprop(solver, vert_integrated_kabs=0._ireals, vert_integrated_ksca=0._ireals)

    call VecGetArrayF90(solver%kabs, xv, ierr); call CHKERR(ierr)
    xxv(1:Nz-1, 1:Ncol) => xv
    if(myid.eq.0) xxv(Nz/2,:) = 10
    nullify(xxv)
    call VecRestoreArrayF90(solver%kabs, xv, ierr); call CHKERR(ierr)

    if(.not.allocated(solver%albedo)) then
      allocate(solver%albedo)
      call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
    endif
    call VecSet(solver%albedo, Ag, ierr); call CHKERR(ierr)

    if(.not.allocated(solver%plck)) then
      allocate(solver%plck)
      call DMCreateGlobalVector(solver%plex%horizface1_dm, solver%plck, ierr); call CHKERR(ierr)
    endif
    call VecSet(solver%plck, 100._ireals, ierr); call CHKERR(ierr)

    call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=sundir)

    call plexrt_get_result(solver, edn, eup, abso)

    if(allocated(abso)) then
      do i=0,numnodes-1
      if(myid.eq.i) then
        print *, ''
        print *, 'Averages on rank'//itoa(myid)
        do k = 1, ubound(abso,1)
        print *, k, meanval(edn(k,:)), meanval(eup(k,:)), meanval(abso(k,:))
        enddo
        print *, k, meanval(edn(k,:)), meanval(eup(k,:))
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
      enddo
    endif

    call destroy_plexrt_solver(solver, lfinalizepetsc=.True.)
  end subroutine

end module
