module test_fish_plex
#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only: &
    mpiint, ireals, iintegers, &
    one, zero, default_str_len, &
    i0, i1, i2, i3, i4, &
    init_mpi_data_parameters

  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: CHKERR

  use m_icon_plex_utils, only: create_2d_fish_plex
  use m_plex_grid, only: print_dmplex, dmplex_set_new_section

  use pfunit_mod
  implicit none

  integer(mpiint) :: comm, myid, numnodes, ierr

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call PetscFinalize(ierr)
    if (myid .eq. 0) print *, 'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes=[1])
  subroutine test_connec_serial_fish_Nx2_Ny3(this)
    class(MpiTestMethod), intent(inout) :: this
    type(tDM) :: dm, dmdist
    integer(iintegers), parameter :: Nx = 2, Ny = 3
    integer(iintegers) :: pStart, pEnd, fStart, fEnd, eStart, eEnd, vStart, vEnd
    integer(iintegers) :: target_closure(14)

    call create_2d_fish_plex(comm, Nx, Ny, dm, dmdist)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
    call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    @assertEqual(0_iintegers, pStart)
    @assertEqual(19_iintegers, pEnd)
    @assertEqual(0_iintegers, fStart)
    @assertEqual(4_iintegers, fEnd)
    @assertEqual(4_iintegers, eStart)
    @assertEqual(13_iintegers, eEnd)
    @assertEqual(13_iintegers, vStart)
    @assertEqual(19_iintegers, vEnd)

    target_closure(:) = i0
    target_closure(1:14:2) = [0, 4, 5, 6, 13, 14, 15]
    call check_transclosure(dm, i0, target_closure)

    target_closure(1:14:2) = [1, 6, 7, 8, 14, 15, 16]
    call check_transclosure(dm, i1, target_closure)

    target_closure(1:14:2) = [2, 9, 10, 12, 15, 17, 18]
    call check_transclosure(dm, i2, target_closure)

    target_closure(1:14:2) = [3, 8, 10, 11, 15, 16, 18]
    call check_transclosure(dm, i3, target_closure)

    call DMDestroy(dm, ierr)
    call DMDestroy(dmdist, ierr)
  end subroutine

  @test(npes=[1])
  subroutine test_connec_serial_fish_Nx4_Ny5(this)
    class(MpiTestMethod), intent(inout) :: this
    type(tDM) :: dm, dmdist
    integer(iintegers), parameter :: Nx = 4, Ny = 5
    integer(iintegers) :: pStart, pEnd, fStart, fEnd, eStart, eEnd, vStart, vEnd
    integer(iintegers) :: target_closure(14)

    call create_2d_fish_plex(comm, Nx, Ny, dm, dmdist)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
    call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
    call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

    @assertEqual(0_iintegers, pStart)
    @assertEqual(61_iintegers, pEnd)
    @assertEqual(0_iintegers, fStart)
    @assertEqual(16_iintegers, fEnd)
    @assertEqual(16_iintegers, eStart)
    @assertEqual(46_iintegers, eEnd)
    @assertEqual(46_iintegers, vStart)
    @assertEqual(61_iintegers, vEnd)

    target_closure(:) = i0
    target_closure(1:14:2) = [0, 16, 18, 19, 46, 47, 49]
    call check_transclosure(dm, target_closure(1), target_closure)

    target_closure(1:14:2) = [1, 19, 20, 23, 47, 49, 50]
    call check_transclosure(dm, target_closure(1), target_closure)

    target_closure(1:14:2) = [2, 17, 20, 21, 47, 48, 50]
    call check_transclosure(dm, target_closure(1), target_closure)

    target_closure(1:14:2) = [11, 35, 36, 38, 54, 56, 57]
    call check_transclosure(dm, target_closure(1), target_closure)

    target_closure(1:14:2) = [12, 39, 40, 44, 55, 58, 59]
    call check_transclosure(dm, target_closure(1), target_closure)

    target_closure(1:14:2) = [15, 38, 42, 43, 56, 57, 60]
    call check_transclosure(dm, target_closure(1), target_closure)

    call DMDestroy(dm, ierr)
    call DMDestroy(dmdist, ierr)
  end subroutine

  subroutine check_transclosure(dm, icell, target_closure)
    type(tDM), intent(in) :: dm
    integer(iintegers), intent(in) :: icell, target_closure(:)
    integer(iintegers) :: i
    PetscCount :: transclosure_size
    integer(iintegers), pointer :: transclosure(:)

    call DMPlexGetTransitiveClosure(dm, icell, PETSC_TRUE, PETSC_NULL_INTEGER, transclosure, ierr); call CHKERR(ierr)
    transclosure_size = size(transclosure(1:size(transclosure):2))
    call PetscSortInt(transclosure_size, &
                      transclosure(1:size(transclosure):2), ierr); call CHKERR(ierr)

    if (size(transclosure) .ne. size(target_closure)) then
      print *, 'cell transclosure ', icell, ':', transclosure(1:size(transclosure):2)
      print *, 'cell targetclosure', icell, ':', target_closure(1:size(target_closure):2)
    end if
    if (.not. all(transclosure .eq. target_closure)) then
      print *, 'cell transclosure ', icell, ':', transclosure(1:size(transclosure):2)
      print *, 'cell targetclosure', icell, ':', target_closure(1:size(target_closure):2)
    end if
    do i = 1, size(transclosure)
      @assertEqual(target_closure(i), transclosure(i))
    end do
    call DMPlexRestoreTransitiveClosure(dm, i1, PETSC_TRUE, PETSC_NULL_INTEGER, transclosure, ierr); call CHKERR(ierr)
  end subroutine
end module
