module test_icon_plex_utils

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    mpiint, ireals, iintegers, &
    one, zero, default_str_len, &
    i0, i1, i2, i3, i4, &
    init_mpi_data_parameters

  use m_helper_functions, only: CHKERR, deg2rad, itoa

  use m_plex_grid, only: create_plex_section

  use m_icon_plex_utils, only: date_to_julian_day, get_sun_vector, &
                               create_2d_regular_plex, rank0_f90vec_to_plex, plex_gVec_toZero, &
                               dmplex_gVec_from_f90_array

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
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    logical :: lpetsc_is_initialized
    call PetscInitialized(lpetsc_is_initialized, ierr)
    if (lpetsc_is_initialized) call PetscFinalize(ierr)
    if (myid .eq. 0) print *, 'Finishing icon_plex_utils tests module'
  end subroutine teardown

  @test(npes=[1])
  subroutine test_gregorian_date_to_julian_day(this)
    class(MpiTestMethod), intent(inout) :: this

    @assertEqual(2451545.0_ireals, date_to_julian_day(2000_iintegers, i1, 1.5_ireals), 'Greenwich Noon')

    @assertEqual(2446113.75_ireals, date_to_julian_day(1985_iintegers, i2, 17.25_ireals))
    @assertEqual(0._ireals, date_to_julian_day(-4712_iintegers, i1, 1.5_ireals))
    @assertEqual(2458365.29749_ireals, date_to_julian_day(2018_iintegers, 9_iintegers, 3._ireals+real(19*3600+8*60+23, ireals)/86400._ireals), 1e-4_ireals)

  end subroutine

  @test(npes=[1])
  subroutine test_sundirection_computation(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: year, month
    real(ireals) :: day

    real(ireals) :: sundir(3)
    integer :: hour, iday
    real(ireals) :: max_z_component
    max_z_component = sin(deg2rad(23.439_ireals))

    year = 2018
    month = 1

    do hour = 6, 18
      day = 5 + hour / 24._ireals
      sundir = get_sun_vector(year, month, day)

      @assertGreaterThan(sundir(1), zero)
      @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
    end do

    do hour = 1, 10
      day = 5 + hour / 24._ireals
      sundir = get_sun_vector(year, month, day)

      @assertGreaterThan(sundir(2), zero, sqrt(epsilon(sundir)))
      @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
    end do

    do hour = 13, 23
      day = 5 + hour / 24._ireals
      sundir = get_sun_vector(year, month, day)

      @assertLessThan(sundir(2), zero, sqrt(epsilon(sundir)))
      @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
    end do

    do iday = 1, 31
      do hour = 1, 24
        do month = 4, 8
          sundir = get_sun_vector(year, month, iday + (hour / 24._ireals))
          @assertGreaterThan(sundir(3), zero)
        end do
      end do

      do month = 1, 2
        sundir = get_sun_vector(year, month, iday + (hour / 24._ireals))
        @assertLessThan(sundir(3), zero)
      end do
      do month = 10, 12
        sundir = get_sun_vector(year, month, iday + (hour / 24._ireals))
        @assertLessThan(sundir(3), zero)
      end do
    end do

    do year = 2000, 2010, 1
      do iday = 1, 31
        do month = 1, 12
          do hour = 1, 24
            sundir = get_sun_vector(year, month, iday + (hour / 24._ireals))
            @assertLessThan(abs(sundir(3)), max_z_component, 'bad z_component '//itoa(year)//':'//itoa(month)//':'//itoa(iday))
          end do
        end do
      end do
    end do

    sundir = get_sun_vector(2018_iintegers, 9_iintegers, 3._ireals + real(19 * 3600 + 8 * 60 + 23, ireals) / 86400._ireals)
 @assertEqual([-0.29155451713720204_ireals, -0.94795805644649656_ireals, 0.12795111080046895_ireals], sundir, sqrt(epsilon(sundir)))
  end subroutine

  @test(npes=[1, 2, 3])
  subroutine test_rank0_to_distributed_plex(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    integer(iintegers), parameter :: Ny_global = 1, Nz = 3
    type(tDM) :: dm2d, dm2d_dist
    type(tPetscSF) :: migration_sf
    type(tPetscSection) :: parSection, r0Section
    real(ireals), allocatable :: arr(:, :)
    real(ireals), pointer :: xarr(:), xxarr(:, :)
    type(tVec) :: gVec
    type(tVec) :: r0Vec

    integer(iintegers) :: i, fStart, fEnd, Nx_global

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    Nx_global = numnodes

    call create_2d_regular_plex(comm, Nx_global + 1, Ny_global + 1, dm2d, dm2d_dist, opt_migration_sf=migration_sf)
    call DMPlexGetDepthStratum(dm2d_dist, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces

    if (myid .eq. 0) then
      allocate (arr(Nz, Nx_global * Ny_global * 2))
      do i = 1, size(arr, dim=2)
        ! this assumes that petsc distribute splits the x axis uniformly...
        ! 2 faces for each proc
        arr(:, i) = real((i - 1) / fEnd, ireals) ! will be mpi_rank after the distribution
      end do
    end if

    call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, &
                              arr, parSection, gVec)

    call VecGetArrayF90(gVec, xarr, ierr); call CHKERR(ierr)
    xxarr(1:Nz, fStart:fEnd - 1) => xarr
    print *, myid, 'xxarr', xxarr(:, :)
    ! check that all are the same
    @assertEqual(xxarr(1, fStart), xxarr, sqrt(epsilon(xxarr)))
    nullify (xxarr)
    call VecRestoreArrayF90(gVec, xarr, ierr); call CHKERR(ierr)

    call plex_gVec_toZero(dm2d_dist, migration_sf, parSection, gVec, &
                          r0Section, r0Vec)

    if (myid .eq. 0) then
      call VecGetArrayF90(r0Vec, xarr, ierr); call CHKERR(ierr)
      xxarr(1:Nz, 1:Nx_global * Ny_global * 2) => xarr
      print *, myid, 'xxarr', xxarr(:, :)
      @assertEqual(arr, xxarr, sqrt(epsilon(xxarr)))
      nullify (xxarr)
      call VecRestoreArrayF90(r0Vec, xarr, ierr); call CHKERR(ierr)

    end if
  end subroutine

  @test(npes=[1])
  subroutine test_dmplex_gVec_from_f90_array(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    integer(iintegers), parameter :: Ny_global = 1, Nz = 3
    type(tDM) :: dm2d, dm2d_dist
    type(tPetscSF) :: migration_sf
    type(tPetscSection) :: parSection, r0Section
    real(ireals), allocatable :: arr(:, :)
    real(ireals), pointer :: xarr(:), xxarr(:, :)
    type(tVec) :: gVec, r0Vec

    real(ireals) :: trg
    integer(iintegers) :: i, k, Nx_global, fStart, fEnd

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    Nx_global = numnodes

    call create_2d_regular_plex(comm, Nx_global + 1, Ny_global + 1, dm2d, dm2d_dist, opt_migration_sf=migration_sf)
    call DMPlexGetDepthStratum(dm2d_dist, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces

    ! Serial Case
    if (myid .eq. 0) then
      allocate (arr(Nz, Nx_global * Ny_global * 2))
      do i = 1, size(arr, dim=2)
        do k = 1, size(arr, dim=1)
          arr(k, i) = real((i - 1) * Nz + k - 1, ireals)
        end do
      end do
    else
      allocate (arr(Nz, 0))
    end if

    call dmplex_gVec_from_f90_array(comm, arr, gVec)

    call VecGetArrayF90(gVec, xarr, ierr); call CHKERR(ierr)
    xxarr(1:Nz, 1:size(arr, dim=2)) => xarr
    print *, myid, 'xxarr', xxarr(:, :)
    @assertEqual(arr, xxarr, epsilon(xxarr))
    nullify (xxarr)
    call VecRestoreArrayF90(gVec, xarr, ierr); call CHKERR(ierr)

    call VecDestroy(gvec, ierr); call CHKERR(ierr)
    deallocate (arr)

    ! All ranks having entries
    allocate (arr(Nz, fEnd - fStart))
    do i = 1, size(arr, dim=2)
      do k = 1, size(arr, dim=1)
        arr(k, i) = real((i - 1) * Nz + k - 1 + myid * 100, ireals)
      end do
    end do

    call dmplex_gVec_from_f90_array(comm, arr, gVec)

    call VecGetArrayF90(gVec, xarr, ierr); call CHKERR(ierr)
    xxarr(1:Nz, fStart:fEnd - 1) => xarr
    print *, myid, 'xxarr', xxarr(:, :)
    @assertEqual(arr, xxarr, epsilon(xxarr))
    nullify (xxarr)
    call VecRestoreArrayF90(gVec, xarr, ierr); call CHKERR(ierr)

    ! Test if we can do scatters with this vec
    call create_plex_section(dm2d_dist, 'face_section', i1, &
                             [i0], [Nz], [i0], [i0], parSection)
    call plex_gVec_toZero(dm2d_dist, migration_sf, parSection, gVec, &
                          r0Section, r0Vec)

    call VecDestroy(gvec, ierr); call CHKERR(ierr)
    deallocate (arr)

    if (myid .eq. 0) then
      call VecGetArrayF90(r0Vec, xarr, ierr); call CHKERR(ierr)
      xxarr(1:Nz, 1:Nx_global * Ny_global * 2) => xarr
      print *, myid, 'xxarr', xxarr(:, :)
      do i = 1, size(xxarr, dim=2)
        do k = 1, size(xxarr, dim=1)
          ! this assumes that petsc distribute splits the x axis uniformly...
          ! 2 faces for each proc
          trg = real((i - 1) * Nz + k - 1 + (i - 1) / fEnd * 100, ireals)
          @assertEqual(trg, xxarr(k, i), epsilon(trg))
        end do
      end do
      nullify (xxarr)
      call VecRestoreArrayF90(r0Vec, xarr, ierr); call CHKERR(ierr)
    end if

  end subroutine
end module
