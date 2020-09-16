module test_buildings
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only : &
    init_mpi_data_parameters,   &
    finalize_mpi,               &
    iintegers, ireals, mpiint,  &
    zero, one, default_str_len
  use m_helper_functions, only : CHKERR

  ! main entry point for solver, and desctructor
!  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

!  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

!  use m_pprts_base, only : t_solver_3_10
  use m_buildings, only: &
    & t_pprts_buildings, &
    & init_buildings,    &
    & clone_buildings,   &
    & destroy_buildings

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call finalize_mpi(this%getMpiCommunicator())
  end subroutine teardown

  @test(npes =[1])
  subroutine test_init_buildings(this)
    class (MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    integer(iintegers), parameter :: Nfaces=6
    integer(iintegers), parameter :: da_sizes(4) = [ integer(iintegers) :: 6, 1, 1, 1 ]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings (B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated (B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    enddo
    call check_iface_data_consistency(B1)

    call destroy_buildings(B1, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B1))

    contains
      subroutine check_iface_data_consistency(B)
        type(t_pprts_buildings), intent(in) :: B
        integer(iintegers) :: i
        do i = 1, Nfaces
          @assertEqual(i, B%iface(i))
          @assertEqual(real(i, ireals), B%albedo(i))
        enddo
      end subroutine
  end subroutine

  @test(npes =[1])
  subroutine test_clone_buildings(this)
    class (MpiTestMethod), intent(inout) :: this

    type(t_pprts_buildings), allocatable :: B1
    type(t_pprts_buildings), allocatable :: B2
    type(t_pprts_buildings), allocatable :: B3
    integer(iintegers), parameter :: Nfaces=6
    integer(iintegers), parameter :: da_sizes(4) = [ integer(iintegers) :: 6, 1, 1, 1 ]
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call init_buildings (B1, da_sizes, Nfaces, ierr); call CHKERR(ierr)
    @assertTrue(allocated (B1%albedo))
    @assertTrue(associated(B1%iface))
    @assertTrue(associated(B1%iface_data%data))
    @assertEqual(Nfaces, size(B1%iface_data%data, kind=iintegers))
    @assertEqual(1_iintegers, B1%iface_data%ref_count)
    do i = 1, Nfaces
      B1%iface(i) = i
      B1%albedo(i) = real(i, ireals)
    enddo
    call check_iface_data_consistency(B1)

    call clone_buildings(B1, B2, .True., ierr); call CHKERR(ierr)
    @assertEqual(2_iintegers, B1%iface_data%ref_count, 'origin ref_count did not increase')

    @assertTrue(allocated (B2%albedo))
    @assertTrue(associated(B2%iface))
    @assertTrue(associated(B2%iface_data%data))
    @assertEqual(Nfaces, size(B2%iface_data%data, kind=iintegers))
    @assertEqual(2_iintegers, B2%iface_data%ref_count)
    @assertEqual(B1%iface, B2%iface)
    @assertTrue(associated(B1%iface, B2%iface), 'cloned iface pointer dont point to same target')
    call check_iface_data_consistency(B2)

    call clone_buildings(B1, B3, .False., ierr); call CHKERR(ierr)
    @assertEqual(3_iintegers, B1%iface_data%ref_count, 'origin ref_count did not increase')

    @assertFalse(allocated(B3%albedo))
    @assertTrue(associated(B3%iface))
    @assertTrue(associated(B3%iface_data%data))
    @assertEqual(Nfaces, size(B3%iface_data%data, kind=iintegers))
    @assertEqual(3_iintegers, B3%iface_data%ref_count)
    @assertEqual(B1%iface, B3%iface)
    @assertTrue(associated(B1%iface, B3%iface), 'cloned iface pointer dont point to same target')
    call check_iface_data_consistency(B3)

    call destroy_buildings(B1, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B1))

    call check_iface_data_consistency(B2)
    call check_iface_data_consistency(B3)

    call destroy_buildings(B2, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B2))
    call destroy_buildings(B3, ierr); call CHKERR(ierr)
    @assertFalse(allocated (B3))

    contains
      subroutine check_iface_data_consistency(B)
        type(t_pprts_buildings), intent(in) :: B
        integer(iintegers) :: i
        do i = 1, Nfaces
          @assertEqual(i, B%iface(i), "iface does not match")
          if(allocated(B%albedo)) then
            @assertEqual(real(i, ireals), B%albedo(i), "albedo does not match")
          endif
        enddo
      end subroutine
  end subroutine
end module
