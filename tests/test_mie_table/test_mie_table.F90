module test_mie_table

  use m_tenstream_options, only : read_commandline_options

  use m_data_parameters, only : &
    & finalize_mpi,             &
    & init_mpi_data_parameters, &
    & mpiint

  use m_mie_tables, only: init_mie_tables, water_table, ice_table

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


  @test(npes = [1, 2])
  subroutine test_load(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    integer(mpiint) :: comm

    comm = this%getMpiCommunicator()

    call init_mie_tables(comm, ierr, lverbose=.True.)
    @assertEqual(0, ierr)
    @assertTrue(allocated(water_table), 'water_table is expected to be allocated')
    @assertTrue(allocated(water_table%wvl ), 'water_table%wvl  is expected to be allocated')
    @assertTrue(allocated(water_table%reff), 'water_table%reff is expected to be allocated')
    @assertTrue(allocated(water_table%qext), 'water_table%qext is expected to be allocated')
    @assertTrue(allocated(water_table%w0  ), 'water_table%w0   is expected to be allocated')
    @assertTrue(allocated(water_table%g   ), 'water_table%g    is expected to be allocated')
    !@assertTrue(allocated(ice_table), 'ice_table is expected to be allocated')
  end subroutine


end module
