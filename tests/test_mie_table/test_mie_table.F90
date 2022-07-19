module test_mie_table

  use m_tenstream_options, only: read_commandline_options

  use m_data_parameters, only:  &
    & finalize_mpi,             &
    & init_mpi_data_parameters, &
    & iintegers,                &
    & ireals,                   &
    & mpiint

  use m_mie_tables, only: &
    & destroy_mie_table, &
    & mie_tables_init, &
    & mie_optprop, &
    & t_mie_table

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm = this%getMpiCommunicator()
    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm = this%getMpiCommunicator()
    call finalize_mpi(comm, lfinalize_mpi=.false., lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes = [1, 2])
  subroutine test_load(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    integer(mpiint) :: comm

    type(t_mie_table), allocatable :: mie_water_table

    comm = this%getMpiCommunicator()

    call mie_tables_init(comm, mie_water_table, ierr, lverbose=.true.)
    @assertEqual(0, ierr)
    @assertTrue(allocated(mie_water_table), 'water_table is expected to be allocated')
    @assertTrue(allocated(mie_water_table%wvl ), 'water_table%wvl  is expected to be allocated')
    @assertTrue(allocated(mie_water_table%reff), 'water_table%reff is expected to be allocated')
    @assertTrue(allocated(mie_water_table%qext), 'water_table%qext is expected to be allocated')
    @assertTrue(allocated(mie_water_table%w0  ), 'water_table%w0   is expected to be allocated')
    @assertTrue(allocated(mie_water_table%g   ), 'water_table%g    is expected to be allocated')
    call destroy_mie_table(mie_water_table, ierr)
    @assertEqual(0, ierr)
  end subroutine

  @test(npes = [1, 2])
  subroutine test_lookup_on_supports(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    integer(mpiint) :: comm

    type(t_mie_table), allocatable :: mie_water_table

    comm = this%getMpiCommunicator()

    call mie_tables_init(comm, mie_water_table, ierr)
    @assertEqual(0, ierr)

    call test_table(mie_water_table)
    call destroy_mie_table(mie_water_table, ierr)
    @assertEqual(0, ierr)
  contains
    subroutine test_table(table)
      type(t_mie_table), intent(in) :: table
      real(ireals) :: qext, w0, g
      integer(iintegers) :: iw, ir

      do iw = lbound(table%wvl, 1), ubound(table%wvl, 1)
        do ir = lbound(table%reff, 1), ubound(table%reff, 1)
          call mie_optprop(table, real(table%wvl(iw), ireals), real(table%reff(ir), ireals), qext, w0, g, ierr)
          @assertEqual(0, ierr)
          @assertEqual(table%qext(ir,iw), qext)
          @assertEqual(table%w0(ir,iw), w0)
          @assertEqual(table%g(ir,iw), g)
        end do
      end do
    end subroutine
  end subroutine

  @test(npes = [1, 2])
  subroutine test_lookup_between_supports(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    integer(mpiint) :: comm

    type(t_mie_table), allocatable :: mie_water_table

    comm = this%getMpiCommunicator()

    call mie_tables_init(comm, mie_water_table, ierr)
    @assertEqual(0, ierr)

    call test_table(mie_water_table)
    call destroy_mie_table(mie_water_table, ierr)
    @assertEqual(0, ierr)
  contains
    subroutine test_table(table)
      type(t_mie_table), intent(in) :: table
      real(ireals) :: wvl, reff, qext, w0, g
      integer(iintegers) :: iw, ir

      do iw = lbound(table%wvl, 1), ubound(table%wvl, 1) - 1
        do ir = lbound(table%reff, 1), ubound(table%reff, 1)
          wvl = (table%wvl(iw) + table%wvl(iw + 1)) / 2
          reff = table%reff(ir)
          call mie_optprop(table, wvl, reff, qext, w0, g, ierr)
          @assertEqual(0, ierr)
          @assertEqual((table%qext(ir,iw)+table%qext(ir,iw+1))/2, qext)
          @assertEqual((table%w0(ir,iw)+table%w0(ir,iw+1))/2, w0)
          @assertEqual((table%g(ir,iw)+table%g(ir,iw+1))/2, g)
        end do
      end do

      do iw = lbound(table%wvl, 1), ubound(table%wvl, 1)
        do ir = lbound(table%reff, 1), ubound(table%reff, 1) - 1
          wvl = table%wvl(iw)
          reff = (table%reff(ir) + table%reff(ir + 1)) / 2
          call mie_optprop(table, wvl, reff, qext, w0, g, ierr)
          @assertEqual(0, ierr)
          @assertEqual((table%qext(ir,iw)+table%qext(ir+1,iw))/2, qext)
          @assertEqual((table%w0(ir,iw)+table%w0(ir+1,iw))/2, w0)
          @assertEqual((table%g(ir,iw)+table%g(ir+1,iw))/2, g)
        end do
      end do
    end subroutine
  end subroutine
end module
