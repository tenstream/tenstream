module test_rayleigh

  use m_data_parameters, only:  &
    & finalize_mpi,             &
    & init_mpi_data_parameters, &
    & iintegers,                &
    & ireals,                   &
    & mpiint

  use m_rayleigh, only: &
    & rayleigh

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call finalize_mpi(this%getMpiCommunicator(), lfinalize_mpi=.false., lfinalize_petsc=.true.)
  end subroutine teardown

  @test(npes = [1])
  subroutine test_load(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: co2 = 400e-6
    real(ireals) :: wvl, ksca, ksca2
    integer(iintegers) :: i
    integer(mpiint) :: ierr
    real(ireals), parameter :: releps = 1

    do i = 100, 5000, 100
      wvl = real(i) * 1e-3_ireals

      call rayleigh(wvl, ksca, ierr)
      if (ierr .eq. 0) then
        @assertTrue(ksca.gt.0, 'should be positive')
      end if

      call rayleigh(wvl, co2, ksca2, ierr)
      if (ierr .eq. 0) then
        @assertTrue(ksca2.gt.0, 'should be positive')
      end if

      print *, 'wvl', wvl, '[mu]', ' bodhaine', ksca, 'bodhaine_co2', ksca2
    end do

    call rayleigh(.1_ireals, ksca, ierr)
    @assertTrue(ierr.ne.0, 'Should throw an error for small wavelengths')

    call rayleigh(.5_ireals, ksca, ierr)
    @assertEqual(6.66140033E-27_ireals, ksca, ksca*releps)

    call rayleigh(.5_ireals, co2, ksca2, ierr)
    @assertEqual(6.66123701E-27_ireals, ksca2, ksca*releps)
  end subroutine

end module
