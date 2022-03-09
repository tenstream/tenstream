module test_rayleigh

  use m_data_parameters, only:  &
    & iintegers,                &
    & ireals, &
    & mpiint

  use m_rayleigh, only: &
    & rayleigh

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
  end subroutine teardown

  @test(npes = [1])
  subroutine test_load(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: co2 = 400e-6
    real(ireals) :: wvl, ksca, ksca2
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    do i = 200, 10000, 100
      wvl = real(i) * 1e-3_ireals

      call rayleigh(wvl, ksca, ierr)
      @assertEqual(0, ierr, 'Should be valid output')

      call rayleigh(wvl, co2, ksca2, ierr)
      @assertEqual(0, ierr, 'Should be valid output')

      print *, 'wvl', wvl, '[mu]', ' bodhaine', ksca, 'bodhaine_co2', ksca2
    end do

    call rayleigh(.1_ireals, ksca, ierr)
    @assertTrue(ierr.ne.0, 'Should throw an error for small wavelengths')
  end subroutine

end module
