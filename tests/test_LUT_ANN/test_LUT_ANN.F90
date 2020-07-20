module test_LUT_ANN
  use m_data_parameters, only : ireals, irealLUT, iintegers, mpiint, default_str_len
  use m_optprop_ANN, only: t_ANN, load_ANN, view_ANN

#include "petsc/finclude/petsc.h"
  use petsc


  use pfunit_mod
  implicit none

  integer(mpiint) :: myid, numnodes, comm, ierr

contains

  @test( npes=[1] )
  subroutine test_LUT_direct_coeff(this)
      class (MpiTestMethod), intent(inout) :: this

      type(t_ANN), allocatable :: ann

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      call load_ANN('test_ANN_diffuse.nc', ann, ierr)
      @assertEqual(0_mpiint, ierr)

      call view_ANN(ann, ierr)
  endsubroutine

end module
