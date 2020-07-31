module test_ANN
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only : ireals, irealLUT, iintegers, mpiint, default_str_len
  use m_optprop_ANN, only: t_ANN, ANN_load, ANN_destroy, ANN_predict
  use m_helper_functions, only: CHKERR, colored_str_by_range
  use m_netcdfio, only: ncwrite
  use m_optprop, only: t_optprop_3_10, t_optprop_3_10_ann

#include "petsc/finclude/petsc.h"
  use petsc


  use pfunit_mod
  implicit none

contains
  @test( npes=[1] )
  subroutine test_load_ANN(this)
      class (MpiTestMethod), intent(inout) :: this

      type(t_ANN), allocatable :: ann
      integer(mpiint) :: myid, numnodes, comm, ierr

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      allocate(ann)
      ann%fname = 'test_ANN_diffuse.nc'
      call ANN_load(ann, ierr)
      @assertEqual(0_mpiint, ierr)

      @assertEqual(3_iintegers, size(ann%fann%layers, kind=iintegers), 'wrong number of layers')
      call ANN_destroy(ann, ierr)
      @assertEqual(0_mpiint, ierr)
  endsubroutine

  @test( npes=[1] )
  subroutine test_compare_to_LUT(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(mpiint) :: myid, numnodes, comm

      type(t_optprop_3_10)     :: OPP_LUT
      type(t_optprop_3_10_ann) :: OPP_ANN

      integer(iintegers) :: idst
      real(irealLUT) :: tauz, w0, g, aspect
      real(irealLUT), target, allocatable :: C_LUT(:)
      real(irealLUT), target, allocatable :: C_ANN(:)
      real(irealLUT), pointer :: pC_LUT(:,:)
      real(irealLUT), pointer :: pC_ANN(:,:)

      real(irealLUT), parameter :: climits(7) = [0., 1e-8, 1e-3, 1e-2, 1e-1, .5, 1.]
      character(len=*), parameter :: colors(6) = [character(len=6) :: 'black', 'green', 'purple', 'peach', 'blue', 'red']

      !real(ireals)   :: BMC_diff2diff(Ndiff*Ndiff), BMC_dir2diff(Ndir*Ndiff), BMC_dir2dir(Ndir*Ndir)

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      call read_commandline_options(comm)

      call OPP_LUT%init(comm)
      call OPP_ANN%init(comm)

      associate( Ndiff => OPP_LUT%LUT%diff_streams )
        allocate(C_LUT(Ndiff**2), C_ANN(Ndiff**2))
        pC_LUT(1:Ndiff, 1:Ndiff) => C_LUT(:)
        pC_ANN(1:Ndiff, 1:Ndiff) => C_ANN(:)

        ! idx 22 5 14 4
        tauz = 1.08083180518_irealLUT
        w0   = .521358613652_irealLUT
        aspect = 1._irealLUT
        g    = .5717_irealLUT

        call OPP_LUT%LUT%get_diff2diff([tauz, w0, aspect, g], C_LUT)
        do idst = 1, Ndiff
          print *,'C_LUT dst', idst, ':', trim(colored_str_by_range(pC_LUT(:, idst), limits=climits, colors=colors))
        enddo
        print *,''

        call OPP_ANN%ANN%get_diff2diff([tauz, w0, aspect, g], C_ANN)
        do idst = 1, Ndiff
          print *,'C_ANN dst', idst, ':', trim(colored_str_by_range(pC_ANN(:, idst), limits=climits, colors=colors))
        enddo

      end associate
      call OPP_LUT%destroy()
      call OPP_ANN%destroy()
  endsubroutine

end module
