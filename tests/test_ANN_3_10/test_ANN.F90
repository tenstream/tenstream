module test_ANN
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only : init_mpi_data_parameters, &
    & ireals, irealLUT, iintegers, mpiint, default_str_len
  use m_helper_functions, only: CHKERR, colored_str_by_range
  use m_netcdfio, only: ncwrite
  use m_optprop, only: t_optprop_3_10, t_optprop_3_10_ann

#include "petsc/finclude/petsc.h"
  use petsc


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
    logical :: lpetsc_is_initialized
    integer(mpiint) :: ierr
    call PetscInitialized(lpetsc_is_initialized, ierr)
    if(lpetsc_is_initialized) call PetscFinalize(ierr)
  end subroutine teardown

  @test( npes=[1] )
  subroutine test_compare_diff2diff_to_LUT(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    type(t_optprop_3_10)     :: OPP_LUT
    type(t_optprop_3_10_ann) :: OPP_ANN

    integer(iintegers) :: isrc
    real(irealLUT) :: tauz, w0, g, aspect
    real(irealLUT), target, allocatable :: C_LUT(:)
    real(irealLUT), target, allocatable :: C_ANN(:)
    real(irealLUT), target, allocatable :: C_rel_err(:)
    real(irealLUT), pointer :: pC_LUT(:,:)
    real(irealLUT), pointer :: pC_ANN(:,:)

    real(irealLUT), parameter :: climits(7) = [0., 1e-8, 1e-3, 1e-2, 1e-1, .5, 1.]
    character(len=*), parameter :: colors(6) = [character(len=6) :: 'black', 'green', 'purple', 'peach', 'blue', 'red']

    real(irealLUT), parameter :: rlimits(6) = [-.5, -.1, -1e-2, 1e-2, .1, .5] * 100
    character(len=*), parameter :: rolors(5) = [character(len=6) :: 'red', 'purple', 'green', 'purple', 'red']

    print *,"Checking ANN for diff2diff coeffs"

    comm     = this%getMpiCommunicator()

    call OPP_LUT%init(comm)
    call OPP_ANN%init(comm)

    associate( Ndiff => OPP_LUT%LUT%diff_streams )
      allocate(C_LUT(Ndiff**2), C_ANN(Ndiff**2), C_rel_err(Ndiff))
      pC_LUT(1:Ndiff, 1:Ndiff) => C_LUT(:)
      pC_ANN(1:Ndiff, 1:Ndiff) => C_ANN(:)

      ! idx 22 5 14 4
      tauz = 1.08083180518_irealLUT
      w0   = .521358613652_irealLUT
      aspect = 1._irealLUT
      g    = .5717_irealLUT

      call OPP_LUT%LUT%get_diff2diff([tauz, w0, aspect, g], C_LUT)
      do isrc = 1, Ndiff
        print *,'C_LUT src', isrc, ':', trim(colored_str_by_range(pC_LUT(isrc,:), limits=climits, colors=colors)), sum(pC_LUT(isrc,:))
      enddo
      print *,''

      call OPP_ANN%ANN%get_diff2diff([tauz, w0, aspect, g], C_ANN)
      do isrc = 1, Ndiff
        print *,'C_ANN src', isrc, ':', trim(colored_str_by_range(pC_ANN(isrc,:), limits=climits, colors=colors)), sum(pC_ANN(isrc,:))
      enddo
      print *,''
      print *,'Relative error: [%] .................. bias'

      call OPP_ANN%ANN%get_diff2diff([tauz, w0, aspect, g], C_ANN)
      do isrc = 1, Ndiff
        where(pC_LUT(isrc,:).gt.0._irealLUT)
          C_rel_err(:) = pC_ANN(isrc,:)/pC_LUT(isrc,:)
        elsewhere
          C_rel_err(:) = pC_ANN(isrc,:)
        end where

        print *,'C_ANN src', isrc, ':', &
          & trim(colored_str_by_range(100._irealLUT - C_rel_err(:)*100, limits=rlimits, colors=rolors)), &
          & sum(pC_LUT(isrc,:)) - sum(pC_ANN(isrc,:))
      enddo

      print *,'RMSE', sqrt(sum((C_LUT-C_ANN)**2)/size(C_LUT))

    end associate
    call OPP_LUT%destroy()
    call OPP_ANN%destroy()
  endsubroutine

  @test( npes=[1] )
  subroutine test_compare_dir2diff_to_LUT(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    type(t_optprop_3_10)     :: OPP_LUT
    type(t_optprop_3_10_ann) :: OPP_ANN

    integer(iintegers) :: isrc
    real(irealLUT) :: tauz, w0, g, aspect, phi, theta, inp(6)
    real(irealLUT), target, allocatable :: C_LUT(:)
    real(irealLUT), target, allocatable :: C_ANN(:)
    real(irealLUT), target, allocatable :: C_rel_err(:)
    real(irealLUT), pointer :: pC_LUT(:,:)
    real(irealLUT), pointer :: pC_ANN(:,:)

    real(irealLUT), parameter :: climits(7) = [0., 1e-8, 1e-3, 1e-2, 1e-1, .5, 1.]
    character(len=*), parameter :: colors(6) = [character(len=6) :: 'black', 'green', 'purple', 'peach', 'blue', 'red']

    real(irealLUT), parameter :: rlimits(6) = [-.5, -.1, -1e-2, 1e-2, .1, .5] * 100
    character(len=*), parameter :: rolors(5) = [character(len=6) :: 'red', 'purple', 'green', 'purple', 'red']

    print *,"Checking ANN for dir2diff coeffs"

    comm     = this%getMpiCommunicator()

    call OPP_LUT%init(comm)
    call OPP_ANN%init(comm)

    associate( Ndir => OPP_LUT%LUT%dir_streams, Ndiff => OPP_LUT%LUT%diff_streams )
      allocate(C_LUT(Ndir*Ndiff), C_ANN(Ndir*Ndiff), C_rel_err(Ndiff))
      pC_LUT(1:Ndir, 1:Ndiff) => C_LUT(:)
      pC_ANN(1:Ndir, 1:Ndiff) => C_ANN(:)

      ! idx 22 5 14 4
      tauz = 1.08083180518_irealLUT
      w0   = .521358613652_irealLUT
      aspect = 1._irealLUT
      g    = .5717_irealLUT
      phi  = 10
      theta= 10
      inp = [tauz, w0, aspect, g, phi, theta]

      call OPP_LUT%LUT%get_dir2diff(inp, C_LUT)
      do isrc = 1, Ndir
        print *,'C_LUT src', isrc, ':', trim(colored_str_by_range(pC_LUT(isrc,:), limits=climits, colors=colors)), sum(pC_LUT(isrc,:))
      enddo
      print *,''

      call OPP_ANN%ANN%get_dir2diff(inp, C_ANN)
      do isrc = 1, Ndir
        print *,'C_ANN src', isrc, ':', trim(colored_str_by_range(pC_ANN(isrc,:), limits=climits, colors=colors)), sum(pC_ANN(isrc,:))
      enddo
      print *,''
      print *,'Relative error: [%] .................. bias'

      call OPP_ANN%ANN%get_dir2diff(inp, C_ANN)
      do isrc = 1, Ndir
        where(pC_LUT(isrc,:).gt.0._irealLUT)
          C_rel_err(:) = pC_ANN(isrc,:)/pC_LUT(isrc,:)
        elsewhere
          C_rel_err(:) = pC_ANN(isrc,:)
        end where
        print *,'C_ANN src', isrc, ':', &
          & trim(colored_str_by_range(100._irealLUT - C_rel_err(:)*100, limits=rlimits, colors=rolors)), &
          & sum(pC_LUT(isrc,:)) - sum(pC_ANN(isrc,:))
      enddo

      print *,'RMSE', sqrt(sum((C_LUT-C_ANN)**2)/size(C_LUT))

    end associate
    call OPP_LUT%destroy()
    call OPP_ANN%destroy()
  endsubroutine

  @test( npes=[1] )
  subroutine test_compare_dir2dir_to_LUT(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm

    type(t_optprop_3_10)     :: OPP_LUT
    type(t_optprop_3_10_ann) :: OPP_ANN

    integer(iintegers) :: isrc
    real(irealLUT) :: tauz, w0, g, aspect, phi, theta, inp(6)
    real(irealLUT), target, allocatable :: C_LUT(:)
    real(irealLUT), target, allocatable :: C_ANN(:)
    real(irealLUT), target, allocatable :: C_rel_err(:)
    real(irealLUT), pointer :: pC_LUT(:,:)
    real(irealLUT), pointer :: pC_ANN(:,:)

    real(irealLUT), parameter :: climits(7) = [0., 1e-8, 1e-3, 1e-2, 1e-1, .5, 1.]
    character(len=*), parameter :: colors(6) = [character(len=6) :: 'black', 'green', 'purple', 'peach', 'blue', 'red']

    real(irealLUT), parameter :: rlimits(6) = [-.5, -.1, -1e-2, 1e-2, .1, .5] * 100
    character(len=*), parameter :: rolors(5) = [character(len=6) :: 'red', 'purple', 'green', 'purple', 'red']

    print *,"Checking ANN for dir2dir coeffs"

    comm     = this%getMpiCommunicator()

    call OPP_LUT%init(comm)
    call OPP_ANN%init(comm)

    associate( Ndir => OPP_LUT%LUT%dir_streams )
      allocate(C_LUT(Ndir**2), C_ANN(Ndir**2), C_rel_err(Ndir))
      pC_LUT(1:Ndir, 1:Ndir) => C_LUT(:)
      pC_ANN(1:Ndir, 1:Ndir) => C_ANN(:)

      ! idx 22 5 14 4
      tauz = 1.08083180518_irealLUT
      w0   = .521358613652_irealLUT
      aspect = 1._irealLUT
      g    = .5717_irealLUT
      phi  = 10
      theta= 10
      inp = [tauz, w0, aspect, g, phi, theta]

      call OPP_LUT%LUT%get_dir2dir(inp, C_LUT)
      do isrc = 1, Ndir
        print *,'C_LUT src', isrc, ':', trim(colored_str_by_range(pC_LUT(isrc,:), limits=climits, colors=colors)), sum(pC_LUT(isrc,:))
      enddo
      print *,''

      call OPP_ANN%ANN%get_dir2dir(inp, C_ANN)
      do isrc = 1, Ndir
        print *,'C_ANN src', isrc, ':', trim(colored_str_by_range(pC_ANN(isrc,:), limits=climits, colors=colors)), sum(pC_ANN(isrc,:))
      enddo
      print *,''
      print *,'Relative error: [%] .................. bias'

      call OPP_ANN%ANN%get_dir2dir(inp, C_ANN)
      do isrc = 1, Ndir
        where(pC_LUT(isrc,:).gt.0._irealLUT)
          C_rel_err(:) = pC_ANN(isrc,:)/pC_LUT(isrc,:)
        elsewhere
          C_rel_err(:) = pC_ANN(isrc,:)
        end where
        print *,'C_ANN src', isrc, ':', &
          & trim(colored_str_by_range(100._irealLUT - C_rel_err(:)*100, limits=rlimits, colors=rolors)), &
          & sum(pC_LUT(isrc,:)) - sum(pC_ANN(isrc,:))
      enddo

      print *,'RMSE', sqrt(sum((C_LUT-C_ANN)**2)/size(C_LUT))

    end associate
    call OPP_LUT%destroy()
    call OPP_ANN%destroy()
  endsubroutine

end module
