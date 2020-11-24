module m_optprop_base

  use m_data_parameters, only : &
    & iintegers, mpiint, &
    & ireals, irealLUT, ireal_dp, &
    & default_str_len, inil, i1, i2, i3

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & linspace, &
    & meanval, &
    & ndarray_offsets, &
    & toStr

  use m_boxmc, only: t_boxmc

  use m_optprop_parameters, only:         &
    & stddev_atol, stddev_rtol

  use m_tenstream_options, only: lLUT_mockup

  use m_optprop_parameters, only:         &
    LUT_MAX_DIM,                          &
    preset_g4,                            &
    preset_g6,                            &
    preset_param_phi6,                    &
    preset_param_phi11,                   &
    preset_param_phi19,                   &
    preset_param_theta4,                  &
    preset_param_theta13,                 &
    preset_aspect5,                       &
    preset_aspect7,                       &
    preset_aspect11,                      &
    preset_aspect13,                      &
    preset_aspect17,                      &
    preset_aspect18,                      &
    preset_aspect23,                      &
    preset_w010,                          &
    preset_w020,                          &
    preset_w021,                          &
    preset_tau15,                         &
    preset_tau20,                         &
    preset_tau31

  implicit none

  type t_optprop_dim
    integer(iintegers) :: N      ! size of dimension
    character(len=default_str_len) :: dimname
    real(irealLUT) :: vrange(2) ! min / max of dimension
    real(irealLUT), allocatable :: v(:) ! sampling points of dimension, size N
  end type

  type t_op_config
    type(t_optprop_dim), allocatable :: dims(:)
    integer(iintegers),allocatable :: offsets(:) ! offsets of respective dimensions (starts at 0 and next one is dim(1)%N ... etc... )
  end type

  type,abstract :: t_optprop_base
    type(t_op_config), allocatable :: dirconfig, diffconfig
    integer(iintegers) :: dir_streams = inil, diff_streams = inil
    class(t_boxmc), allocatable :: bmc
    contains
      procedure :: bmc_wrapper
      procedure(get_dir2dir  ), deferred :: get_dir2dir
      procedure(get_dir2diff ), deferred :: get_dir2diff
      procedure(get_diff2diff), deferred :: get_diff2diff
  end type

  abstract interface
    subroutine get_dir2dir(OPP, sample_pts, C)
      import t_optprop_base, irealLUT
      class(t_optprop_base), intent(in) :: OPP
      real(irealLUT),intent(in) :: sample_pts(:)
      real(irealLUT), target, intent(out):: C(:)
    end subroutine
  end interface
  abstract interface
    subroutine get_dir2diff(OPP, sample_pts, C)
      import t_optprop_base, irealLUT
      class(t_optprop_base), intent(in) :: OPP
      real(irealLUT),intent(in) :: sample_pts(:)
      real(irealLUT), target, intent(out):: C(:)
    end subroutine
  end interface
  abstract interface
    subroutine get_diff2diff(OPP, sample_pts, C)
      import t_optprop_base, irealLUT
      class(t_optprop_base), intent(in) :: OPP
      real(irealLUT),intent(in) :: sample_pts(:)
      real(irealLUT), target, intent(out):: C(:)
    end subroutine
  end interface

  logical, parameter :: ldebug=.False.

  interface populate_op_dim
    module procedure populate_op_dim_vrange
    module procedure populate_op_dim_preset
  end interface
contains

  pure subroutine populate_op_dim_vrange(dimname, N, op_dim, vrange)
    character(len=*),intent(in) :: dimname
    integer(iintegers), intent(in) :: N
    type(t_optprop_dim),intent(out) :: op_dim
    real(irealLUT), intent(in) :: vrange(:)
    integer(iintegers) :: k
    if(allocated(op_dim%v)) return ! already done
    allocate(op_dim%v(N))
    do k=1,N
      op_dim%v(k) = linspace(k, vrange, N)
    enddo
    op_dim%vrange = [op_dim%v(1), op_dim%v(N)]
    op_dim%dimname = trim(dimname)
    op_dim%N = size(op_dim%v)
  end subroutine

  pure subroutine populate_op_dim_preset(dimname, op_dim, preset)
    character(len=*),intent(in) :: dimname
    type(t_optprop_dim),intent(out) :: op_dim
    real(irealLUT), intent(in) :: preset(:)
    if(allocated(op_dim%v)) return ! already done

    op_dim%N = size(preset)
    allocate(op_dim%v(op_dim%N))
    op_dim%v = preset
    op_dim%vrange = [op_dim%v(1), op_dim%v(op_dim%N)]
    op_dim%dimname = trim(dimname)
  end subroutine


  subroutine print_op_config(config)
    type(t_op_config), allocatable, intent(in) :: config
    integer(iintegers) :: i

    print *,'LUT Config initialized', allocated(config)
    if(.not.allocated(config)) return

    print *,'LUT Config Ndim', size(config%dims)
    do i = 1, size(config%dims)
      print *,'Dimension '//toStr(i)//' '//trim(config%dims(i)%dimname)// &
        ' size '//toStr(config%dims(i)%N), '(', config%dims(i)%vrange, ')'
    enddo
  end subroutine

  subroutine set_op_param_space(OPP, cname, ierr)
  class(t_optprop_base) :: OPP
    character(len=*) :: cname
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    allocate(OPP%dirconfig)
    allocate(OPP%diffconfig)

    if(.not.lLUT_mockup) then
      select case(cname)
      case('LUT_1_2')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       i1, OPP%dirconfig%dims(5), vrange=real([0], irealLUT)) ! azimithally average in 1D
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_8_10')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_8_12')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       i2, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     i3, OPP%dirconfig%dims(6), vrange=real([55,60,65], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_8_16')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_8_18')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_3_10')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_3_10_for_ANN')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau15)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w010)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect13)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g4)
        call populate_op_dim('phi',       7_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     7_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_3_16')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_3_6')
        allocate(OPP%dirconfig%dims(6))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('phi',       19_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
        call populate_op_dim('theta',     19_iintegers, OPP%dirconfig%dims(6), vrange=real([0,90], irealLUT))
        allocate(OPP%diffconfig%dims(4))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w020)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)

      case('LUT_wedge_5_8')
        allocate(OPP%dirconfig%dims(8))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau15)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w010)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect18)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('wedge_coord_Cx', 7_iintegers, OPP%dirconfig%dims(5), &
          vrange=real([.35,.65], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 7_iintegers, OPP%dirconfig%dims(6), &
          vrange=real([0.7760254, 0.9560254], irealLUT))

        call populate_op_dim('param_phi', &
          OPP%dirconfig%dims(7), preset=preset_param_phi19)
        call populate_op_dim('param_theta', &
          OPP%dirconfig%dims(8), preset=preset_param_theta13)

        allocate(OPP%diffconfig%dims(6))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w021)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)
        call populate_op_dim('wedge_coord_Cx', 7_iintegers, OPP%diffconfig%dims(5), &
          vrange=real([.35,.65], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 7_iintegers, OPP%diffconfig%dims(6), &
          vrange=real([0.7760254, 0.9560254], irealLUT))

      case('LUT_rectilinear_wedge_5_8')
        allocate(OPP%dirconfig%dims(8))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w021)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect18)
        call populate_op_dim('g',         OPP%dirconfig%dims(4), preset=preset_g6)
        call populate_op_dim('wedge_coord_Cx', 3_iintegers, OPP%dirconfig%dims(5), &
          vrange=real([-0.000001,1.000001], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 2_iintegers, OPP%dirconfig%dims(6), &
          vrange=real([.499999,1.000001], irealLUT))

        call populate_op_dim('param_phi', &
          OPP%dirconfig%dims(7), preset=preset_param_phi19)
        call populate_op_dim('param_theta', &
          OPP%dirconfig%dims(8), preset=preset_param_theta13)

        allocate(OPP%diffconfig%dims(6))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau31)
        call populate_op_dim('w0',        OPP%diffconfig%dims(2), preset=preset_w021)
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect23)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)
        call populate_op_dim('wedge_coord_Cx', 3_iintegers, OPP%diffconfig%dims(5), &
          vrange=real([-0.000001,1.000001], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 2_iintegers, OPP%diffconfig%dims(6), &
          vrange=real([.499999,1.000001], irealLUT))

      case('LUT_wedge_18_8')
        allocate(OPP%dirconfig%dims(7))
        call populate_op_dim('tau',       OPP%dirconfig%dims(1), preset=preset_tau15)
        call populate_op_dim('w0',        OPP%dirconfig%dims(2), preset=preset_w010)
        call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=preset_aspect18)
        call populate_op_dim('wedge_coord_Cx', 7_iintegers, OPP%dirconfig%dims(4), &
          vrange=real([.35,.65], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 7_iintegers, OPP%dirconfig%dims(5), &
          vrange=real([0.7760254, 0.9560254], irealLUT))

        call populate_op_dim('param_phi', &
          OPP%dirconfig%dims(6), preset=preset_param_phi19)
        call populate_op_dim('param_theta', &
          OPP%dirconfig%dims(7), preset=preset_param_theta13)

        allocate(OPP%diffconfig%dims(6))
        call populate_op_dim('tau',       OPP%diffconfig%dims(1), preset=preset_tau20)
        call populate_op_dim('w0',        i2, OPP%diffconfig%dims(2), vrange=real([.0,.99999], irealLUT))
        call populate_op_dim('aspect_zx', OPP%diffconfig%dims(3), preset=preset_aspect7)
        call populate_op_dim('g',         OPP%diffconfig%dims(4), preset=preset_g6)
        call populate_op_dim('wedge_coord_Cx', 7_iintegers, OPP%diffconfig%dims(5), &
          vrange=real([.35,.65], irealLUT))
        call populate_op_dim('wedge_coord_Cy', 7_iintegers, OPP%diffconfig%dims(6), &
          vrange=real([0.7760254, 0.9560254], irealLUT))

      case default
        call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop object name: '//trim(cname))
      end select

    else ! do a mockup of a LUT
      select case(cname)
      case('LUT_1_2')
        call setup_pprts_mockup()
      case('LUT_8_10')
        call setup_pprts_mockup()
      case('LUT_8_12')
        call setup_pprts_mockup()
      case('LUT_8_16')
        call setup_pprts_mockup()
      case('LUT_8_18')
        call setup_pprts_mockup()
      case('LUT_3_10')
        call setup_pprts_mockup()
      case('LUT_3_16')
        call setup_pprts_mockup()
      case('LUT_3_6')
        call setup_pprts_mockup()
      case('LUT_wedge_5_8')
        call setup_wedge_mockup([real(irealLUT) :: .35,.65], [real(irealLUT) :: 0.7760254, 0.9560254])
      case('LUT_rectilinear_wedge_5_8')
        call setup_wedge_mockup([real(irealLUT) :: -0.000001,1.000001], [real(irealLUT) :: .499999,1.000001])
      case('LUT_wedge_18_8')
        call setup_wedge_mockup([real(irealLUT) :: .35,.65], [real(irealLUT) :: 0.7760254, 0.9560254])
      case default
        call CHKERR(1_mpiint, 'set_parameter space: unexpected type for optprop object name: '//trim(cname))
      end select
    endif

    if(size(OPP%dirconfig%dims).gt.LUT_MAX_DIM) call CHKERR(1_mpiint, 'Parameter LUT_MAX_DIM is too small '//&
      'for the LUT you are using... please increase it to '// &
      toStr(size(OPP%dirconfig%dims))//' in src/optprop_parameters.f90')
    if(size(OPP%diffconfig%dims).gt.LUT_MAX_DIM) call CHKERR(1_mpiint, 'Parameter LUT_MAX_DIM is too small '//&
      'for the LUT you are using... please increase it to '// &
      toStr(size(OPP%diffconfig%dims))//' in src/optprop_parameters.f90')

    ! Determine offsets
    allocate(OPP%dirconfig%offsets(size(OPP%dirconfig%dims)))
    call ndarray_offsets(OPP%dirconfig%dims(:)%N, OPP%dirconfig%offsets)

    allocate(OPP%diffconfig%offsets(size(OPP%diffconfig%dims)))
    call ndarray_offsets(OPP%diffconfig%dims(:)%N, OPP%diffconfig%offsets)

    !    if(ldebug.and..False.) then
    !      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr); call CHKERR(ierr)
    !      print *,myid,'set_parameter space dims:', size(OPP%diffconfig%dims), size(OPP%dirconfig%dims)
    !      do k=1,size(OPP%dirconfig%dims)
    !        print *,myid,'dim ',trim(OPP%dirconfig%dims(k)%dimname), OPP%dirconfig%offsets(k), &
    !          OPP%dirconfig%dims(k)%vrange, ':', OPP%dirconfig%dims(k)%v
    !      enddo
    !    endif
  contains
    subroutine setup_pprts_mockup()
      allocate(OPP%dirconfig%dims(6))
      call populate_op_dim('tau',       2_iintegers, OPP%dirconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_op_dim('w0',        2_iintegers, OPP%dirconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_op_dim('aspect_zx', OPP%dirconfig%dims(3), preset=[real(irealLUT) :: 1e-2, 1., 2.])
      call populate_op_dim('g',         2_iintegers, OPP%dirconfig%dims(4), vrange=real([0,1], irealLUT))
      call populate_op_dim('phi',       2_iintegers, OPP%dirconfig%dims(5), vrange=real([0,90], irealLUT))
      call populate_op_dim('theta',     OPP%dirconfig%dims(6), preset=[ real(irealLUT) :: 0, 45, 90])
      allocate(OPP%diffconfig%dims(4))
      call populate_op_dim('tau',       2_iintegers, OPP%diffconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_op_dim('w0',        2_iintegers, OPP%diffconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_op_dim('aspect_zx', 2_iintegers, OPP%diffconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_op_dim('g',         2_iintegers, OPP%diffconfig%dims(4), vrange=real([0,1], irealLUT))
    end subroutine
    subroutine setup_wedge_mockup(Cx, Cy)
      real(irealLUT), intent(in) :: Cx(2), Cy(2)
      allocate(OPP%dirconfig%dims(7))
      call populate_op_dim('tau',       2_iintegers, OPP%dirconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_op_dim('w0',        2_iintegers, OPP%dirconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_op_dim('aspect_zx', 2_iintegers, OPP%dirconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_op_dim('wedge_coord_Cx', 3_iintegers, OPP%dirconfig%dims(4), vrange=Cx)
      call populate_op_dim('wedge_coord_Cy', 2_iintegers, OPP%dirconfig%dims(5), vrange=Cy)
      call populate_op_dim('param_phi'  , OPP%dirconfig%dims(6), preset=preset_param_phi6)
      call populate_op_dim('param_theta', OPP%dirconfig%dims(7), preset=preset_param_theta4)
      allocate(OPP%diffconfig%dims(5))
      call populate_op_dim('tau',       2_iintegers, OPP%diffconfig%dims(1), vrange=real([1e-30,1e3], irealLUT))
      call populate_op_dim('w0',        2_iintegers, OPP%diffconfig%dims(2), vrange=real([0,1], irealLUT))
      call populate_op_dim('aspect_zx', 2_iintegers, OPP%diffconfig%dims(3), vrange=real([1e-2,2.], irealLUT))
      call populate_op_dim('wedge_coord_Cx', 3_iintegers, OPP%diffconfig%dims(4), vrange=Cx)
      call populate_op_dim('wedge_coord_Cy', 2_iintegers, OPP%diffconfig%dims(5), vrange=Cy)
    end subroutine
  end subroutine

  ! return the integer in config%dims that corresponds to the given dimension
  function find_op_dim_by_name(config, dimname) result(kdim)
    type(t_op_config), intent(in) :: config
    character(len=*), intent(in) :: dimname
    integer(iintegers) :: kdim

    integer(iintegers) :: k
    do k=1,size(config%dims)
      if(trim(dimname).eq.trim(config%dims(k)%dimname)) then
        kdim = k
        return
      endif
    enddo
    kdim=-1
  end function

  subroutine get_sample_pnt_by_name_and_index(config, dimname, index_1d, sample_pnt, ierr)
    type(t_op_config), intent(in) :: config
    character(len=*), intent(in) :: dimname
    integer(iintegers), intent(in) :: index_1d
    real(irealLUT), intent(inout) :: sample_pnt
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: kdim, nd_indices(size(config%dims))

    call ind_1d_to_nd(config%offsets, index_1d, nd_indices)

    kdim = find_op_dim_by_name(config, trim(dimname))
    if(kdim.lt.i1) then ! could not find the corresponding dimension
      ierr = 1
      return
    endif
    if(nd_indices(kdim).gt.size(config%dims(kdim)%v)) then
      print *,index_1d,'nd_indices', nd_indices
      call CHKERR(1_mpiint, 'wrong indices in kdim')
    endif
    sample_pnt = config%dims(kdim)%v(nd_indices(kdim))
    ierr = 0
  end subroutine

  subroutine check_if_samplepts_in_bounds(sample_pts, config)
    real(irealLUT),intent(in) :: sample_pts(:)
    type(t_op_config), allocatable, intent(in) :: config
    integer(mpiint) :: ierr, kdim

    ierr = 0
    do kdim = 1,size(sample_pts)
      if(sample_pts(kdim).lt.config%dims(kdim)%vrange(1).or.sample_pts(kdim).gt.config%dims(kdim)%vrange(2)) then
        print *,'ERROR value in dimension '//trim(config%dims(kdim)%dimname)// &
          ' ('//toStr(kdim)//') is outside of LUT range', &
          sample_pts(kdim), 'not in:', config%dims(kdim)%vrange
        ierr = ierr +1
      endif
    enddo

    if(ierr.ne.0) then
      call print_op_config(config)
      call CHKERR(ierr, 'Out of Bounds ERROR in LUT retrieval')
    endif
  end subroutine

  subroutine bmc_wrapper(OPP, src, vertices, tauz, w0, g, dir, phi, theta, comm, &
      S_diff, T_dir, S_tol, T_tol, inp_atol, inp_rtol)
    class(t_optprop_base), intent(in) :: OPP
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: dir
    integer(mpiint),intent(in) :: comm
    real(irealLUT), intent(in) :: tauz, w0, g, phi, theta
    real(ireal_dp) :: mean_height, vertices(:)

    real(irealLUT),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(irealLUT),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)
    real(irealLUT),intent(in),optional :: inp_atol, inp_rtol
    real(ireals) :: rS_diff(OPP%diff_streams),rT_dir(OPP%dir_streams)
    real(ireals) :: rS_tol (OPP%diff_streams),rT_tol(OPP%dir_streams)

    real(ireal_dp) :: bg(3), dz, atol, rtol

    atol = get_arg(stddev_atol, inp_atol)
    rtol = get_arg(stddev_rtol, inp_rtol)

    atol = atol-epsilon(atol)*10
    rtol = rtol-epsilon(rtol)*10

    mean_height = meanval(vertices(3:size(vertices):3))
    dz = real(meanval(2*abs(mean_height - vertices(3:size(vertices):3))), kind(dz))

    bg(1) = tauz / dz * (1._irealLUT-w0)
    bg(2) = tauz / dz * w0
    bg(3) = g

    !print *,comm,'BMC :: calling bmc_get_coeff tauz',tauz,'w0,g',w0,g,phi,theta
    !print *,comm,'BMC :: calling bmc_get_coeff dz bg',vertices(size(vertices)),bg, &
    !  & '=>', sum(bg(1:2))*vertices(size(vertices)),'/',tauz
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(1:3) ,':', vertices(10:12)
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(4:6) ,':', vertices(13:15)
    !print *,comm,'BMC :: calling bmc_get_coeff verts', vertices(7:9) ,':', vertices(16:18)
    !print *,'area bot', triangle_area_by_vertices(vertices(1:3), vertices(4:6), vertices(7:9))
    !print *,'area top', triangle_area_by_vertices(vertices(10:12), vertices(13:15), vertices(16:18))
    !print *,'area_ratio:', triangle_area_by_vertices(vertices(1:3), vertices(4:6), vertices(7:9)) &
    !  / triangle_area_by_vertices(vertices(10:12), vertices(13:15), vertices(16:18))
    !call CHKERR(1_mpiint, 'DEBUG')

    call OPP%bmc%get_coeff(comm, bg, &
      src, dir, &
      real(phi, ireal_dp), real(theta, ireal_dp), &
      vertices,   &
      rS_diff, rT_dir, &
      rS_tol, rT_tol, &
      inp_atol=atol, inp_rtol=rtol)

    S_diff = real(rS_diff, irealLUT)
    T_dir  = real(rT_dir, irealLUT)
    S_tol  = real(rS_tol, irealLUT)
    T_tol  = real(rT_tol, irealLUT)
  end subroutine
end module
