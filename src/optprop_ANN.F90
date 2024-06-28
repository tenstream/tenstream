!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_optprop_ANN
  use m_data_parameters, only: ireals, irealLUT, iintegers, i1, mpiint, default_str_len
  use m_optprop_parameters, only: ldebug_optprop, lut_basename
  use m_optprop_base, only: t_optprop_base, t_op_config, set_op_param_space, &
    & print_op_config, check_if_samplepts_in_bounds
  use m_netcdfio, only: get_attribute, ncload
  use mpi
  use m_helper_functions, only: imp_bcast, CHKERR, CHKWARN, toStr, char_to_upper, expm1, get_petsc_opt
  use m_search, only: find_real_location
  use m_eddington, only: eddington_coeff_ec

  use m_boxmc, only: t_boxmc, t_boxmc_3_10

#ifdef HAVE_FORNADO
  use m_fornado_base, only: fr, fi, ferr, func_name_to_id, ANN_view
  use m_fornado_base, only: t_fornado_ANN => t_ANN, t_fornado_ANN_layer => t_ANN_layer
  use m_fornado, only: fornado_inference
#endif

  implicit none
  private
  public t_optprop_ANN, t_optprop_ANN_3_10

#ifdef HAVE_FORNADO
  type t_ANN
    type(t_fornado_ANN), allocatable :: fann
    character(len=default_str_len) :: fname
    logical :: lphysical_input
  end type
#endif

  type, abstract, extends(t_optprop_base) :: t_optprop_ANN
    character(default_str_len) :: basename
#ifdef HAVE_FORNADO
    type(t_ANN), allocatable :: ann_dir2dir
    type(t_ANN), allocatable :: ann_dir2diff
    type(t_ANN), allocatable :: ann_diff2diff
#endif
    logical :: initialized = .false.
  contains
    procedure :: init
    procedure :: destroy
    procedure :: get_dir2dir => ANN_get_dir2dir
    procedure :: get_dir2diff => ANN_get_dir2diff
    procedure :: get_diff2diff => ANN_get_diff2diff
  end type

  type, extends(t_optprop_ANN) :: t_optprop_ANN_3_10
  end type

  logical, parameter :: lrenormalize = .false.
  real(irealLUT), parameter :: renorm_eps = 1e-5_ireallut

#ifdef HAVE_FORNADO

contains

  subroutine init(ANN, comm, ierr)
    class(t_optprop_ANN), intent(inout) :: ANN
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    integer(mpiint) :: myid
    logical :: lflg

    if (ANN%initialized) return

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    select type (ANN)
    class is (t_optprop_ANN_3_10)
      ANN%dir_streams = 3
      ANN%diff_streams = 10
      !call set_op_param_space(ANN, 'LUT_3_10_for_ANN', ierr); call CHKERR(ierr)
      call set_op_param_space(ANN, 'LUT_3_10', ierr); call CHKERR(ierr)

    class default
      call CHKERR(1_mpiint, 'initialize ANN: unexpected type for optprop_ANN object!')
    end select

    ANN%basename = &
      & trim(lut_basename)!//&
    !& trim(gen_ann_basename("_ANN_", ANN%dirconfig))

    call get_petsc_opt('', "-ann_basename", ANN%basename, lflg, ierr); call CHKERR(ierr)

    print *, 'Init ANN for basename: ', trim(ANN%basename)
    print *, 'ANN dirconfig'
    call print_op_config(ANN%dirconfig)
    print *, 'ANN diffconfig'
    call print_op_config(ANN%diffconfig)

    allocate (ANN%ann_dir2dir)
    ANN%ann_dir2dir%fname = trim(ANN%basename)//"_dir2dir"//&
      & trim(gen_ann_basename("_ANN_", ANN%dirconfig))//&
      & '.c'//toStr(ANN%dir_streams)//'.nc'

    allocate (ANN%ann_dir2diff)
    ANN%ann_dir2diff%fname = trim(ANN%basename)//"_dir2diff"//&
      & trim(gen_ann_basename("_ANN_", ANN%dirconfig))//&
      & ".c"//toStr(ANN%dir_streams)// &
      & "_"//toStr(ANN%diff_streams)//'.nc'

    allocate (ANN%ann_diff2diff)
    ANN%ann_diff2diff%fname = trim(ANN%basename)//"_diff2diff"//&
      & trim(gen_ann_basename("_ANN_", ANN%diffconfig))//&
      & ".c"//toStr(ANN%diff_streams)//'.nc'

    call ANN_load(comm, ANN%ann_dir2dir, ierr); call CHKERR(ierr)
    call ANN_load(comm, ANN%ann_dir2diff, ierr); call CHKERR(ierr)
    call ANN_load(comm, ANN%ann_diff2diff, ierr); call CHKERR(ierr)

    ANN%initialized = .true.
    ierr = 0
  end subroutine

  function gen_ann_basename(prefix, config) result(lutname)
    character(len=:), allocatable :: lutname
    character(len=*), intent(in) :: prefix
    type(t_op_config), intent(in) :: config
    integer(iintegers) :: k
    lutname = trim(prefix)
    do k = 1, size(config%dims)
      lutname = trim(lutname)//'.'//trim(config%dims(k)%dimname)//toStr(config%dims(k)%N)
    end do
  end function

  subroutine destroy(ANN, ierr)
    class(t_optprop_ANN), intent(inout) :: ann
    integer(mpiint), intent(out) :: ierr
    if (allocated(ann%ann_dir2dir)) deallocate (ann%ann_dir2dir)
    if (allocated(ann%ann_dir2diff)) deallocate (ann%ann_dir2diff)
    if (allocated(ann%ann_diff2diff)) deallocate (ann%ann_diff2diff)
    ann%basename = 'deallocated'
    ann%dir_streams = -1
    ann%diff_streams = -1
    ann%initialized = .false.
    ierr = 0
  end subroutine

  subroutine ANN_load(comm, ann, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_ANN), intent(inout) :: ann
    integer(mpiint), intent(out) :: ierr

    character(len=default_str_len) :: groups(2), act_name
    integer(iintegers) :: Nlayer, i
    integer(mpiint) :: myid, iphys
    integer(ferr) :: fe

    if (allocated(ann%fann)) call CHKERR(1_mpiint, 'ANN already allocated')
    allocate (ann%fann)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    if (myid .eq. 0) then
      print *, 'loading ann from file:', ann%fname

      call get_attribute(ann%fname, "global", "physical_input", iphys, ierr)
      call CHKWARN(ierr, "Could not load global attribute if network was trained on physical_input or with indices")
      ann%lphysical_input = iphys .ne. 0_mpiint

      call get_attribute(ann%fname, "global", "Nlayer", Nlayer, ierr); call CHKERR(ierr)
      allocate (ann%fann%layers(Nlayer))

      groups(1) = trim(ann%fname)
      do i = 1, size(ann%fann%layers)
        groups(2) = 'w'//toStr(i - 1); call ncload(groups, ann%fann%layers(i)%wgt, ierr)
        call CHKERR(ierr, "Could not load ANN weights: "//trim(groups(2)))
        groups(2) = 'b'//toStr(i - 1); call ncload(groups, ann%fann%layers(i)%bias, ierr)
        call CHKERR(ierr, "Could not load ANN bias: "//trim(groups(2)))
        call get_attribute(ann%fname, 'w'//toStr(i - 1), "activation", act_name, ierr); call CHKERR(ierr)
        call char_to_upper(act_name)
        ann%fann%layers(i)%activation_id = func_name_to_id(act_name)
        if (ann%fann%layers(i)%activation_id .lt. 0) then
          call CHKERR(int(ann%fann%layers(i)%activation_id, mpiint), 'Could not find func id for activation_name '//trim(act_name))
        end if
      end do
      call ANN_view(ann%fann, fe); ierr = ierr + fe
    end if
    call scatter_net()
    ierr = 0
  contains
    subroutine scatter_net()
      call imp_bcast(comm, ann%lphysical_input, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, Nlayer, 0_mpiint, ierr); call CHKERR(ierr)
      if (.not. allocated(ann%fann%layers)) allocate (ann%fann%layers(Nlayer))
      do i = 1, size(ann%fann%layers)
        call imp_bcast(comm, ann%fann%layers(i)%wgt, 0_mpiint, ierr); call CHKERR(ierr)
        call imp_bcast(comm, ann%fann%layers(i)%bias, 0_mpiint, ierr); call CHKERR(ierr)
        call imp_bcast(comm, ann%fann%layers(i)%activation_id, 0_mpiint, ierr); call CHKERR(ierr)
      end do
    end subroutine
  end subroutine

  pure subroutine ANN_predict(ann, inp, C, ierr)
    type(t_fornado_ANN), intent(in) :: ann
    real(irealLUT), intent(in) :: inp(:)
    real(irealLUT), intent(out) :: C(:)
    integer(mpiint), intent(out) :: ierr

    real(fr) :: Cfr(size(C))
    integer(ferr) :: iferr

    call fornado_inference(ann, real(inp, fr), Cfr, iferr)
    C = Cfr
    ierr = int(iferr, mpiint)
  end subroutine

  subroutine ANN_get_dir2dir(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%dir_streams**2)

    integer(iintegers) :: src, kdim
    real(irealLUT) :: pti_buffer(size(sample_pts)+2), maxtrans
    real(irealLUT), pointer :: pC(:, :) ! dim(src, dst)
    integer(mpiint) :: ierr

    if (OPP%ann_dir2dir%lphysical_input) then
      pti_buffer = sample_pts
      !pti_buffer(1) = exp(-pti_buffer(1))
      call ANN_predict(OPP%ann_dir2dir%fann, pti_buffer, C, ierr); call CHKERR(ierr)
    else
      do kdim = 1, size(sample_pts)
        pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
      end do
      pti_buffer = pti_buffer - 1 ! because the network is trained on C indices
      call ANN_predict(OPP%ann_dir2dir%fann, pti_buffer, C, ierr); call CHKERR(ierr)
    end if

    if (ldebug_optprop) then
      !Check for energy conservation:
      ierr = 0
      pC(1:OPP%dir_streams, 1:OPP%dir_streams) => C(:)
      maxtrans = 0
      do src = 1, OPP%dir_streams
        maxtrans = max(maxtrans, sum(pC(src, :)))
        if (sum(pC(src, :)) .gt. 1._ireallut + 1e-2_ireallut) ierr = ierr + 1
      end do
      if (ierr .ne. 0) then
        do src = 1, OPP%dir_streams
          print *, 'SUM dir2dir coeff for src '//toStr(src)//' :: sum', sum(pC(src, :)), ':: C', pC(src, :)
        end do
        call CHKWARN(1_mpiint, 'Check for energy conservation failed: '//toStr(maxtrans))
      end if
    end if

    if (lrenormalize) then
      pC(1:OPP%dir_streams, 1:OPP%dir_streams) => C(:)
      pC = min(1._ireallut, max(0._ireallut, pC))
      do src = 1, OPP%dir_streams
        maxtrans = sum(pC(src, :))
        if (maxtrans .gt. 1._ireallut) pC(src, :) = pC(src, :) / (maxtrans + renorm_eps)
      end do
    end if
  end subroutine

  subroutine ANN_get_dir2diff(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%dir_streams*ANN%diff_streams)

    integer(iintegers) :: src, kdim
    real(irealLUT) :: pti_buffer(size(sample_pts)), maxtrans
    real(irealLUT), pointer :: pC(:, :) ! dim(src, dst)
    integer(mpiint) :: ierr

    if (OPP%ann_dir2diff%lphysical_input) then
      pti_buffer = sample_pts
      !pti_buffer(1) = exp(-pti_buffer(1))
      call ANN_predict(OPP%ann_dir2diff%fann, pti_buffer, C, ierr); call CHKERR(ierr)
    else
      do kdim = 1, size(sample_pts)
        pti_buffer(kdim) = find_real_location(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
      end do
      pti_buffer = pti_buffer - 1 ! because the network is trained on C indices
      call ANN_predict(OPP%ann_dir2diff%fann, pti_buffer, C, ierr); call CHKERR(ierr)
    end if

    if (ldebug_optprop) then
      !Check for energy conservation:
      ierr = 0
      pC(1:OPP%dir_streams, 1:OPP%diff_streams) => C(:)
      maxtrans = 0
      do src = 1, OPP%dir_streams
        maxtrans = max(maxtrans, sum(pC(src, :)))
        if (sum(pC(src, :)) .gt. 1._ireallut + 1e-2_ireallut) ierr = ierr + 1
      end do
      if (ierr .ne. 0) then
        do src = 1, OPP%dir_streams
          print *, 'SUM dir2diff coeff for src '//toStr(src)//' :: sum', sum(pC(src, :)), ':: C', pC(src, :)
        end do
        call CHKWARN(1_mpiint, 'Check for energy conservation failed: '//toStr(maxtrans))
      end if
    end if

    if (lrenormalize) then
      pC(1:OPP%dir_streams, 1:OPP%diff_streams) => C(:)
      pC = min(1._ireallut, max(0._ireallut, pC))
      do src = 1, OPP%dir_streams
        maxtrans = sum(pC(src, :))
        if (maxtrans .gt. 1._ireallut) pC(src, :) = pC(src, :) / (maxtrans + renorm_eps)
      end do
    end if
  end subroutine

  subroutine ANN_get_diff2diff(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%diff_streams**2)

    real(irealLUT) :: pti_buffer(15), tauz, w0, dz, g, tanx, tany, maxtrans
    real(irealLUT), parameter :: theta = 0, phi = 0, mu = cos(theta)
    real(irealLUT), target :: Cp(110)
    real(ireals) :: t, r, rdir, sdir, tdir

    integer(iintegers) :: src, kdim
    real(irealLUT), pointer :: pC(:, :) ! dim(src, dst)
    real(irealLUT), pointer :: pCp(:, :) ! dim(dst, src)
    integer(mpiint) :: ierr

    if (OPP%ann_diff2diff%lphysical_input) then
      !pti_buffer = sample_pts
      !pti_buffer(1) = exp(-pti_buffer(1))

      tauz = sample_pts(1)
      w0 = sample_pts(2)
      dz = sample_pts(3)
      g = sample_pts(4)
      tanx = 0
      tany = 0

      call eddington_coeff_ec(real(tauz, ireals), real(w0, ireals), real(g, ireals), real(mu, ireals), t, r, rdir, sdir, tdir)

      pti_buffer(:) = [real(irealLUT) :: &
        & (tauz / dz) * (1.-w0),  & ! kabs
        & (tauz / dz) * w0,       & ! ksca
        & -expm1(-tauz * (1.-w0)), &
        & -expm1(-tauz * w0),      &
        & g,                      &
        & dz,                     &
        & phi,                    &
        & theta,                  &
        & tanx,                   &
        & tany,                   &
        & real(t, irealLUT),   &
        & real(r, irealLUT),   &
        & real(rdir, irealLUT),   &
        & real(sdir, irealLUT),   &
        & real(tdir, irealLUT)]
      !print *, 'inp ANN', sample_pts
      !print *, 'pti ANN', pti_buffer
      call ANN_predict(OPP%ann_diff2diff%fann, pti_buffer, Cp, ierr); call CHKERR(ierr)
      !print *,'coeff ANN', Cp
      pCp(1:OPP%diff_streams, 1:OPP%diff_streams + 1) => Cp
      pC(1:OPP%diff_streams, 1:OPP%diff_streams) => C(:)
      !do src = 1, OPP%diff_streams+1
      !  print *, 'src', src, pCp(:,src), 'norm', sum(pCp(:, src))
      !enddo
      do src = 1, OPP%diff_streams
        pC(src, :) = max(0._ireallut, min(1._ireallut, pCp(:, src)))
        pC(src, :) = pC(src, :) / sum(pC(src, :)) * max(0._ireallut, min(1._ireallut, pCp(src, OPP%diff_streams + 1)))
        !print *, 'src', src, pC(src,:), 'norm', sum(pC(src,:))
      end do
    else

      if (ldebug_optprop) then
        if (size(sample_pts) .ne. size(OPP%diffconfig%dims)) then
          call print_op_config(OPP%diffconfig)
          call CHKERR(1_mpiint, 'size of sample_pts array ne number of dimensions in ANN ' &
                      //toStr(size(sample_pts, kind=iintegers))//'/'//toStr(size(OPP%diffconfig%dims)))
        end if
        call check_if_samplepts_in_bounds(sample_pts, OPP%diffconfig)
      end if

      do kdim = 1, size(sample_pts)
        pti_buffer(kdim) = find_real_location(OPP%diffconfig%dims(kdim)%v, sample_pts(kdim))
      end do
      pti_buffer = pti_buffer - 1 ! because the network is trained on C indices
      call ANN_predict(OPP%ann_diff2diff%fann, pti_buffer, C, ierr); call CHKERR(ierr)
    end if

    if (ldebug_optprop) then
      !Check for energy conservation:
      ierr = 0
      pC(1:OPP%diff_streams, 1:OPP%diff_streams) => C(:)
      !pC = min(1._ireallut, max(0._ireallut, pC))
      maxtrans = 0
      do src = 1, OPP%diff_streams
        maxtrans = max(maxtrans, sum(pC(src, :)))
        if (sum(pC(src, :)) .gt. 1._ireallut + 1e-2_ireallut) ierr = ierr + 1
      end do
      if (ierr .ne. 0) then
        do src = 1, OPP%diff_streams
          print *, 'SUM diff2diff coeff for src '//toStr(src)//' :: sum', sum(pC(src, :)), ':: C', pC(src, :)
        end do
        call CHKWARN(1_mpiint, 'Check for energy conservation failed: '//toStr(maxtrans))
      end if
    end if

    if (lrenormalize) then
      pC(1:OPP%diff_streams, 1:OPP%diff_streams) => C(:)
      do src = 1, OPP%diff_streams
        maxtrans = sum(pC(src, :))
        if (maxtrans .gt. 1._ireallut) pC(src, :) = pC(src, :) / (maxtrans + renorm_eps)
      end do
    end if
  end subroutine

! else HAVE_FORNADO
#else

  character(len=*), parameter :: no_avail_errmsg = &
    & "TenStream package is not configured with FORNADO, you cannot use ANN as LUT replacement."// &
    & "Try to compile with -DBUILD_FORNADO"

contains
  subroutine init(ANN, comm, ierr)
    class(t_optprop_ANN), intent(inout) :: ANN
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    ierr = 10 + comm
    call CHKERR(ierr, no_avail_errmsg)
    if (ANN%initialized) return
  end subroutine
  subroutine destroy(ANN, ierr)
    class(t_optprop_ANN), intent(inout) :: ann
    integer(mpiint), intent(out) :: ierr
    ierr = 1
    call CHKERR(ierr, no_avail_errmsg)
    ANN%initialized = .false.
  end subroutine
  subroutine ANN_get_dir2dir(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%dir_streams**2)
    C(:) = sample_pts(1) + real(OPP%dir_streams)
    call CHKERR(1_mpiint, no_avail_errmsg)
  end subroutine

  subroutine ANN_get_dir2diff(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%dir_streams*ANN%diff_streams)
    C(:) = sample_pts(1) + real(OPP%dir_streams)
    call CHKERR(1_mpiint, no_avail_errmsg)
  end subroutine

  subroutine ANN_get_diff2diff(OPP, sample_pts, C)
    class(t_optprop_ANN), intent(in) :: OPP
    real(irealLUT), intent(in) :: sample_pts(:)
    real(irealLUT), target, intent(out) :: C(:) ! dimension(ANN%diff_streams**2)
    C(:) = sample_pts(1) + real(OPP%dir_streams)
    call CHKERR(1_mpiint, no_avail_errmsg)
  end subroutine

! endif HAVE_FORNADO
#endif
end module
