module test_LUT_wedge_18_8
  use m_boxmc, only : t_boxmc_wedge_18_8
  use m_data_parameters, only : ireals, irealLUT, ireal_params, ireal_dp, &
    iintegers, mpiint, &
    init_mpi_data_parameters, default_str_len, &
    i1, i2, i3, i4, i5
  use m_optprop_LUT, only : t_optprop_LUT_wedge_18_8, find_lut_dim_by_name, &
    azimuth_from_param_phi, param_phi_from_azimuth
  use m_optprop, only : t_optprop_wedge_18_8
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse, CHKERR, get_arg, itoa, &
    ind_nd_to_1d, ind_1d_to_nd, rad2deg, deg2rad
  use m_search, only: find_real_location
  use m_boxmc_geometry, only : setup_default_wedge_geometry

#include "petsc/finclude/petsc.h"
  use petsc

  use pfunit_mod
  implicit none

  integer(iintegers),parameter :: Ndir=18, Ndiff=8

  real(irealLUT) :: bg(3)
  real(irealLUT) :: S(Ndiff),T(Ndir), Stol(Ndiff), Ttol(Ndir)
  real(ireals)   :: S_target(Ndiff), T_target(Ndir)
  real(ireals)   :: S_tol(Ndiff),T_tol(Ndir)

  real(ireals)   :: BMC_diff2diff(Ndiff**2), BMC_dir2diff(Ndir*Ndiff), BMC_dir2dir(Ndir**2)
  real(irealLUT) :: LUT_diff2diff(Ndiff**2), LUT_dir2diff(Ndir*Ndiff), LUT_dir2dir(Ndir**2)
  real(irealLUT) :: OPP_diff2diff(Ndiff**2), OPP_dir2diff(Ndir*Ndiff), OPP_dir2dir(Ndir**2)

  type(t_boxmc_wedge_18_8) :: bmc
  type(t_optprop_wedge_18_8) :: OPP
  type(t_optprop_LUT_wedge_18_8) :: OPPLUT

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(irealLUT),parameter :: atol=1e-3, rtol=1e-1, one=1, zero=0
  real(irealLUT),parameter :: sigma = 2 ! normal test range for coefficients

  integer(iintegers), parameter :: NSLICE = 3


  integer(mpiint) :: ierr

contains

  @before
  subroutine setup(this)
      class (MpiTestMethod), intent(inout) :: this

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      if(myid.eq.0) &
        print *,'Testing LUT coefficients against BOXMC with tolerances atol/rtol ::',atol,rtol,':: on',numnodes,'ranks'

      PETSC_COMM_WORLD = comm
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

      call init_mpi_data_parameters(comm)
      call read_commandline_options(comm)

      call bmc%init(comm)

      S_target = zero
      T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      call OPPLUT%destroy()
      call OPP%destroy()
      call PetscFinalize(ierr)
  end subroutine teardown

  !@test(npes=[1])
  subroutine test_LUT_wedge_custom1(this)
    class (MpiTestMethod), intent(inout) :: this
    real(irealLUT) :: tau, w0, g, aspect, phi, param_phi, theta, Cx, Cy
    real(irealLUT), dimension(Ndir**2) :: d2d1, d2d4
    !real(irealLUT), dimension(Ndir**2) :: d2d2, d2d3

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call OPP%init(comm)
    call OPPLUT%init(comm)

    tau = 1e-3_irealLUT
    w0  = 0._irealLUT
    g   = 0._irealLUT
    aspect = .5_irealLUT
    theta = 10._irealLUT
    Cx = .5_irealLUT
    Cy = 0.8660254037844386_irealLUT

    phi = -0._irealLUT
    param_phi = real(param_phi_from_azimuth(deg2rad(real(phi, ireal_params)), real([Cx, Cy], ireal_params)), irealLUT)
    call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, param_phi, theta], d2d1)
    call print_dir2dir(d2d1)

    !call compute_bmc_coeff(d2d2)
    !call print_dir2dir(d2d2)

    !call compute_bmc_coeff(d2d3)
    !call print_dir2dir(d2d3)

    call bmc_through_sample_pts([tau, w0, aspect, Cx, Cy, phi, theta], d2d4)
    call print_dir2dir(d2d4)

    contains

      subroutine bmc_through_sample_pts(sample_pts, d2d)
        real(irealLUT), intent(in) :: sample_pts(:)
        real(irealLUT), intent(out) :: d2d(:)

        real(ireals) :: pti(size(sample_pts))
        integer(iintegers) :: rev_pti(size(sample_pts))
        integer(iintegers) :: kdim, ind1d

        real(ireal_dp), allocatable :: vertices(:)
        real(irealLUT) :: tauz, w0, g, phi, theta
        logical, parameter :: dir=.True.
        logical :: lvalid_entry
        integer(iintegers) :: isrc

        do kdim = 1, size(sample_pts)
          pti(kdim) = find_real_location(OPPLUT%dirconfig%dims(kdim)%v, sample_pts(kdim))
        enddo
        ind1d = ind_nd_to_1d(OPPLUT%dirconfig%offsets, nint(pti, kind=iintegers))
        call ind_1d_to_nd(OPPLUT%dirconfig%offsets, ind1d, rev_pti)
        print *,'pti', pti, 'ind1d', ind1d, 'revpti', rev_pti


        call OPPLUT%LUT_bmc_wrapper_determine_sample_pts(OPPLUT%dirconfig, ind1d, dir, &
          vertices, tauz, w0, g, phi, theta, lvalid_entry)

        do isrc = 1, Ndir
          call OPPLUT%bmc_wrapper(isrc, vertices, &
            tauz, w0, g, &
            dir, phi, theta, &
            this%getMpiCommunicator(), &
            S, T, Stol, Ttol)
          d2d(isrc:Ndir**2:Ndir) = real(T, irealLUT)
        enddo
      end subroutine

      !subroutine compute_bmc_wrapper_coeff(d2d)
      !  real(irealLUT), intent(out) :: d2d(:)
      !  integer(iintegers) :: isrc
      !  real(ireals), allocatable :: vertices(:)

      !  call setup_default_wedge_geometry(&
      !    [0._ireals, 0._ireals], &
      !    [1._ireals, 0._ireals], &
      !    real([Cx, Cy], ireals), &
      !    real(aspect, ireals), vertices)

      !  do isrc = 1, Ndir
      !    call OPPLUT%bmc_wrapper(isrc, vertices, &
      !      real(tau, ireals), real(w0, ireals), 0._ireals, &
      !      .True., real(phi, ireals), real(theta,ireals), &
      !      this%getMpiCommunicator(), S, T, Stol, Ttol, real(atol, ireals), real(rtol,ireals))
      !    d2d(isrc:Ndir**2:Ndir) = real(T, irealLUT)
      !  enddo
      !end subroutine

      !subroutine compute_bmc_coeff(d2d)
      !  real(irealLUT), intent(out) :: d2d(:)
      !  integer(iintegers) :: isrc
      !  real(ireals), allocatable :: vertices(:)
      !  real(irealLUT) :: kabs, ksca

      !  call setup_default_wedge_geometry(&
      !    [0._ireals, 0._ireals], &
      !    [1._ireals, 0._ireals], &
      !    real([Cx, Cy], ireals), &
      !    real(aspect, ireals), vertices)

      !  do isrc = 1, Ndir
      !    kabs = tau*(one-w0)
      !    ksca = tau*w0

      !    call bmc%get_coeff(this%getMpiCommunicator(),&
      !      real([kabs,ksca,g], ireals), &
      !      isrc, .True., &
      !      real(phi, ireals), &
      !      real(theta, ireals), &
      !      real(vertices,ireals), &
      !      S_target,T_target,S_tol,T_tol, &
      !      inp_atol=real(atol, ireals), inp_rtol=real(rtol, ireals))
      !    d2d(isrc:Ndir**2:Ndir) = real(T_target, irealLUT)
      !  enddo
      !end subroutine

      subroutine print_dir2dir(d2d)
        real(irealLUT), intent(in) :: d2d(:)
        integer(iintegers) :: isrc
        print *,'INPUT:', [tau, w0, aspect, Cx, Cy, phi, theta]
        do isrc = i1, Ndir
          print *,': src '//itoa(isrc), d2d(isrc:Ndir**2:Ndir)
        enddo
      end subroutine
  end subroutine

  @test(npes=[2,1])
  subroutine test_LUT_wedge_direct_coeff_onsamplepts(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: isrc
      integer(iintegers) :: idim_tau, idim_w0, idim_g, idim_aspect, idim_param_phi, idim_theta, idim_Cx, idim_Cy
      integer(iintegers) :: itau, iw0, ig, iaspect, iparam_phi, itheta, iCx, iCy
      real(irealLUT) :: tau, w0, g, aspect, param_phi, phi, theta, Cx, Cy

      real(irealLUT) :: kabs, ksca, dz, err(2)
      real(ireals), allocatable :: vertices(:)
      real(irealLUT), parameter :: dx = 911
      real(irealLUT),allocatable :: g_dim(:)

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      call OPP%init(comm)
      call OPPLUT%init(comm)

      associate( LUTconfig => OPPLUT%dirconfig )

      idim_tau    = find_lut_dim_by_name(LUTconfig, 'tau')
      idim_w0     = find_lut_dim_by_name(LUTconfig, 'w0')
      idim_g      = find_lut_dim_by_name(LUTconfig, 'g')
      if(idim_g.eq.-1) then
        allocate(g_dim(1), source=zero)
      else
        allocate(g_dim(LUTconfig%dims(idim_g)%N), source=LUTconfig%dims(idim_g     )%v(:))
      endif
      idim_aspect = find_lut_dim_by_name(LUTconfig, 'aspect_zx')
      idim_param_phi    = find_lut_dim_by_name(LUTconfig, 'param_phi')
      idim_theta  = find_lut_dim_by_name(LUTconfig, 'theta')
      idim_Cx     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cx')
      idim_Cy     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cy')

      do itau = 1, LUTconfig%dims(idim_tau)%N, NSLICE
        do iw0  = 1, LUTconfig%dims(idim_w0)%N, NSLICE
          do ig   = 1, size(g_dim), NSLICE
            do iaspect = 1, LUTconfig%dims(idim_aspect)%N, NSLICE
              do iparam_phi = 1, LUTconfig%dims(idim_param_phi)%N, NSLICE
                do itheta = 1, LUTconfig%dims(idim_theta)%N, NSLICE
                  do iCx = 1, LUTconfig%dims(idim_Cx)%N, NSLICE
                    do iCy = 1, LUTconfig%dims(idim_Cy)%N, NSLICE
                      tau    = LUTconfig%dims(idim_tau   )%v(itau)
                      w0     = LUTconfig%dims(idim_w0    )%v(iw0)
                      g      = g_dim(ig)
                      aspect = LUTconfig%dims(idim_aspect)%v(iaspect)
                      param_phi = LUTconfig%dims(idim_param_phi   )%v(iparam_phi)
                      theta  = LUTconfig%dims(idim_theta )%v(itheta)
                      Cx     = LUTconfig%dims(idim_Cx    )%v(iCx)
                      Cy     = LUTconfig%dims(idim_Cy    )%v(iCy)
                      phi    =real(rad2deg(azimuth_from_param_phi(&
                        real(param_phi, ireal_params), real([Cx,Cy], ireal_params))), irealLUT)

                      call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, param_phi, theta], LUT_dir2dir)
                      call OPPLUT%LUT_get_dir2diff([tau, w0, aspect, Cx, Cy, param_phi, theta], LUT_dir2diff)

                      call setup_default_wedge_geometry(&
                        [0._ireals, 0._ireals], &
                        [1._ireals, 0._ireals], &
                        real([Cx, Cy], ireals), &
                        real(aspect, ireals), vertices)

                      vertices = vertices * dx
                      do isrc = 1, Ndir
                        dz = real(vertices(12)-vertices(3), irealLUT)
                        kabs = (one-w0) * tau / dz
                        ksca = w0 * tau / dz

                        call bmc%get_coeff(comm,&
                          real([kabs,ksca,g], ireal_dp), &
                          isrc, .True., &
                          real(phi, ireal_dp), &
                          real(theta, ireal_dp), &
                          real(vertices,ireal_dp), &
                          S_target,T_target,S_tol,T_tol, &
                          inp_atol=real(atol, ireal_dp), inp_rtol=real(rtol, ireal_dp))

                        err = rmse(LUT_dir2dir(isrc:Ndir**2:Ndir), real(T_target, irealLUT))
                        if(err(1).ge.sigma*atol .or. err(2).ge.sigma*rtol) then
                          print *,'Testing: ', tau, w0, g, aspect, param_phi, phi, theta, Cx, Cy,':: RMSE', err
                          print *,'LUT :::', isrc, LUT_dir2dir(isrc:Ndir**2:Ndir)
                          print *,'dir2dir', isrc, T_target
                        endif
                        call check(S_target, T_target, LUT_dir2diff(isrc:Ndir*Ndiff:Ndir), LUT_dir2dir(isrc:Ndir**2:Ndir), &
                          msg='test_LUT_wedge_direct_coeff_onsamplepts')
                      enddo !isrc
                    enddo !Cy
                  enddo !Cx
                enddo !theta
              enddo !param_phi
            enddo !aspect
          enddo !g
        enddo !w0
      enddo !tau

      end associate


  endsubroutine

  subroutine check(S_target,T_target, S,T, msg, opt_atol, opt_rtol)
      real(ireals),intent(in),dimension(:) :: S_target,T_target
      real(irealLUT),intent(in),dimension(:) :: S,T
      character(len=*),optional :: msg
      real(irealLUT), intent(in), optional :: opt_atol, opt_rtol

      character(default_str_len) :: local_msgS, local_msgT
      logical, parameter :: ldetail=.False.

      real(irealLUT) :: arg_atol, arg_rtol, Terr(2), Serr(2)

      arg_atol = get_arg(atol*sigma, opt_atol)
      arg_rtol = get_arg(rtol*sigma, opt_rtol)

      Terr = rmse(T,real(T_target, irealLUT))
      Serr = rmse(S,real(S_target, irealLUT))

      if(myid.eq.0) then
        if(ldetail) then
          print*,''

          if( present(msg) ) then
            write(local_msgS,*) trim(msg),':: Diffuse boxmc coefficient not as '
            write(local_msgT,*) trim(msg),':: Direct  boxmc coefficient not as '
            print *,msg
          else
            write(local_msgS,*) 'Diffuse boxmc coefficient not as '
            write(local_msgT,*) 'Direct  boxmc coefficient not as '
          endif

          print*,'---------------------'
          write(*, FMT='( " diffuse :: ",10(f12.5) )' ) S
          print*,''
          write(*, FMT='( " target  :: ",10(f12.5) )' ) S_target
          print*,''
          write(*, FMT='( " diff    :: ",10(f12.5) )' ) S_target-S
          print*,'RMSE ::: ',Serr*[one, 100._irealLUT],'%'
          print*,''
          write(*, FMT='( " direct  :: ", 8(f12.5) )' ) T
          print*,''
          write(*, FMT='( " target  :: ", 8(f12.5) )' ) T_target
          print*,''
          write(*, FMT='( " diff    :: ", 8(f12.5) )' ) T_target-T
          print*,'RMSE ::: ',Terr*[one, 100._irealLUT],'%'
          print*,'---------------------'
          print*,''
        else
          if(Terr(1).gt.arg_atol .or. Serr(1).gt.arg_atol) then
            print*,'RMSE ::: ',rmse(S,real(S_target, irealLUT))*[one, 100._irealLUT],'% ', &
                   'direct:::',rmse(T,real(T_target, irealLUT))*[one, 100._irealLUT],'%'
            print *,''
          endif
        endif

        @assertEqual(real(S_target, irealLUT), S, arg_atol, local_msgS )
        @assertLessThanOrEqual   (zero, S)
        @assertGreaterThanOrEqual(one , S)

        @assertEqual(real(T_target, irealLUT), T, arg_atol, local_msgT )
        @assertLessThanOrEqual   (zero, T)
        @assertGreaterThanOrEqual(one , T)
      endif
  end subroutine
end module
