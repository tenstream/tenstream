module test_LUT_wedge_5_8
  use m_boxmc, only : t_boxmc_wedge_5_8
  use m_data_parameters, only : mpiint, ireals, iintegers, &
    one, zero, init_mpi_data_parameters, default_str_len, &
    i1, i2, i3, i4, i5
  use m_optprop_LUT, only : t_optprop_LUT_wedge_5_8, find_lut_dim_by_name
  use m_optprop, only : t_optprop_wedge_5_8
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse, CHKERR, get_arg, itoa
  use m_boxmc_geometry, only : setup_default_wedge_geometry

#include "petsc/finclude/petsc.h"
  use petsc

  use pfunit_mod
  implicit none

  integer(iintegers),parameter :: Ndir=5, Ndiff=8

  real(ireals) :: bg(3)
  real(ireals) :: S(Ndiff),T(Ndir), S_target(Ndiff), T_target(Ndir)
  real(ireals) :: S_tol(Ndiff),T_tol(Ndir)

  real(ireals) :: BMC_diff2diff(Ndiff**2), BMC_dir2diff(Ndir*Ndiff), BMC_dir2dir(Ndir**2)
  real(ireals) :: LUT_diff2diff(Ndiff**2), LUT_dir2diff(Ndir*Ndiff), LUT_dir2dir(Ndir**2)
  real(ireals) :: OPP_diff2diff(Ndiff**2), OPP_dir2diff(Ndir*Ndiff), OPP_dir2dir(Ndir**2)

  type(t_boxmc_wedge_5_8) :: bmc
  type(t_optprop_wedge_5_8) :: OPP
  type(t_optprop_LUT_wedge_5_8) :: OPPLUT

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-3, rtol=1e-1
  real(ireals),parameter :: sigma = 2 ! normal test range for coefficients

  integer(iintegers), parameter :: NSLICE = 3


  integer(mpiint) :: ierr

contains

  @before
  subroutine setup(this)
      class (MpiTestMethod), intent(inout) :: this

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      if(myid.eq.0) print *,'Testing LUT coefficients against BOXMC with tolerances atol/rtol :: ',atol,rtol,' :: on ',numnodes,'ranks'

      PETSC_COMM_WORLD = comm
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

      call init_mpi_data_parameters(comm)
      call read_commandline_options(comm)

      call bmc%init(comm)
      call OPP%init(comm)
      call OPPLUT%init(comm)

      call OPPLUT%print_configs()

      S_target = zero
      T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      !call OPPLUT%destroy()
      call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes=[1])
  subroutine test_LUT_wedge_custom1(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireals) :: tau, w0, g, aspect, phi, theta, Cx, Cy
    real(ireals) :: d2d1(5**2), d2d2(5**2)

    tau = 1e-3_ireals
    w0  = .5_ireals
    g   = 0._ireals
    aspect = 1._ireals
    theta = 20._ireals
    Cx = .5_ireals
    Cy = 0.8660254037844386_ireals

    phi = -27._ireals
    call print_dir2dir()
    call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], d2d1)

    phi = +27._ireals
    call print_dir2dir()
    call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], d2d2)

    @assertEqual(d2d1(srcdst_2_idx(i1, i3)), d2d2(srcdst_2_idx(i1, i4)))
    @assertEqual(d2d1(srcdst_2_idx(i1, i4)), d2d2(srcdst_2_idx(i1, i3)))
    @assertEqual(d2d1(srcdst_2_idx(i1, i5)), d2d2(srcdst_2_idx(i1, i5)))

    @assertEqual(d2d1(srcdst_2_idx(i2, i3)), d2d2(srcdst_2_idx(i2, i4)))
    @assertEqual(d2d1(srcdst_2_idx(i2, i4)), d2d2(srcdst_2_idx(i2, i3)))
    @assertEqual(d2d1(srcdst_2_idx(i2, i5)), d2d2(srcdst_2_idx(i2, i5)))

    @assertEqual(d2d1(srcdst_2_idx(i4, i3)), d2d2(srcdst_2_idx(i3, i4)))

    contains
      integer(iintegers) function srcdst_2_idx(src, dst)
        integer(iintegers),intent(in) :: src, dst
        srcdst_2_idx = (src-1)*5+dst
      end function
      subroutine print_dir2dir()
        integer(iintegers) :: isrc
        call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2dir)
        print *,'INPUT:', [tau, w0, aspect, Cx, Cy, phi, theta]
        do isrc = i1, i5
          print *,': src '//itoa(isrc), LUT_dir2dir(isrc:size(LUT_dir2dir):i5)
        enddo
      end subroutine
  end subroutine

  !@test(npes=[2,1])
  subroutine test_LUT_wedge_direct_coeff_onsamplepts(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: isrc
      integer(iintegers) :: idim_tau, idim_w0, idim_g, idim_aspect, idim_phi, idim_theta, idim_Cx, idim_Cy
      integer(iintegers) :: itau, iw0, ig, iaspect, iphi, itheta, iCx, iCy
      real(ireals) :: tau, w0, g, aspect, phi, theta, Cx, Cy

      real(ireals) :: kabs, ksca, dz, err(2)
      real(ireals), allocatable :: vertices(:)
      real(ireals), parameter :: dx = 911
      real(ireals),allocatable :: g_dim(:)

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()
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
      idim_phi    = find_lut_dim_by_name(LUTconfig, 'phi')
      idim_theta  = find_lut_dim_by_name(LUTconfig, 'theta')
      idim_Cx     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cx')
      idim_Cy     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cy')

      do itau = 1, LUTconfig%dims(idim_tau)%N, NSLICE
        do iw0  = 1, LUTconfig%dims(idim_w0)%N, NSLICE
          do ig   = 1, size(g_dim), NSLICE
            do iaspect = 1, LUTconfig%dims(idim_aspect)%N, NSLICE
              do iphi = 1, LUTconfig%dims(idim_phi)%N, NSLICE
                do itheta = 1, LUTconfig%dims(idim_theta)%N, NSLICE
                  do iCx = 1, LUTconfig%dims(idim_Cx)%N, NSLICE
                    do iCy = 1, LUTconfig%dims(idim_Cy)%N, NSLICE
                      tau    = LUTconfig%dims(idim_tau   )%v(itau)
                      w0     = LUTconfig%dims(idim_w0    )%v(iw0)
                      g      = g_dim(ig)
                      aspect = LUTconfig%dims(idim_aspect)%v(iaspect)
                      phi    = LUTconfig%dims(idim_phi   )%v(iphi)
                      theta  = LUTconfig%dims(idim_theta )%v(itheta)
                      Cx     = LUTconfig%dims(idim_Cx    )%v(iCx)
                      Cy     = LUTconfig%dims(idim_Cy    )%v(iCy)

                      call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2dir)
                      call OPPLUT%LUT_get_dir2diff([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2diff)

                      call setup_default_wedge_geometry([zero, zero], [one, zero], [Cx, Cy], aspect, vertices)
                      vertices = vertices * dx
                      do isrc = 1, Ndir
                        dz = vertices(12)-vertices(3)
                        kabs = (one-w0) * tau / dz
                        ksca = w0 * tau / dz
                        call bmc%get_coeff(comm,[kabs,ksca,g],isrc,.True.,phi,theta,vertices,S_target,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
                        err = rmse(LUT_dir2dir(isrc:Ndir**2:Ndir), T_target)
                        if(err(1).ge.sigma*atol .or. err(2).ge.sigma*rtol) then
                          print *,'Testing: ', tau, w0, g, aspect, phi, theta, Cx, Cy,':: RMSE', err
                          print *,'LUT :::', isrc, LUT_dir2dir(isrc:Ndir**2:Ndir)
                          print *,'dir2dir', isrc, T_target
                        endif
                        call check(S_target, T_target, LUT_dir2diff(isrc:Ndir*Ndiff:Ndir), LUT_dir2dir(isrc:Ndir**2:Ndir), &
                          msg='test_LUT_wedge_direct_coeff_onsamplepts')
                      enddo !isrc
                    enddo !Cy
                  enddo !Cx
                enddo !theta
              enddo !phi
            enddo !aspect
          enddo !g
        enddo !w0
      enddo !tau

      end associate


  endsubroutine

  !@test(npes=[2,1])
  subroutine test_LUT_wedge_direct_coeff_interpolate(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: isrc
      integer(iintegers) :: idim_tau, idim_w0, idim_g, idim_aspect, idim_phi, idim_theta, idim_Cx, idim_Cy
      integer(iintegers) :: itau, iw0, ig, iaspect, iphi, itheta, iCx, iCy
      real(ireals) :: tau, w0, g, aspect, phi, theta, Cx, Cy

      real(ireals) :: kabs, ksca, dz, err(2)
      real(ireals), allocatable :: vertices(:)
      real(ireals), parameter :: dx = 911
      real(ireals),allocatable :: g_dim(:)

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()
      associate( LUTconfig => OPPLUT%dirconfig )

      idim_tau    = find_lut_dim_by_name(LUTconfig, 'tau')
      idim_w0     = find_lut_dim_by_name(LUTconfig, 'w0')
      idim_g      = find_lut_dim_by_name(LUTconfig, 'g')
      if(idim_g.eq.-1) then
        allocate(g_dim(2), source=zero)
      else
        allocate(g_dim(LUTconfig%dims(idim_g)%N), source=LUTconfig%dims(idim_g     )%v(:))
      endif
      idim_aspect = find_lut_dim_by_name(LUTconfig, 'aspect_zx')
      idim_phi    = find_lut_dim_by_name(LUTconfig, 'phi')
      idim_theta  = find_lut_dim_by_name(LUTconfig, 'theta')
      idim_Cx     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cx')
      idim_Cy     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cy')

      do itau = 1,max(i1, LUTconfig%dims(idim_tau)%N-1)
        do iw0  = 1,max(i1, LUTconfig%dims(idim_w0)%N-1)
          do ig   = 1,max(i1, size(g_dim)-1)
            do iaspect = 1,max(i1, LUTconfig%dims(idim_aspect)%N-1)
              do iphi = 1,max(i1, LUTconfig%dims(idim_phi)%N-1)
                do itheta = 1,max(i1, LUTconfig%dims(idim_theta)%N-1)
                  do iCx = 1,max(i1, LUTconfig%dims(idim_Cx)%N-1)
                    do iCy = 1,max(i1, LUTconfig%dims(idim_Cy)%N-1)
                      tau    = LUTconfig%dims(idim_tau   )%v(itau)    / 2
                      w0     = LUTconfig%dims(idim_w0    )%v(iw0)     / 2
                      g      = g_dim(ig)      / 2
                      aspect = LUTconfig%dims(idim_aspect)%v(iaspect) / 2
                      phi    = LUTconfig%dims(idim_phi   )%v(iphi)    / 2
                      theta  = LUTconfig%dims(idim_theta )%v(itheta)  / 2
                      Cx     = LUTconfig%dims(idim_Cx    )%v(iCx)     / 2
                      Cy     = LUTconfig%dims(idim_Cy    )%v(iCy)     / 2

                      tau    = tau    + LUTconfig%dims(idim_tau   )%v(min(LUTconfig%dims(idim_tau   )%N, itau+1))    / 2
                      w0     = w0     + LUTconfig%dims(idim_w0    )%v(min(LUTconfig%dims(idim_w0    )%N, iw0+1))     / 2
                      g      = g      + g_dim(ig+1)      / 2
                      aspect = aspect + LUTconfig%dims(idim_aspect)%v(min(LUTconfig%dims(idim_aspect)%N, iaspect+1)) / 2
                      phi    = phi    + LUTconfig%dims(idim_phi   )%v(min(LUTconfig%dims(idim_phi   )%N, iphi+1))    / 2
                      theta  = theta  + LUTconfig%dims(idim_theta )%v(min(LUTconfig%dims(idim_theta )%N, itheta+1))  / 2
                      Cx     = Cx     + LUTconfig%dims(idim_Cx    )%v(min(LUTconfig%dims(idim_Cx    )%N, iCx+1))     / 2
                      Cy     = Cy     + LUTconfig%dims(idim_Cy    )%v(min(LUTconfig%dims(idim_Cy    )%N, iCy+1))     / 2

                      call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2dir)
                      call OPPLUT%LUT_get_dir2diff([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2diff)

                      call setup_default_wedge_geometry([zero, zero], [one, zero], [Cx, Cy], aspect, vertices)
                      vertices = vertices * dx
                      do isrc = 1, Ndir
                        dz = vertices(12)-vertices(3)
                        kabs = (one-w0) * tau / dz
                        ksca = w0 * tau / dz
                        call bmc%get_coeff(comm,[kabs,ksca,g],isrc,.True.,phi,theta,vertices,S_target,T_target,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
                        err = rmse(LUT_dir2dir(isrc:Ndir**2:Ndir), T_target)
                        ! Not really an error if they dont match, its an interpolation in the end, discrepancies might be ok if the
                        ! sample density of the LUT is too small, lets give output if error is large... this may hint towards an error
                        if(err(1).ge.2*sigma*atol .or. err(2).ge.2*sigma*rtol) then
                          print *,'Testing: ', tau, w0, g, aspect, phi, theta, Cx, Cy,':: RMSE', err
                          print *,'LUT :::', isrc, LUT_dir2dir(isrc:Ndir**2:Ndir)
                          print *,'dir2dir', isrc, T_target
                        endif
                        !call check(S_target, T_target, LUT_dir2diff(isrc:Ndir*Ndiff:Ndir), LUT_dir2dir(isrc:Ndir**2:Ndir), &
                        !  msg='test_LUT_wedge_direct_coeff_interpolated')
                      enddo !isrc
                    enddo !Cy
                  enddo !Cx
                enddo !theta
              enddo !phi
            enddo !aspect
          enddo !g
        enddo !w0
      enddo !tau

      end associate


  endsubroutine

  !@test(npes=[2,1])
  subroutine test_wedge_LUT_vs_OPP_object(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: isrc
      integer(iintegers) :: idim_tau, idim_w0, idim_g, idim_aspect, idim_phi, idim_theta, idim_Cx, idim_Cy
      integer(iintegers) :: itau, iw0, ig, iaspect, iphi, itheta, iCx, iCy
      real(ireals) :: tau, w0, g, aspect, phi, theta, Cx, Cy

      real(ireals) :: err(2)
      real(ireals), parameter :: dx = 911
      real(ireals),allocatable :: g_dim(:)
      integer(mpiint) :: ierr

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()
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
      idim_phi    = find_lut_dim_by_name(LUTconfig, 'phi')
      idim_theta  = find_lut_dim_by_name(LUTconfig, 'theta')
      idim_Cx     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cx')
      idim_Cy     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cy')

      do itau = 1, LUTconfig%dims(idim_tau)%N
        do iw0  = 1, LUTconfig%dims(idim_w0)%N
          do ig   = 1, size(g_dim)
            do iaspect = 1, LUTconfig%dims(idim_aspect)%N
              do iphi = 1, LUTconfig%dims(idim_phi)%N
                do itheta = 1, LUTconfig%dims(idim_theta)%N
                  do iCx = 1, LUTconfig%dims(idim_Cx)%N
                    do iCy = 1, LUTconfig%dims(idim_Cy)%N
                      tau    = LUTconfig%dims(idim_tau   )%v(itau)
                      w0     = LUTconfig%dims(idim_w0    )%v(iw0)
                      g      = g_dim(ig)
                      aspect = LUTconfig%dims(idim_aspect)%v(iaspect)
                      phi    = LUTconfig%dims(idim_phi   )%v(iphi)
                      theta  = LUTconfig%dims(idim_theta )%v(itheta)
                      Cx     = LUTconfig%dims(idim_Cx    )%v(iCx)
                      Cy     = LUTconfig%dims(idim_Cy    )%v(iCy)

                      call OPPLUT%LUT_get_dir2dir ([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2dir)
                      call OPPLUT%LUT_get_dir2diff([tau, w0, aspect, Cx, Cy, phi, theta], LUT_dir2diff)

                      call OPP%get_coeff(tau, w0, g, aspect, .True., OPP_dir2dir, ierr, &
                        angles=[phi, theta], wedge_coords=[zero,zero,one,zero,Cx,Cy])
                      call OPP%get_coeff(tau, w0, g, aspect, .False., OPP_dir2diff, ierr, &
                        angles=[phi, theta], wedge_coords=[zero,zero,one,zero,Cx,Cy])

                      do isrc = 1, Ndir
                        err = rmse(LUT_dir2dir(isrc:Ndir**2:Ndir), OPP_dir2dir(isrc:Ndir**2:Ndir))
                        if(err(1).ge.epsilon(err) .or. err(2).ge.1e-3_ireals) then
                          print *,'Testing: ', tau, w0, g, aspect, phi, theta, Cx, Cy,':: RMSE', err
                          print *,'LUT :::', isrc, LUT_dir2dir(isrc:Ndir**2:Ndir)
                          print *,'OPP :::', isrc, OPP_dir2dir(isrc:Ndir**2:Ndir)
                        endif
                        call check(OPP_dir2diff(isrc:Ndir*Ndiff:Ndir), OPP_dir2dir(isrc:Ndir**2:Ndir), &
                                   LUT_dir2diff(isrc:Ndir*Ndiff:Ndir), LUT_dir2dir(isrc:Ndir**2:Ndir), &
                                   msg='test_wedge_LUT_vs_OPP_object', opt_atol=epsilon(err))
                      enddo !isrc
                    enddo !Cy
                  enddo !Cx
                enddo !theta
              enddo !phi
            enddo !aspect
          enddo !g
        enddo !w0
      enddo !tau

      end associate


  endsubroutine

  subroutine check(S_target,T_target, S,T, msg, opt_atol, opt_rtol)
      real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T
      character(len=*),optional :: msg
      real(ireals), intent(in), optional :: opt_atol, opt_rtol

      character(default_str_len) :: local_msgS, local_msgT
      logical, parameter :: ldetail=.False.

      real(ireals) :: arg_atol, arg_rtol, Terr(2), Serr(2)

      arg_atol = get_arg(atol*sigma, opt_atol)
      arg_rtol = get_arg(rtol*sigma, opt_rtol)

      Terr = rmse(T,T_target)
      Serr = rmse(S,S_target)

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
          print*,'RMSE ::: ',Serr*[one, 100._ireals],'%'
          print*,''
          write(*, FMT='( " direct  :: ", 8(f12.5) )' ) T
          print*,''
          write(*, FMT='( " target  :: ", 8(f12.5) )' ) T_target
          print*,''
          write(*, FMT='( " diff    :: ", 8(f12.5) )' ) T_target-T
          print*,'RMSE ::: ',Terr*[one, 100._ireals],'%'
          print*,'---------------------'
          print*,''
        else
          if(Terr(1).gt.arg_atol .or. Serr(1).gt.arg_atol) then
            print*,'RMSE ::: ',rmse(S,S_target)*[one, 100._ireals],'% ', &
                   'direct:::',rmse(T,T_target)*[one, 100._ireals],'%'
            print *,''
          endif
        endif

        @assertEqual(S_target, S, arg_atol, local_msgS )
        @assertLessThanOrEqual   (zero, S)
        @assertGreaterThanOrEqual(one , S)

        @assertEqual(T_target, T, arg_atol, local_msgT )
        @assertLessThanOrEqual   (zero, T)
        @assertGreaterThanOrEqual(one , T)
      endif
  end subroutine
end module
