module test_LUT_wedge_5_8
  use m_boxmc, only : t_boxmc_wedge_5_8
  use m_data_parameters, only : mpiint, ireals, iintegers, &
    one, zero, init_mpi_data_parameters, i1, default_str_len
  use m_optprop_LUT, only : t_optprop_LUT_wedge_5_8, find_lut_dim_by_name
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: rmse
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

  type(t_boxmc_wedge_5_8) :: bmc
  type(t_optprop_LUT_wedge_5_8) :: OPP

  integer(mpiint) :: myid,mpierr,numnodes,comm

  real(ireals),parameter :: atol=1e-4, rtol=1e-1
  real(ireals),parameter :: sigma = 3 ! normal test range for coefficients


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
      call read_commandline_options()

      call bmc%init(comm)
      call OPP%init(comm)

      call OPP%print_configs()

      S_target = zero
      T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      !call OPP%destroy()
      call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes=[1])
  subroutine test_LUT_direct_coeff(this)
      class (MpiTestMethod), intent(inout) :: this

      integer(iintegers) :: isrc
      integer(iintegers) :: idim_tau, idim_w0, idim_g, idim_aspect, idim_phi, idim_theta, idim_Cx, idim_Cy
      integer(iintegers) :: itau, iw0, ig, iaspect, iphi, itheta, iCx, iCy
      real(ireals) :: tau, w0, g, aspect, phi, theta, Cx, Cy

      real(ireals) :: kabs, ksca, dz, err(2)
      real(ireals), allocatable :: vertices(:)
      real(ireals), parameter :: dx = 911

      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()
      associate( LUTconfig => OPP%dirconfig )

      idim_tau    = find_lut_dim_by_name(LUTconfig, 'tau')
      idim_w0     = find_lut_dim_by_name(LUTconfig, 'w0')
      idim_g      = find_lut_dim_by_name(LUTconfig, 'g')
      idim_aspect = find_lut_dim_by_name(LUTconfig, 'aspect_zx')
      idim_phi    = find_lut_dim_by_name(LUTconfig, 'phi')
      idim_theta  = find_lut_dim_by_name(LUTconfig, 'theta')
      idim_Cx     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cx')
      idim_Cy     = find_lut_dim_by_name(LUTconfig, 'wedge_coord_Cy')

      do itau = 1, LUTconfig%dims(idim_tau)%N
        do iw0  = 1, LUTconfig%dims(idim_w0)%N
          do ig   = 1, LUTconfig%dims(idim_g)%N
            do iaspect = 1, LUTconfig%dims(idim_aspect)%N
              do iphi = 2, LUTconfig%dims(idim_phi)%N
                do itheta = 3, LUTconfig%dims(idim_theta)%N
                  do iCx = 3, LUTconfig%dims(idim_Cx)%N
                    do iCy = 1, LUTconfig%dims(idim_Cy)%N
                      tau    = LUTconfig%dims(idim_tau   )%v(itau)
                      w0     = LUTconfig%dims(idim_w0    )%v(iw0)
                      g      = LUTconfig%dims(idim_g     )%v(ig)
                      aspect = LUTconfig%dims(idim_aspect)%v(iaspect)
                      phi    = LUTconfig%dims(idim_phi   )%v(iphi)
                      theta  = LUTconfig%dims(idim_theta )%v(itheta)
                      Cx     = LUTconfig%dims(idim_Cx    )%v(iCx)
                      Cy     = LUTconfig%dims(idim_Cy    )%v(iCy)

                      call OPP%LUT_get_dir2dir ([tau, w0, g , aspect, Cx, Cy, phi, theta], LUT_dir2dir)
                      call OPP%LUT_get_dir2diff([tau, w0, g , aspect, Cx, Cy, phi, theta], LUT_dir2diff)

                      call setup_default_wedge_geometry([zero, zero], [one, zero], [Cx, Cy], aspect, vertices)
                      vertices = vertices * dx
                      do isrc = 1, Ndir
                        !print *,'Testing: ', tau, w0, g, aspect, phi, theta, Cx, Cy,'::', err
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
                        call check(S_target, T_target, LUT_dir2diff(isrc:Ndir*Ndiff:Ndir), LUT_dir2dir(isrc:Ndir**2:Ndir), msg='test_LUT_direct_coeffs')
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

  subroutine check(S_target,T_target, S,T, msg)
      real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

      character(len=*),optional :: msg
      character(default_str_len) :: local_msgS, local_msgT
      logical, parameter :: ldetail=.True.

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
          print*,'RMSE ::: ',rmse(S,S_target)*[one, 100._ireals],'%'
          print*,''
          write(*, FMT='( " direct  :: ", 8(f12.5) )' ) T
          print*,''
          write(*, FMT='( " target  :: ", 8(f12.5) )' ) T_target
          print*,''
          write(*, FMT='( " diff    :: ", 8(f12.5) )' ) T_target-T
          print*,'RMSE ::: ',rmse(T,T_target)*[one, 100._ireals],'%'
          print*,'---------------------'
          print*,''
        else
          print*,'RMSE ::: ',rmse(S,S_target)*[one, 100._ireals],'% ', &
                 'direct:::',rmse(T,T_target)*[one, 100._ireals],'%'
          print *,''
        endif

        @assertEqual(S_target, S, atol*sigma, local_msgS )
        @assertLessThanOrEqual   (zero, S)
        @assertGreaterThanOrEqual(one , S)

        @assertEqual(T_target, T, atol*sigma, local_msgT )
        @assertLessThanOrEqual   (zero, T)
        @assertGreaterThanOrEqual(one , T)
      endif
  end subroutine
end module
