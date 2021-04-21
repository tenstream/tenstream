module test_geometric_coeffs
  use m_boxmc, only : t_boxmc, t_boxmc_3_10
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, irealLUT, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_geometric_coeffs, only : dir2dir3_geometric_coeffs
  use m_helper_functions, only : spherical_2_cartesian, cstr, toStr

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(3), S_target(10), T_target(3)
  real(ireals) :: S_tol(10),T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_10) :: bmc_3_10

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireal_dp),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-2
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()


    call init_mpi_data_parameters(comm)

    call bmc_3_10%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if(myid.eq.0) print *,numnodes,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 50
    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown


  @test(npes =[1])
  subroutine test_boxmc_dir2dir3_geometric_coeffs_vs_online_monte_carlo(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), allocatable :: verts(:), verts_dtd(:)
    real( ireals), parameter :: dx=1, dy= dx, dz=1 * dx
    real(irealLUT) :: v(9), v_mc(9)
    real(ireals) :: sundir(3)
    integer(iintegers) :: itheta, iphi

    bg = [real(ireal_dp) ::  -log(5e-1)/dz, 0e-0/dz, 1./2 ]
    S_target = zero
    iphi=30
    itheta=5

    do iphi=181,360,30
      do itheta=10,50,20
        phi = real(iphi, ireals)
        theta = real(itheta, ireals)

        call setup_default_unit_cube_geometry(dx, dy, dz, verts)
        verts_dtd = verts
        verts_dtd([9,12,21,24]) = verts_dtd([9,12,21,24]) + dz / 4

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v(src:3**2:3) = real(T, irealLUT)
        enddo

        sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals))

        print *, cstr('regular not corrected', 'red')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        call dir2dir3_geometric_coeffs(verts_dtd, sundir, real(bg(1) + bg(2), ireals), v)

        print *, cstr('regular corrected', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
      enddo
    enddo

  end subroutine

  !@test(npes =[1])
  subroutine test_boxmc_dir2dir3_geometric_coeffs(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), dimension(24) :: verts_dtd
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), v_mc(9), v_mc_dir2diff_reg(30), v_mc_dir2diff_dst(30)
    real(ireals) :: sundir(3), extinction_coeff
    integer(iintegers) :: itheta, iphi
    real( ireals), parameter :: dx=1, dy= dx, dz=dx

    !                        k_abs (=1e-6), k_scatter (=k_abs * 100), 1./2 sym parameter
    bg = [real(ireal_dp) :: 1e-6, 1, 0.85]

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    verts_dtd = verts


    verts_dtd([9,12,21,24]) = verts_dtd([9,12,21,24]) + dz / 2

    S_target = zero

    iphi=181
    itheta=40
    !do iphi=0,360,30
    !  do itheta=10,50,20
        phi = real(iphi, ireals)
        theta = real(itheta, ireals)

        sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals))

        call dir2dir3_geometric_coeffs(verts_dtd, sundir, real(bg(1) + bg(2), irealLUT), v)

        print *, cstr('regular corrected', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc_dir2diff_reg(src:30:3) = real(S, irealLUT)
        enddo

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
          v_mc_dir2diff_dst(src:30:3) = real(S, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        print *, cstr('dir2diff', 'red')
        print *, 'srz z', v_mc_dir2diff_reg(1:30:3)
        print *, 'srz x', v_mc_dir2diff_reg(2:30:3)
        print *, 'srz y', v_mc_dir2diff_reg(3:30:3)

        print *, cstr('dir2diff', 'red')
        print *, 'srz z', v_mc_dir2diff_dst(1:30:3)
        print *, 'srz x', v_mc_dir2diff_dst(2:30:3)
        print *, 'srz y', v_mc_dir2diff_dst(3:30:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
!        @assertEqual(v_mc([2,3,4]), v([2,3,4]), max(maxval(v_mc([2,3,4]))*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
    !  enddo
    !enddo
    end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

    character(len=*),optional :: msg
    character(default_str_len) :: local_msgS, local_msgT

    if(myid.eq.0) then
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
      write(*, FMT='( " diffuse ::  :: ",10(es12.5) )' ) S
      write(*, FMT='( " target  ::  :: ",10(es12.5) )' ) S_target
      write(*, FMT='( " diff    ::  :: ",10(es12.5) )' ) S_target-S
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
      print*,'---------------------'
      print*,''

      @assertEqual(S_target, S, real(atol*sigma, ireals), local_msgS )
      @assertLessThanOrEqual   (zero, S)
      @assertGreaterThanOrEqual(one , S)

      @assertEqual(T_target, T, real(atol*sigma, ireals), local_msgT )
      @assertLessThanOrEqual   (zero, T)
      @assertGreaterThanOrEqual(one , T)
    endif
  end subroutine

end module
