program main
      use m_boxmc, only : t_boxmc,t_boxmc_8_10,t_boxmc_1_2
      use m_optprop, only : t_optprop_1_2,t_optprop_8_10
      use mpi
      use m_data_parameters, only : mpiint,ireals,iintegers,one, init_mpi_data_parameters,i1
      use m_eddington, only : eddington_coeff_fab
      use m_helper_functions, only: deg2rad,delta_scale

      implicit none

      real(ireals) :: bg(3),delta_bg(3),delta1_bg(3),S(10),T(8),phi0,theta0,dx,dy,dz

      integer(iintegers),parameter :: Niter=10
      real(ireals) :: Siter(Niter,size(S) ), Titer(Niter,size(T) )
      real(ireals),allocatable,dimension(:) :: dir2dir_coeff,dir2diff_coeff,diff2diff_coeff

      integer(iintegers) :: src=1,iter
      integer(mpiint) :: myid,ierr,numnodes

      real(ireals) :: tau,w,g
      real(ireals) :: a11,a12,a13,a23,a33

      type(t_boxmc_8_10) :: bmc_8_10

      type(t_optprop_8_10) OPP_8_10
      type(t_optprop_1_2) OPP_1_2


      phi0 = 0.
      theta0=0.

      dx=70
      dy=70
      dz=40

      print *,'Testing boxmc...'
      call MPI_Init(ierr)
      call MPI_Comm_Rank(MPI_COMM_WORLD, myid,ierr)
      call mpi_comm_size(MPI_COMM_WORLD, numnodes,ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)

      call bmc_8_10%init(MPI_COMM_WORLD)
      call OPP_8_10%init(dx,dy,[phi0],[theta0],MPI_COMM_WORLD)


      allocate( dir2dir_coeff(OPP_8_10%OPP_LUT%dir_streams**2) )
      allocate( dir2diff_coeff(OPP_8_10%OPP_LUT%dir_streams*OPP_8_10%OPP_LUT%diff_streams) )
      allocate( diff2diff_coeff(OPP_8_10%OPP_LUT%diff_streams**2) )


      bg = [1e-6, 1e-3, .8 ]
      delta_bg = bg
      delta1_bg = bg
      call delta_scale(delta1_bg(1),delta1_bg(2),delta1_bg(3),delta1_bg(3)**1.5 )
      call delta_scale(delta_bg(1),delta_bg(2),delta_bg(3) )
      print *,'unscaled optprop',bg
      print *,'  scaled optprop',delta_bg
      print *,'1.5 aled optprop',delta1_bg

      if(myid.eq.0) then
        call bmc_8_10%get_coeff(MPI_COMM_WORLD,bg,i1,S,T,.True.,phi0,theta0,dx,dy,dz)
        write(*, FMT='( " direct normal     ", 8(es10.3), "  ::  ",10(es10.3)  )' ) T,S

        call bmc_8_10%get_coeff(MPI_COMM_WORLD,delta_bg,i1,S,T,.True.,phi0,theta0,dx,dy,dz)
        write(*, FMT='( " direct deltascaled", 8(es10.3), "  ::  ",10(es10.3)  )' ) T,S

        call bmc_8_10%get_coeff(MPI_COMM_WORLD,delta1_bg,i1,S,T,.True.,phi0,theta0,dx,dy,dz)
        write(*, FMT='( " direct 1.5  scaled", 8(es10.3), "  ::  ",10(es10.3)  )' ) T,S

        call OPP_8_10%get_coeff(dz,delta_bg(1),delta_bg(2),delta_bg(3),.True.,dir2dir_coeff,[phi0,theta0])
        call OPP_8_10%get_coeff(dz,delta_bg(1),delta_bg(2),delta_bg(3),.False.,dir2diff_coeff,[phi0,theta0])
        write(*, FMT='( " direct LUT        ", 8(es10.3), "  ::  ",10(es10.3)  )' ) dir2dir_coeff(1:OPP_8_10%OPP_LUT%dir_streams),dir2diff_coeff(1:OPP_8_10%OPP_LUT%diff_streams)

        tau = (delta_bg(1)+delta_bg(2))*dz
        w   = delta_bg(2)/(delta_bg(1)+delta_bg(2))
        g   = delta_bg(3)
        call eddington_coeff_fab ( tau , w, g, cos(deg2rad(theta0)), &
            a11,          &
            a12,          &
            a13,          &
            a23,          &
            a33 )

        write(*, FMT='( "directeddington    ", 1(es10.3),"                                                                       ::   ",2(es10.3) )' ) a33,a13,a23
      endif

      if(.False.) then
        src=6
        do iter=1,Niter
          tau = 1e-6
!          w = dble(iter)/10._ireals-1e-3_ireals
          w = .9
          g = .0_ireals
          bg = [tau*(one-w)/dz, tau*w/dz, g ]
          bg = [1e-6,1e-6,0.]
          call bmc_8_10%get_coeff(MPI_COMM_WORLD,bg,src,S,T,.True.,phi0,theta0,dx,dy,dz)
          if(myid.le.1) write(*, FMT='( "iter ",I2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,T,S

          Siter(iter,:) = S
          Titer(iter,:) = T
        enddo

      if(myid.eq.0) then
        print *,''
        print *,''

        write(*, FMT='( "    mean    :  ",8(f10.5) ,"::",10(f10.5) )' ) sum(Titer,dim=1 )/Niter  , sum(Siter,dim=1 )/Niter
        write(*, FMT='( "    stddev  :  ",8(f10.5) ,"::",10(f10.5) )' ) sqrt(one*Niter*numnodes)* (sqrt(sum(Titer**2,dim=1 )/Niter - (sum(Titer,dim=1 )/Niter)**2) ),  &
                                                                        sqrt(one*Niter*numnodes)* (sqrt(sum(Siter**2,dim=1 )/Niter - (sum(Siter,dim=1 )/Niter)**2) )
        write(*, FMT='( "rel stddev  :  ",8(f10.5) ,"::",10(f10.5) )' ) sqrt(one*Niter*numnodes)* (sqrt(sum(Titer**2,dim=1 )/Niter - (sum(Titer,dim=1 )/Niter)**2) / (sum(Titer,dim=1 )/Niter)), &
                                                                        sqrt(one*Niter*numnodes)* (sqrt(sum(Siter**2,dim=1 )/Niter - (sum(Siter,dim=1 )/Niter)**2) / (sum(Siter,dim=1 )/Niter))
        print *,''

      endif
      call mpi_barrier(mpi_comm_world,ierr)
      if(myid.eq.1) then
        print *,''
        print *,''

        write(*, FMT='( "    mean    :  ",8(f10.5) ,"::",10(f10.5) )' ) sum(Titer,dim=1 )/Niter  , sum(Siter,dim=1 )/Niter
        write(*, FMT='( "    stddev  :  ",8(f10.5) ,"::",10(f10.5) )' ) sqrt(one*Niter)*(sqrt(sum(Titer**2,dim=1 )/Niter - (sum(Titer,dim=1 )/Niter)**2) )  , &
                                                                        sqrt(one*Niter)*(sqrt(sum(Siter**2,dim=1 )/Niter - (sum(Siter,dim=1 )/Niter)**2) )
        write(*, FMT='( "rel stddev  :  ",8(f10.5) ,"::",10(f10.5) )' ) sqrt(one*Niter)*(sqrt(sum(Titer**2,dim=1 )/Niter - (sum(Titer,dim=1 )/Niter)**2) / (sum(Titer,dim=1 )/Niter)), &
                                                                        sqrt(one*Niter)*(sqrt(sum(Siter**2,dim=1 )/Niter - (sum(Siter,dim=1 )/Niter)**2) / (sum(Siter,dim=1 )/Niter))
        print *,''

      endif

      endif

      !        od_sca = tau*w
!        op_bg = [tau*(one-w)/dz, od_sca*(one-g)/dz, zero ]
!        call bmc_get_coeff_8_10(MPI_COMM_WORLD,op_bg,src,S_out,Sdir_out,.True.,delta_scale,phi0,theta0,dx,dy,dz)
!        if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!      enddo
!      if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!
!      if(.False.) then
!        print *,'Testing optprop'
!        call optprop_8_10_init(dx,dy,[phi0],[theta0],MPI_COMM_WORLD)
!
!        optprop_8_10_debug = .True.
!        
!        tau = 1e-0_ireals/dz
!        w = .9_ireals
!        g = .9_ireals
!        allocate(coeff(dir_streams*diff_streams) )
!        call optprop_8_10_lookup_coeff(dz,tau ,w,g,direct,coeff,[zero,zero])
!
!      endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Finalize(ierr)

end program
