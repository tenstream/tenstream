program main
      use m_boxmc, only : t_boxmc,t_boxmc_8_10,t_boxmc_1_2
      use mpi
      use m_data_parameters, only : mpiint,ireals,iintegers,one, init_mpi_data_parameters

      implicit none

      real(ireals) :: bg(3),S(10),T(8),phi0,theta0,dx,dy,dz

      integer(iintegers),parameter :: Niter=10
      real(ireals) :: Siter(Niter,size(S) ), Titer(Niter,size(T) )

      integer(iintegers) :: src=1,iter
      integer(mpiint) :: myid,ierr,numnodes

      real(ireals) :: tau,w,g

      type(t_boxmc_8_10) :: bmc_8_10


      print *,'Testing boxmc...'
      call MPI_Init(ierr)
      call MPI_Comm_Rank(MPI_COMM_WORLD, myid,ierr)
      call mpi_comm_size(MPI_COMM_WORLD, numnodes,ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)


      phi0 = 270.
      theta0=0.

      dx=70
      dy=70
      dz=40

       call bmc_8_10%init(MPI_COMM_WORLD)

!       bg = [1e-6, 0.13,.89 ]
!       call bmc_8_10%get_coeff(MPI_COMM_WORLD,bg,i1,S,T,.False.,.True.,phi0,theta0,dx,dy,dz)
!       if(myid.eq.0) write(*, FMT='( " direct ", 8(es15.5), "::",10(es15.5)  )' ) T,S

      if(.True.) then
        src=6
        do iter=1,Niter
          tau = 1e-6
!          w = dble(iter)/10._ireals-1e-3_ireals
          w = .9
          g = .0_ireals
          bg = [tau*(one-w)/dz, tau*w/dz, g ]
          bg = [1e-6,1e-6,0.]
          call bmc_8_10%get_coeff(MPI_COMM_WORLD,bg,src,S,T,.True.,.True.,phi0,theta0,dx,dy,dz)
          if(myid.le.1) write(*, FMT='( "iter ",I2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,T,S

          Siter(iter,:) = S
          Titer(iter,:) = T
        enddo
      endif

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
        call MPI_Finalize(ierr)

end program
