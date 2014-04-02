program main
      use tenstream_optprop_8_10
      use boxmc_8_10
      use mpi
      use data_parameters, only : mpiint,ireals,iintegers,zero,one
      use boxmc_parameters_8_10,only : dir_streams,diff_streams,delta_scale
      implicit none

      double precision :: op_bg(3),S_out(diff_streams),Sdir_out(dir_streams),phi0,theta0,dx,dy,dz

      integer(iintegers) :: src=1,itau,iter
      integer(mpiint) :: myid,ierr
      logical  :: direct

      double precision :: tau,w,od_sca,g
      double precision,allocatable :: coeff(:)


      print *,'Testing boxmc...'
      call MPI_Init(ierr)
      call MPI_Comm_Rank(MPI_COMM_WORLD, myid,ierr)

      phi0 = 270.
      theta0=0.

      dx=2800.
      dy=2800.
      dz=1000.

      if(.True.) then
        src=1
        do iter=1,10
          tau = 20
          w = dble(iter)/10._ireals-1e-3_ireals
          g = .9_ireals
          op_bg = [tau*(one-w)/dz, tau*w/dz, g ]
          call bmc_get_coeff(MPI_COMM_WORLD,op_bg,src,S_out,Sdir_out,.True.,delta_scale,phi0,theta0,dx,dy,dz)
          if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
        enddo

        od_sca = tau*w
        op_bg = [tau*(one-w)/dz, od_sca*(one-g)/dz, zero ]
        call bmc_get_coeff(MPI_COMM_WORLD,op_bg,src,S_out,Sdir_out,.True.,delta_scale,phi0,theta0,dx,dy,dz)
        if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!      enddo
!      if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
      endif

      if(.True.) then
        print *,'Testing optprop'
        call init_optprop(dx,dy,[phi0],[theta0],MPI_COMM_WORLD)

        optprop_debug = .True.
        
        tau = 1e-0_ireals/dz
        w = .9_ireals
        g = .9_ireals
        allocate(coeff(dir_streams*diff_streams) )
        call optprop_lookup_coeff(dz,tau ,w,g,direct,coeff,[zero,zero])

      endif
        call MPI_Finalize(ierr)

end program
