program main
      use boxmc, only : t_boxmc,t_boxmc_8_10,t_boxmc_1_2
      use mpi
      use data_parameters, only : mpiint,ireals,iintegers,zero,one,i0,i1
      implicit none

      real(ireals) :: bg(3),S(10),T(8),phi0,theta0,dx,dy,dz

      integer(iintegers) :: src=1,itau,iter
      integer(mpiint) :: myid,ierr
      logical  :: direct

      real(ireals) :: tau,w,od_sca,g
      real(ireals),allocatable :: coeff(:)

      type(t_boxmc_8_10) :: bmc_8_10


      print *,'Testing boxmc...'
      call MPI_Init(ierr)
      call MPI_Comm_Rank(MPI_COMM_WORLD, myid,ierr)

      phi0 = 270.
      theta0=0.

      dx=2800.
      dy=2800.
      dz=1000.

       call bmc_8_10%init(MPI_COMM_WORLD)

       bg = [1e-6, 0.013,.89 ]
       call bmc_8_10%get_coeff(MPI_COMM_WORLD,bg,i1,S,T,.False.,.True.,phi0,theta0,dx,dy,dz)
       if(myid.eq.0) write(*, FMT='( " direct ", 8(es10.5), "::",10(es10.5)  )' ) T,S

!      if(.True.) then
!        src=1
!        do iter=1,10
!          tau = 20
!          w = dble(iter)/10._ireals-1e-3_ireals
!          g = .9_ireals
!          op_bg = [tau*(one-w)/dz, tau*w/dz, g ]
!          call bmc_get_coeff_8_10(MPI_COMM_WORLD,op_bg,src,S_out,Sdir_out,.True.,delta_scale,phi0,theta0,dx,dy,dz)
!          if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!        enddo
!
!        od_sca = tau*w
!        op_bg = [tau*(one-w)/dz, od_sca*(one-g)/dz, zero ]
!        call bmc_get_coeff_8_10(MPI_COMM_WORLD,op_bg,src,S_out,Sdir_out,.True.,delta_scale,phi0,theta0,dx,dy,dz)
!        if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!!      enddo
!!      if(myid.eq.0) write(*, FMT='( i2," direct ", 8(f10.5), "::",10(f10.5)  )' ) iter,Sdir_out,S_out
!      endif
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
