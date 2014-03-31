program main
      use tenstream_optprop_1_2 ,only : optprop_1_2_init, optprop_1_2_lookup_coeff,optprop_1_2_debug
      use boxmc_1_2,only : bmc_get_coeff
      use mpi, only : MPI_Init,MPI_Comm_Rank,MPI_Finalize,MPI_COMM_WORLD
      use data_parameters, only : ireals,iintegers,i1,zero,one
      use boxmc_parameters_1_2, only : dir_streams, diff_streams, delta_scale
      use eddington, only : rodents
      implicit none

      double precision :: op_bg(3),S_out(diff_streams),Sdir_out(dir_streams),phi0,theta0,dx,dy,dz

      integer(iintegers)  :: src=i1,itau,iter
      integer :: ierr,myid,comm
      logical :: direct

      double precision :: tau,w,od_sca,g
      double precision,allocatable :: coeff(:)
      double precision :: twostr(5)


      print *,'Testing two stream boxmc...'
      call MPI_Init(ierr)
      comm = MPI_COMM_WORLD
      call MPI_Comm_Rank(comm, myid,ierr)

      phi0 = 270.
      theta0=0.

      dx=28000.
      dy=28000.
      dz=100.

      if(.True.) then
        direct=.True.
        src=i1
        do iter=1,10
          tau = iter*one
          w = .5_ireals !dble(iter)/10._ireals-1e-3_ireals
          g = .9_ireals
          op_bg = [tau*(one-w)/dz, tau*w/dz, g ]

          print *,'Testing direct Solar coefficients, no scattering: optprop: ',op_bg

          call bmc_get_coeff(comm,op_bg,src,S_out,Sdir_out,direct,.False.,phi0,theta0,dx,dy,dz)

          if(myid.eq.0) write(*, FMT='( i2," direct ", 1(f10.5), "::",2(f10.5)  )' ) iter,Sdir_out,S_out

          call rodents( tau, w, g, cos(theta0), twostr)

          if(myid.eq.0) write(*, FMT='( i2," Rodents ", 1(f10.5), "::",2(f10.5)  )' ) iter,twostr(5),twostr([3,4])

          print *,'lambert Beer absorption',one-exp(-tau*(one-w)),'scatter',one-exp(-tau*w),'Transmission',exp(-tau)
          print *,''
        enddo

      endif

      if(.False.) then
        print *,'Testing optprop'
        call optprop_1_2_init(dx,dy,comm)

        optprop_1_2_debug = .True.
        
        tau = 1e-0_ireals/dz
        w = .9_ireals
        g = .9_ireals
        allocate(coeff(dir_streams*diff_streams) )
        call optprop_1_2_lookup_coeff(dz,tau ,w,g,direct,coeff,[zero,zero])

      endif
        call MPI_Finalize(ierr)

end program
