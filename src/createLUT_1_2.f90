program main
      use m_data_parameters, only: mpiint, ireals, init_mpi_data_parameters
      use mpi
      use m_optprop_LUT, only : t_optprop_LUT_1_2

      integer(mpiint) :: myid,comm

      character(len=32) :: arg
      real(ireals) :: dx,user_sza
      real(ireals) :: azis(2),szas(5)

      type(t_optprop_LUT_1_2) :: OPP

      call mpi_init(ierr)
      comm = MPI_COMM_WORLD
      call mpi_comm_rank(comm,myid,ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)

      azis = [0,90]
      szas = [0,20,40,60,80]

      call get_command_argument(1, arg)
      if(len_trim(arg) == 0) call exit
      read (arg,*) dx

      call get_command_argument(2, arg)
      if(len_trim(arg) .gt. 0) then
        read (arg,*) user_sza
        szas=user_sza
      endif

      print *,'calculating coeffs for dx',dx
      call OPP%init(dx,dx,azis,szas,comm)
      print *,'loaded 1_2 coeffs for dx',dx,'szas',szas,'azis',azis


      call mpi_finalize(ierr)
end program
