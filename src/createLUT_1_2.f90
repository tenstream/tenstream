program main
      use data_parameters, only: mpiint, ireals
      use mpi
      use optprop_LUT, only : t_optprop_LUT_1_2

      integer(mpiint) :: myid,comm

      character(len=32) :: arg
      real(ireals) :: dx
      real(ireals) :: azis(2),szas(5)

      type(t_optprop_LUT_1_2) :: OPP

      call mpi_init(ierr)
      comm = MPI_COMM_WORLD
      call mpi_comm_rank(comm,myid,ierr)

      azis = [90,0]
      szas = [0,20,40,60,80]

      call get_command_argument(1, arg)
      if(len_trim(arg) == 0) call exit
      read (arg,*) dx

      print *,'calculating coeffs for dx',dx
      call OPP%init(dx,dx,azis,szas,comm)

      call mpi_finalize(ierr)
end program
