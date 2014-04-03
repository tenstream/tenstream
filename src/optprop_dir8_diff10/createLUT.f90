program main
      use mpi
      use tenstream_optprop_LUT_8_10

      integer myid

      character(len=32) :: arg
      double precision :: dx
      double precision :: azis(2),szas(5)

      azis = [0,90]
      szas = [0,20,40,60,80]

      call get_command_argument(1, arg)
      if(len_trim(arg) == 0) call exit
      read (arg,*) dx

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
      print *,'calculating coeffs for dx',dx
      call init_LUT(dx,dx,azis,szas,MPI_COMM_WORLD)

      call mpi_finalize(ierr)
end program
