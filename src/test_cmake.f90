program main

  implicit none
  include 'mpif.h'

  integer :: myid,ierr

  call MPI_Init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myid, ierr)
  print *,'This is a cmake test',myid
  call MPI_FINALIZE(ierr)

end program
