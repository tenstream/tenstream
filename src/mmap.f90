module m_mmap
use iso_c_binding

use m_data_parameters, only : iintegers, ireals, mpiint
use m_helper_functions, only : CHKERR, imp_bcast

implicit none
! Definitions from mman-linux.h

integer, parameter :: &
  PROT_READWRITE = 3, &
  PROT_READ = 1,      &
  MAP_SHARED = 1

interface
  type(c_ptr) function c_mmap(addr, length, prot, &
      flags, fildes, off) result(result) bind(c, name='mmap')
    use iso_c_binding
    integer(c_int), value :: addr
    integer(c_size_t), value :: length
    integer(c_int), value :: prot
    integer(c_int), value :: flags
    integer(c_int), value :: fildes
    integer(c_size_t), value :: off
  end function
end interface
interface
  integer(c_int) function c_munmap(addr, length) bind(c,name='munmap')
    use iso_c_binding
    type(c_ptr), value :: addr
    integer(c_size_t), value :: length
  end function
end interface

interface arr_to_binary_datafile
  module procedure arr_to_binary_datafile_2d
end interface

contains

  subroutine arr_to_binary_datafile_2d(arr, fname, ierr)
    real(ireals), dimension(:,:), intent(in) :: arr
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr

    integer :: funit
    logical :: lexists

    inquire(file=trim(fname), exist=lexists)
    if(lexists) then
      ierr = 1
      return
    endif

    !print *,'Dumping to binary file ...'
    open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='new')
    write (unit=funit) arr
    close(funit)
    ierr = 0
    !print *,'Dumping to binary file ... finished'
  end subroutine

  subroutine binary_file_to_mmap(fname, bytesize, mmap_c_ptr)
    character(len=*), intent(in) :: fname
    integer(c_size_t), intent(in) :: bytesize
    type(c_ptr), intent(out) :: mmap_c_ptr
    integer(c_size_t),parameter :: offset=0

    integer :: funit
    logical :: lexists

    inquire(file=trim(fname), exist=lexists)
    if(.not.lexists) call CHKERR(1_mpiint, 'Tried to create a mmap from a file that does not exist: '//fname)

    open(newunit=funit, file=fname, form='unformatted', access='stream')
    mmap_c_ptr = c_mmap(0, bytesize, PROT_READ, MAP_SHARED, fnum(funit), offset)
    close(funit)
  end subroutine

  subroutine arr_to_mmap(comm, fname, inp_arr, mmap_ptr, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: fname
    real(ireals), dimension(:,:), intent(in) :: inp_arr
    real(ireals), pointer, intent(out) :: mmap_ptr(:,:)
    integer(mpiint), intent(out) :: ierr
    integer(iintegers), allocatable :: arrshape(:)

    type(c_ptr) :: mmap_c_ptr
    integer(c_size_t) :: bytesize
    integer(mpiint) :: myid
    integer(iintegers) :: size_of_inp_arr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      call arr_to_binary_datafile(inp_arr, fname, ierr)
      if(ierr.ne.0) print *,'arr_to_mmap::binary file already exists, I did not overwrite it... YOU have to make sure that the file is as expected or delete it...'
      size_of_inp_arr = int(sizeof(inp_arr), kind=iintegers)
      allocate(arrshape(size(shape(inp_arr))), source=shape(inp_arr))
    endif

    call imp_bcast(comm, arrshape, 0_mpiint)
    call imp_bcast(comm, size_of_inp_arr, 0_mpiint)

    bytesize = size_of_inp_arr

    call binary_file_to_mmap(fname, bytesize, mmap_c_ptr)

    call c_f_pointer(mmap_c_ptr, mmap_ptr, arrshape)
  end subroutine

end module
