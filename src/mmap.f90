module m_mmap
  use iso_fortran_env, only: INT64
  use iso_c_binding

  use m_data_parameters, only : iintegers, irealLUT, mpiint
  use m_helper_functions, only : CHKERR, imp_bcast, itoa, resize_arr
  use m_netcdfIO, only: acquire_file_lock, release_file_lock

  implicit none
  !Posix Standard Definition for sysconf
  integer, parameter :: &
    SC_PAGESIZE = 30

  ! Definitions from mman-linux.h
  integer, parameter :: &
    PROT_READWRITE = 3, &
    PROT_READ = 1,      &
    MAP_SHARED = 1,     &
    MAP_PRIVATE = 2,    &
    MAP_NORESERVE = 16384, & ! X'04000'
    O_RDONLY = 0

interface
  type(c_ptr) function c_mmap(addr, length, prot, &
      flags, filedes, off) result(result) bind(c, name='mmap')
    use iso_c_binding
    integer(c_int), value :: addr
    integer(c_size_t), value :: length
    integer(c_int), value :: prot
    integer(c_int), value :: flags
    integer(c_int), value :: filedes
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
interface
  integer(c_int) function c_open(pathname, flags) bind(c,name='open')
    use iso_c_binding
    character(kind=c_char, len=1) :: pathname(*)
    integer(c_int), value :: flags
  end function
end interface
interface
  integer(c_int) function c_close(fd) bind(c,name='close')
    use iso_c_binding
    integer(c_int), value :: fd
  end function
end interface
interface
  integer(c_long) function c_sysconf(name) bind(c,name='sysconf')
    use iso_c_binding
    integer(c_int), value :: name
  end function
end interface


interface load_mmap_array
  module procedure load_mmap_array_2d
end interface
interface arr_to_binary_datafile
  module procedure arr_to_binary_datafile_2d
end interface

contains

  subroutine arr_to_binary_datafile_2d(arr, fname, ierr)
    real(irealLUT), dimension(:,:), intent(in) :: arr
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr

    integer :: funit, flock_unit
    logical :: lexists

    integer(c_size_t), allocatable :: header(:)
    integer(c_size_t), parameter :: bytesize_header=c_sizeof(bytesize_header)
    integer(c_long) :: c_pagesize

    integer(c_size_t) :: size_of_inp_arr, bytesize
    integer(c_size_t), parameter :: dtype_size=c_sizeof(arr(1,1))

    size_of_inp_arr = size(arr, kind=c_size_t)
    bytesize = dtype_size * size_of_inp_arr

    c_pagesize = c_sysconf(SC_PAGESIZE)

    allocate(header(c_pagesize/bytesize_header))
    header(:) = 0
    header(1) = dtype_size
    header(2) = size_of_inp_arr
    header(3) = bytesize
    header(4) = size(arr,dim=1)
    header(5) = size(arr,dim=2)

    call acquire_file_lock(fname, flock_unit, ierr); call CHKERR(ierr)
    ierr = 0
    inquire(file=trim(fname), exist=lexists)
    if(lexists) then
      open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='replace')
    else
      open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='new')
    endif
    write (unit=funit) header, arr
    close(funit)
    call release_file_lock(flock_unit, ierr); call CHKERR(ierr)
  end subroutine

  subroutine binary_file_to_mmap(fname, mmap_c_ptr, dtype_size, arr_shape)
    character(len=*), intent(in) :: fname
    type(c_ptr), intent(out) :: mmap_c_ptr
    integer(c_size_t), intent(out) :: dtype_size
    integer(c_size_t), allocatable, intent(out), optional :: arr_shape(:)
    integer(c_size_t) :: offset

    integer :: funit, i
    integer(c_int) :: c_fd

    logical :: lexists

    integer(c_size_t) :: bytesize_data
    integer(c_long) :: c_pagesize
    integer(c_size_t), allocatable :: header(:)
    integer(c_size_t), parameter :: bytesize_header=c_sizeof(bytesize_header)

    c_pagesize = c_sysconf(SC_PAGESIZE)
    offset = c_pagesize

    inquire(file=trim(fname), exist=lexists)
    if(.not.lexists) call CHKERR(1_mpiint, 'Tried to create a mmap from a file that does not exist: '//fname)

    allocate(header(c_pagesize/bytesize_header))

    open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='old')
    read (unit=funit) header
    close(funit)
    dtype_size = header(1)
    bytesize_data = header(3)

    if(present(arr_shape)) then
      allocate(arr_shape(count(header(4:size(header)).ne.0)))
      do i = 1, size(arr_shape)
        arr_shape(i) = header(3+i)
      enddo
    endif

    c_fd = c_open(trim(fname)//C_NULL_CHAR, O_RDONLY)
    if(c_fd.le.0) call CHKERR(c_fd, 'Could not open mmap file')
    mmap_c_ptr = c_mmap(0_c_int, bytesize_data, PROT_READ, IOR(MAP_PRIVATE, MAP_NORESERVE), c_fd, offset)
    c_fd = c_close(c_fd)
    if(c_fd.le.0) call CHKERR(c_fd, 'Could not close file descriptor to mmap file')
  end subroutine

  subroutine load_mmap_array_2d(fname, arr)
    character(len=*), intent(in) :: fname
    real(irealLUT), pointer, intent(out) :: arr(:,:)

    type(c_ptr) :: mmap_c_ptr
    integer(c_size_t), allocatable :: arrshape(:)
    integer(c_size_t) :: dtype_size

    call binary_file_to_mmap(fname, mmap_c_ptr, dtype_size, arrshape)

    if(dtype_size.ne.c_sizeof(arr(1,1))) &
      call CHKERR(1_mpiint, 'size type of binary data '//itoa(dtype_size)//'does not match inp_arr_dtype'//itoa(c_sizeof(arr(1,1))))

    if(size(arrshape).ne.2) &
      call CHKERR(size(arrshape), 'mmap in '//trim(fname)//' has '//itoa(size(arrshape))//' dimensions... expected 2')

    call c_f_pointer(mmap_c_ptr, arr, arrshape)
  end subroutine

  subroutine arr_to_mmap(comm, fname, mmap_ptr, ierr, inp_arr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: fname
    real(irealLUT), dimension(:,:), intent(in), optional :: inp_arr
    real(irealLUT), pointer, intent(out) :: mmap_ptr(:,:)
    integer(mpiint), intent(out) :: ierr

    character(len=len_trim(fname)+2) :: fname_fpsuffix
    integer(iintegers), allocatable :: arrshape(:)
    type(c_ptr) :: mmap_c_ptr
    integer(mpiint) :: myid
    integer(c_size_t), parameter :: dtype_size=c_sizeof(1._irealLUT) ! size of the supplied irealLUT type
    integer(c_size_t) :: size_of_inp_arr, dtype_size_mmap_data, bytesize

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    fname_fpsuffix = trim(fname)//itoa(irealLUT)

    if(myid.eq.0) then
      if(.not.present(inp_arr)) call CHKERR(1_mpiint, 'rank 0 has to provide an input array!')
      call arr_to_binary_datafile(inp_arr, fname_fpsuffix, ierr)
      if(ierr.ne.0) print *,'arr_to_mmap::binary file already exists, I did not overwrite it... YOU have to make sure that the file is as expected or delete it...'
      size_of_inp_arr = size(inp_arr, kind=c_size_t)
      bytesize = dtype_size * size_of_inp_arr
      if(size_of_inp_arr.ge.huge(size_of_inp_arr)/dtype_size) then
        print *,'dtype_size', dtype_size
        print *,'size_of_inp_arr', size_of_inp_arr
        print *,'bytesize', bytesize
        call CHKERR(1_mpiint, 'size_of_inp_arr too large for c_size_t')
      endif
      allocate(arrshape(size(shape(inp_arr))), source=shape(inp_arr, kind=iintegers))
    endif
    call mpi_barrier(comm, ierr)

    call imp_bcast(comm, arrshape, 0_mpiint)
    call imp_bcast(comm, bytesize, 0_mpiint)

    if(bytesize.le.0_c_size_t) then
      if(myid.eq.0) print *,'huge(c_size_t)', huge(size_of_inp_arr), 'shape inp_arr', shape(inp_arr)
      call mpi_barrier(comm, ierr)
      call CHKERR(1_mpiint, 'bytesize of mmap is wrong!'//itoa(int(bytesize,iintegers)))
    endif

    call binary_file_to_mmap(fname_fpsuffix, mmap_c_ptr, dtype_size_mmap_data)
    if(dtype_size.ne.dtype_size_mmap_data) &
      call CHKERR(1_mpiint, 'size type of binary data '//itoa(dtype_size)// &
        ' does not match inp_arr_dtype '//itoa(dtype_size_mmap_data))

    call c_f_pointer(mmap_c_ptr, mmap_ptr, arrshape)
  end subroutine

  subroutine munmap_mmap_ptr(mmap_ptr, ierr)
    real(irealLUT), pointer, intent(inout) :: mmap_ptr(:,:)

    type(c_ptr) :: mmap_c_ptr
    integer(c_int) :: cerr
    integer(c_size_t) :: bytesize

    integer(mpiint) :: ierr

    bytesize = int(sizeof(mmap_ptr), kind=c_size_t)
    mmap_c_ptr = c_loc(mmap_ptr(1,1))

    mmap_ptr=>NULL()

    cerr = c_munmap(mmap_c_ptr, bytesize); ierr = cerr
    call CHKERR(ierr, 'Error Unmapping the memory map')
  end subroutine

end module
