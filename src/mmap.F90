module m_mmap
  use iso_fortran_env, only: REAL32, INT8, INT16, INT32, INT64
  use iso_c_binding

  use m_data_parameters, only : iintegers, irealLUT, mpiint
  use m_helper_functions, only : CHKWARN, CHKERR, imp_bcast, itoa, resize_arr
  use m_netcdfIO, only: acquire_file_lock, release_file_lock

  use m_c_syscall_wrappers, only: c_sysconf, c_mmap, c_munmap, &
    c_open, c_close, c_lockf, &
    MAP_NORESERVE, MAP_PRIVATE, PROT_READ, SC_PAGESIZE, O_RDONLY

  implicit none

interface load_mmap_array
  module procedure load_mmap_array_2d
end interface
interface arr_to_binary_datafile
  module procedure arr_to_binary_datafile_2d
end interface

contains

  elemental function swap_endianness(realin) result (realout)
    real(REAL32), intent(in) :: realin
    real(REAL32)             :: realout
    integer                  :: i_element
    integer                  :: i_element_br
    i_element = transfer( realin, 0 )
    ! reverse order of 4 bytes in 32 bit integer space:
    call mvbits( i_element, 24, 8, i_element_br, 0  )
    call mvbits( i_element, 16, 8, i_element_br, 8  )
    call mvbits( i_element,  8, 8, i_element_br, 16 )
    call mvbits( i_element,  0, 8, i_element_br, 24 )
    realout = transfer( i_element_br, 0.0 )
  end function

  logical function get_system_endianness() result(littleendian)
    character(4) :: bytes !may not work on some processors

    bytes = transfer(1_int32,bytes)
    littleendian = .not. (ichar(bytes(4:4))==1)
  end function

  logical function binary_io_is_little_endian() result (littleendian)
    integer(int8) :: byte1, byte2
    integer(int16) :: twobytes

    integer :: funit

    twobytes = 1_int16

    open(newunit=funit, status='scratch', form='unformatted')
    write(funit) twobytes
    rewind(funit)
    read(funit) byte1, byte2
    close(funit)

    littleendian = .not. (byte1 == 0 .and. byte2 == 1)
  end function

  subroutine arr_to_binary_datafile_2d(arr, fname, ierr)
    real(irealLUT), dimension(:,:), intent(in) :: arr
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr

    integer :: funit, flock_unit
    logical :: lexists

    integer(c_size_t), allocatable :: header(:)
    integer(c_size_t) :: bytesize_header
    integer(c_long) :: c_pagesize

    integer(c_size_t) :: size_of_inp_arr, bytesize, dtype_size
    logical :: IO_is_little_endian, SYS_is_little_endian

    dtype_size=c_sizeof(arr(1,1))
    bytesize_header = c_sizeof(bytesize_header)
    size_of_inp_arr = size(arr, kind=c_size_t)
    bytesize = dtype_size * size_of_inp_arr

    c_pagesize = c_sysconf(SC_PAGESIZE)

    allocate(header(c_pagesize/bytesize_header))
    if(size(header).lt.3+size(shape(arr))) &
      call CHKERR(1_mpiint, 'Pagesize of this system is so small that we cant fit mmap header in first page element')
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
      call release_file_lock(flock_unit, ierr); call CHKERR(ierr)
      return
      open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='replace')
    else
      open(newunit=funit, file=trim(fname), form='unformatted', access='stream', status='new')
    endif

    ! Using TenStream in WRF needs to cope with binary data written with non native endianness.
    ! This comes from compiler option: -fconvert=big-endian
    IO_is_little_endian = binary_io_is_little_endian()
    SYS_is_little_endian = get_system_endianness()

    print *, 'ENDIANNESS little?', IO_is_little_endian, SYS_is_little_endian
    ! Note: I currently dont have a big endian system available to test it. So this may actually not work.
    if (IO_is_little_endian.eqv.SYS_is_little_endian) then
      write (unit=funit) header, arr
    else
      write (unit=funit) header, swap_endianness(arr)
    endif
    close(funit)
    call release_file_lock(flock_unit, ierr); call CHKERR(ierr)
  end subroutine

  subroutine binary_file_to_mmap(fname, mmap_c_ptr, dtype_size, arr_shape, ierr)
    character(len=*), intent(in) :: fname
    type(c_ptr), intent(out) :: mmap_c_ptr
    integer(c_size_t), intent(out) :: dtype_size
    integer(c_size_t), allocatable, intent(out), optional :: arr_shape(:)
    integer(mpiint), intent(out) :: ierr
    integer(c_size_t) :: offset

    integer :: funit, i
    integer(c_int) :: c_fd

    logical :: lexists

    integer(c_size_t) :: bytesize_data
    integer(c_long) :: c_pagesize
    integer(c_size_t), allocatable :: header(:)
    integer(c_size_t) :: bytesize_header

    bytesize_header = c_sizeof(bytesize_header)
    c_pagesize = c_sysconf(SC_PAGESIZE)
    offset = c_pagesize

    inquire(file=trim(fname), exist=lexists)
    if(.not.lexists) then
      ierr = 1
      call CHKWARN(1_mpiint, 'Tried to create a mmap from a file that does not exist: '//fname)
      return
    endif

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

    ierr = 0
  end subroutine

  subroutine load_mmap_array_2d(fname, arr, ierr)
    character(len=*), intent(in) :: fname
    real(irealLUT), pointer, intent(out) :: arr(:,:)
    integer(mpiint) :: ierr

    type(c_ptr) :: mmap_c_ptr
    integer(c_size_t), allocatable :: arrshape(:)
    integer(c_size_t) :: dtype_size

    call binary_file_to_mmap(fname, mmap_c_ptr, dtype_size, arrshape, ierr)
    if(ierr.ne.0) return

    if(dtype_size.ne.c_sizeof(arr(1,1))) &
      call CHKERR(1_mpiint, 'Could not correctly setup a memory map because'//new_line('')// &
      & 'size type of binary data '//itoa(dtype_size)//' does not match inp_arr_dtype '//itoa(c_sizeof(arr(1,1)))//new_line('')// &
      & '')

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
    integer(mpiint) :: myid
    integer(c_size_t) :: dtype_size ! size of the supplied irealLUT type
    integer(c_size_t) :: size_of_inp_arr, bytesize

    dtype_size = c_sizeof(1._irealLUT)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    fname_fpsuffix = trim(fname)//itoa(irealLUT)

    if(myid.eq.0.and.present(inp_arr)) then
      call arr_to_binary_datafile(inp_arr, fname_fpsuffix, ierr)
      if(ierr.ne.0) print *,'arr_to_mmap::binary file already exists, I did not overwrite it...'// &
                      ' YOU have to make sure that the file is as expected or delete it...'
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

    call load_mmap_array_2d(fname_fpsuffix, mmap_ptr, ierr)
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
