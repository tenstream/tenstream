module m_c_syscall_wrappers
  use iso_c_binding
  use m_data_parameters, only: mpiint
  use m_helper_functions, only: CHKERR, itoa, get_arg

  integer(c_int), bind(C, name="c_SC_PAGESIZE") :: SC_PAGESIZE
  integer(c_int), bind(C, name="cF_ULOCK") :: F_ULOCK
  integer(c_int), bind(C, name="cF_LOCK") :: F_LOCK

  integer(c_int), bind(C, name="cPROT_READ") :: PROT_READ
  integer(c_int), bind(C, name="cMAP_PRIVATE") :: MAP_PRIVATE
  integer(c_int), bind(C, name="cMAP_NORESERVE") :: MAP_NORESERVE

  integer(c_int), bind(C, name="cO_RDONLY") :: O_RDONLY
  integer(c_int), bind(C, name="cO_CREAT") :: O_CREAT
  integer(c_int), bind(C, name="cO_WRONLY") :: O_WRONLY
  integer(c_int), bind(C, name="cO_APPEND") :: O_APPEND
  integer(c_int), bind(C, name="cO_RDWR") :: O_RDWR
  integer(c_int), bind(C, name="cS_IRWXU") :: S_IRWXU
  integer(c_int), bind(C, name="cS_IRUSR") :: S_IRUSR
  integer(c_int), bind(C, name="cS_IWUSR") :: S_IWUSR
  integer(c_int), bind(C, name="cS_IROTH") :: S_IROTH
  integer(c_int), bind(C, name="default_user_wrmode") :: default_user_wrmode

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
    integer(c_int) function c_munmap(addr, length) bind(c, name='munmap')
      use iso_c_binding
      type(c_ptr), value :: addr
      integer(c_size_t), value :: length
    end function
  end interface

  interface
    integer(c_int) function c_open(pathname, flags) bind(c, name='open')
      use iso_c_binding
      character(kind=c_char, len=1) :: pathname(*)
      integer(c_int), value :: flags
    end function
  end interface
  interface
    integer(c_int) function c_creat(pathname, mode_t) bind(c, name='creat')
      use iso_c_binding
      character(kind=c_char, len=1) :: pathname(*)
      integer(c_int), value :: mode_t
    end function
  end interface

  interface
    integer(c_int) function c_close(fd) bind(c, name='close')
      use iso_c_binding
      integer(c_int), value :: fd
    end function
  end interface

  interface
    integer(c_long) function c_sysconf(name) bind(c, name='sysconf')
      use iso_c_binding
      integer(c_int), value :: name
    end function
  end interface

  interface
    integer(c_int) function c_lockf(fd, cmd, len) bind(c, name='lockf')
      use iso_c_binding
      integer(c_int), value :: fd, cmd
      integer(c_size_t), value :: len
    end function
  end interface

  interface
    function fsync(fd) bind(c, name="fsync")
      use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int) :: fsync
    end function fsync
  end interface

contains
  subroutine acquire_flock_lock(fname, ierr, l_throw_err)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: l_throw_err
    logical :: lexists
    integer(c_int) :: c_fd, cerr

    inquire (file=trim(fname), exist=lexists)

    if (.not. lexists) then
      ierr = 0
      return
      c_fd = c_creat(trim(fname)//c_null_char, default_user_wrmode)
    else
      c_fd = c_open(trim(fname)//c_null_char, O_RDONLY)
    end if
    if (c_fd .eq. -1_c_int) &
      call CHKERR(int(c_fd, mpiint), 'error opening file '//trim(fname)//' with c_open :: ')

    cerr = c_lockf(c_fd, F_LOCK, 0_c_size_t)
    if (get_arg(.false., l_throw_err)) then
      call CHKERR(cerr, 'Could not obtain File Lock'//itoa(cerr))
    end if
    ierr = cerr
  end subroutine

  subroutine release_flock_lock(fname, ierr, l_throw_err)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: l_throw_err
    logical :: lexists
    integer(c_int) :: c_fd, cerr

    inquire (file=trim(fname), exist=lexists)

    if (.not. lexists) then
      ierr = 0
    else
      c_fd = c_open(trim(fname)//c_null_char, O_RDONLY)
      if (c_fd .eq. -1_c_int) &
        call CHKERR(int(c_fd, mpiint), 'error opening file '//trim(fname)// &
                    ' with c_open')
      cerr = c_lockf(c_fd, F_ULOCK, 0_c_size_t)
      if (get_arg(.false., l_throw_err)) then
        call CHKERR(cerr, 'Could not release file Lock'//itoa(cerr))
      end if
      ierr = cerr
    end if
  end subroutine

end module
