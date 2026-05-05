!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

!> Standalone key-value options database with no external dependencies.
!> Parses -key [value] pairs from command-line arguments and options files.
!> Mirrors the PETSc options database API as used by TenStream:
!>   opts_init        → PetscInitialize + PetscOptionsInsertFile
!>   opts_has         → PetscOptionsHasName
!>   opts_get_raw     → basis for get_petsc_opt wrappers
!>   opts_get_logical → PetscOptionsGetBool
!>   opts_get_integer → PetscOptionsGetInt
!>   opts_get_real    → PetscOptionsGetReal
!>   opts_get_string  → PetscOptionsGetString
!>   opts_get_*_array → PetscOptionsGet*Array
module m_options_database
  use iso_fortran_env, only: int32, real64, error_unit
  implicit none
  private

  public :: opts_init
  public :: opts_insert_file
  public :: opts_has
  public :: opts_get_raw
  public :: opts_get_logical
  public :: opts_get_integer
  public :: opts_get_real
  public :: opts_get_string
  public :: opts_get_integer_array
  public :: opts_get_real_array

  integer, parameter :: KEY_LEN  = 256
  integer, parameter :: VAL_LEN  = 4096
  integer, parameter :: MAX_OPTS = 2048

  integer           :: nopt = 0
  character(KEY_LEN) :: db_keys(MAX_OPTS)
  character(VAL_LEN) :: db_vals(MAX_OPTS)
  logical           :: db_initialized = .false.

contains

  !> Parse argv and auto-load tenstream.options. Idempotent.
  subroutine opts_init()
    integer :: i, argc
    character(VAL_LEN) :: arg, nxt
    logical :: file_exists

    if (db_initialized) return
    db_initialized = .true.

    argc = command_argument_count()
    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      if (is_key(arg)) then
        if (i < argc) then
          call get_command_argument(i + 1, nxt)
          if (.not. is_key(nxt)) then
            call store_opt(trim(arg(2:)), trim(nxt))
            i = i + 2
            cycle
          end if
        end if
        call store_opt(trim(arg(2:)), '')
      end if
      i = i + 1
    end do

    inquire(file='tenstream.options', exist=file_exists)
    if (file_exists) call opts_insert_file('tenstream.options')
  end subroutine

  !> Load options from file; silently skips missing files.
  !> Format: one option per line, -key [value]. # and ! start comments.
  subroutine opts_insert_file(filename)
    character(len=*), intent(in) :: filename
    character(VAL_LEN) :: line
    integer :: u, ios, cpos, sp

    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      cpos = scan(line, '#!')
      if (cpos > 0) line = line(:cpos - 1)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (.not. is_key(line)) cycle
      line = line(2:)  ! strip leading '-'
      sp = scan(line, ' '//char(9)//'=')
      if (sp == 0) then
        call store_opt(trim(line), '')
      else
        call store_opt(trim(line(:sp - 1)), trim(adjustl(line(sp + 1:))))
      end if
    end do
    close(u)
  end subroutine

  !> True if option exists. Tries prefix+name first, then bare name as fallback.
  subroutine opts_has(prefix, name, found)
    character(len=*), intent(in) :: prefix, name
    logical, intent(out) :: found
    found = find_opt(build_key(prefix, name)) > 0
    if (.not. found .and. len_trim(prefix) > 0) &
      found = find_opt(build_key('', name)) > 0
  end subroutine

  !> Raw string lookup. val is unchanged if option not found.
  subroutine opts_get_raw(prefix, name, val, found)
    character(len=*), intent(in)    :: prefix, name
    character(len=*), intent(inout) :: val
    logical, intent(out) :: found
    integer :: idx
    idx = find_opt(build_key(prefix, name))
    if (idx == 0 .and. len_trim(prefix) > 0) idx = find_opt(build_key('', name))
    if (idx > 0) then
      found = .true.
      val = trim(db_vals(idx))
    else
      found = .false.
    end if
  end subroutine

  !> Logical option. val unchanged if not found.
  !> A bare flag (-opt with no value) is treated as .true.
  subroutine opts_get_logical(prefix, name, val, found)
    character(len=*), intent(in) :: prefix, name
    logical, intent(inout) :: val
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    raw = ''
    call opts_get_raw(prefix, name, raw, found)
    if (found) val = parse_bool(raw)
  end subroutine

  !> Integer(int32) option. val unchanged if not found.
  subroutine opts_get_integer(prefix, name, val, found)
    character(len=*), intent(in) :: prefix, name
    integer(int32), intent(inout) :: val
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    raw = ''
    call opts_get_raw(prefix, name, raw, found)
    if (found) read(raw, *) val
  end subroutine

  !> Real(real64) option. val unchanged if not found.
  subroutine opts_get_real(prefix, name, val, found)
    character(len=*), intent(in) :: prefix, name
    real(real64), intent(inout) :: val
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    raw = ''
    call opts_get_raw(prefix, name, raw, found)
    if (found) read(raw, *) val
  end subroutine

  !> String option. val unchanged if not found.
  subroutine opts_get_string(prefix, name, val, found)
    character(len=*), intent(in) :: prefix, name
    character(len=*), intent(inout) :: val
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    raw = val
    call opts_get_raw(prefix, name, raw, found)
    if (found) val = trim(raw)
  end subroutine

  !> Integer(int32) array from comma-separated value.
  !> n: in=array capacity, out=count of values parsed.
  subroutine opts_get_integer_array(prefix, name, arr, n, found)
    character(len=*), intent(in) :: prefix, name
    integer(int32), intent(inout) :: arr(:)
    integer(int32), intent(inout) :: n
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    integer :: p, prev, ios
    integer(int32) :: tmp
    raw = ''
    call opts_get_raw(prefix, name, raw, found)
    if (.not. found) return
    raw = trim(raw)//','
    n = 0
    prev = 1
    do p = 1, len_trim(raw)
      if (raw(p:p) == ',') then
        read(raw(prev:p - 1), *, iostat=ios) tmp
        if (ios == 0 .and. n < int(size(arr), int32)) then
          n = n + 1
          arr(n) = tmp
        end if
        prev = p + 1
      end if
    end do
  end subroutine

  !> Real(real64) array from comma-separated value.
  !> n: in=array capacity, out=count of values parsed.
  subroutine opts_get_real_array(prefix, name, arr, n, found)
    character(len=*), intent(in) :: prefix, name
    real(real64), intent(inout) :: arr(:)
    integer(int32), intent(inout) :: n
    logical, intent(out) :: found
    character(VAL_LEN) :: raw
    integer :: p, prev, ios
    real(real64) :: tmp
    raw = ''
    call opts_get_raw(prefix, name, raw, found)
    if (.not. found) return
    raw = trim(raw)//','
    n = 0
    prev = 1
    do p = 1, len_trim(raw)
      if (raw(p:p) == ',') then
        read(raw(prev:p - 1), *, iostat=ios) tmp
        if (ios == 0 .and. n < int(size(arr), int32)) then
          n = n + 1
          arr(n) = tmp
        end if
        prev = p + 1
      end if
    end do
  end subroutine

  ! ---- private helpers ----

  pure function build_key(prefix, name) result(key)
    character(len=*), intent(in) :: prefix, name
    character(KEY_LEN) :: key
    integer :: n
    n = len_trim(name)
    if (n > 0 .and. name(1:1) == '-') then
      key = to_lower(trim(prefix)//name(2:n))
    else
      key = to_lower(trim(prefix)//name(:n))
    end if
  end function

  function find_opt(key) result(idx)
    character(len=*), intent(in) :: key
    integer :: idx, i
    character(KEY_LEN) :: k
    k = trim(key)
    idx = 0
    do i = 1, nopt
      if (db_keys(i) == k) then
        idx = i
        return
      end if
    end do
  end function

  subroutine store_opt(key, val)
    character(len=*), intent(in) :: key, val
    integer :: idx
    character(KEY_LEN) :: k
    k = to_lower(trim(key))
    if (len_trim(k) == 0) return
    idx = find_opt(k)
    if (idx == 0) then
      if (nopt >= MAX_OPTS) then
        write(error_unit, '(A)') 'opts_db: database full, ignoring -'//trim(k)
        return
      end if
      nopt = nopt + 1
      idx = nopt
    end if
    db_keys(idx) = k
    db_vals(idx) = trim(val)
  end subroutine

  !> True if arg looks like an option key: starts with '-' followed by a letter.
  pure function is_key(arg) result(res)
    character(len=*), intent(in) :: arg
    logical :: res
    integer :: c
    res = .false.
    if (len_trim(arg) < 2) return
    if (arg(1:1) /= '-') return
    c = iachar(arg(2:2))
    res = (c >= iachar('a') .and. c <= iachar('z')) .or. &
          (c >= iachar('A') .and. c <= iachar('Z'))
  end function

  pure function to_lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len(str)) :: out
    integer :: i, c
    out = str
    do i = 1, len(str)
      c = iachar(str(i:i))
      if (c >= 65 .and. c <= 90) out(i:i) = achar(c + 32)
    end do
  end function

  !> Parse boolean string. Empty string or 'true'/'yes'/'1' → .true.
  !> 'false'/'no'/'0' → .false.
  pure function parse_bool(str) result(val)
    character(len=*), intent(in) :: str
    logical :: val
    select case (trim(to_lower(adjustl(str))))
    case ('false', 'no', '0')
      val = .false.
    case default
      val = .true.
    end select
  end function

end module m_options_database
