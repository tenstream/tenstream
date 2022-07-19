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

logical, intent(in), optional :: lverbose
integer, dimension(:), optional, intent(in) :: ostart, ocount, ostride, omap

character(len=*) :: groups(:)

character(len=default_str_len) :: fname, varname
integer :: i, ncid, varid
integer, dimension(size(shape(arr))) :: dimIDs, dimsize

logical :: file_exists, var_exists

integer :: ncerr, flock_unit
integer(int64) :: bytesize
real(kind(arr)) :: bytesize_dummy
integer(mpiint), intent(out) :: ierr
ierr = 0

if (associated(arr)) then
  print *, 'Array which should have been read in is already allocated... I will not overwrite it...'
  ierr = 100
  return
end if

fname = trim(groups(1))

inquire (file=trim(fname), exist=file_exists)
if (.not. file_exists) then
  print *, ''
  print *, 'ERROR! NetCDF File ', trim(fname), ' could not be found!'
  print *, ''
  ierr = 101
  return
end if

call acquire_file_lock(fname, flock_unit, ierr); call CHKERR(ierr, 'Could not obtain file lock')

varname = trim(groups(2))
do i = 3, size(groups)
  varname = trim(varname)//'.'//trim(groups(i))
end do
if (ldebug .or. get_arg(.false., lverbose)) print *, 'reading from file: ', trim(fname)
if (ldebug .or. get_arg(.false., lverbose)) print *, 'variable: ', trim(varname), ' :: Ndims', size(shape(arr))

ncerr = nf90_open(trim(fname), NF90_NOWRITE, ncid); call nccheck(ncerr); ierr = ierr + ncerr
if (ncerr .ne. NF90_NOERR) goto 100
ncerr = nf90_inq_varid(ncid, trim(varname), varid); call nccheck(ncerr); ierr = ierr + ncerr
if (ncerr .ne. NF90_NOERR) goto 100

var_exists = ncerr .eq. NF90_NOERR
if (ldebug .or. get_arg(.false., lverbose)) print *, 'Reading NetCDF :: varname: ', trim(varname), ' :: exists?', var_exists

bytesize = c_sizeof(bytesize_dummy)
if (var_exists) then ! How big are the data dimensions?
  ncerr = nf90_inquire_variable(ncid, varid, dimids=dimIDs); call nccheck(ncerr); ierr = ierr + ncerr
  do i = 1, size(shape(arr))
    ncerr = nf90_inquire_dimension(ncid, dimIDs(i), len=dimsize(i)); call nccheck(ncerr); ierr = ierr + ncerr
    bytesize = bytesize * dimsize(i)
  end do
end if

if (ldebug .or. get_arg(.false., lverbose)) &
  print *, 'Dimension of variable: '//trim(varname), dimsize, '( '//toStr(real(bytesize) / 1024.**3)//' Gb )'

