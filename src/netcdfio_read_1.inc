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


character(len=*) :: groups(:)

character(len=default_str_len) :: fname,varname
integer :: i,ncid,varid
integer, dimension(size(shape(arr))) :: dimIDs, dimsize

logical :: file_exists, var_exists

integer :: ncerr, flock_unit, ndims
integer(mpiint),intent(out) :: ierr
ierr=0

if(allocated(arr)) then
  print *,'Array which should have been read in is already allocated... I will not overwrite it...'
  ierr = 100
  return
endif

fname=trim(groups(1))

inquire( file=trim(fname), exist=file_exists )
if(.not.file_exists) then
  print *,''
  print *,'ERROR! NetCDF File ',trim(fname),' could not be found!'
  print *,''
  ierr = 101
  return
endif

call acquire_file_lock(fname, flock_unit, ierr); call CHKERR(ierr, 'Could not obtain file lock')

varname=trim(groups(2))
do i=3,size(groups)
  varname=trim(varname)//'.'//trim(groups(i))
enddo
if(ldebug) print *,'reading from file: ',trim(fname)
if(ldebug) print *,'variable: ',trim(varname),' :: Ndims',size(shape(arr))

ncerr = nf90_open(trim(fname), NF90_NOWRITE, ncid) ; call nccheck(ncerr) ; ierr = ierr+ncerr
if(ncerr.ne.NF90_NOERR) goto 100
ncerr = nf90_inq_varid(ncid, trim(varname), varid) ; call nccheck(ncerr) ; ierr = ierr+ncerr
if(ncerr.ne.NF90_NOERR) goto 100

var_exists = ncerr.eq.NF90_NOERR
if(ldebug) print *,'Reading NetCDF :: varname: ',trim(varname), ' :: exists?',var_exists

if(var_exists) then ! How big are the data dimensions?
  ncerr = nf90_inquire_variable(ncid, varid, ndims = ndims) ; call nccheck(ncerr) ; ierr = ierr+ncerr
  if(ndims.ne.size(shape(arr))) then
    if(ldebug) print *,'Input and NC variable have different number dimensions '//itoa(size(shape(arr)))//' vs '//itoa(ndims)
    ierr = ierr + 100 * size(shape(arr)) + ndims
    goto 100
  endif

  ncerr = nf90_inquire_variable(ncid, varid, dimids = dimIDs) ; call nccheck(ncerr) ; ierr = ierr+ncerr
  do i=1,size(shape(arr))
    ncerr = nf90_inquire_dimension(ncid, dimIDs(i), len = dimsize(i)) ; call nccheck(ncerr) ;  ierr = ierr+ncerr
  enddo
endif

if(ldebug) print *,'Dimension of variable:',dimsize

