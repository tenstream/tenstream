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

integer, intent(in), optional :: arr_shape(:), startp(:), countp(:), stride(:), map(:)

character(len=*) :: groups(:)

character(len=default_str_len) :: fname,varname
integer :: Ndims,ncid,i,varid

integer,allocatable :: dims(:),dim_ids(:)
character(len=default_str_len) :: dimname

logical :: file_exists, var_exists, def_mode

integer(mpiint),intent(out) :: ierr
integer(mpiint) :: ncerr, flock_unit
fname=trim(groups(1))

call acquire_file_lock(fname, flock_unit, ierr); call CHKERR(ierr, 'Could not obtain file lock')

varname=trim(groups(2))
do i=3,size(groups)
varname=trim(varname)//'.'//trim(groups(i))
enddo

if(ldebug) print *,'***********************************************************************************************'
if(ldebug) print *,'writing to file: ',trim(fname)
if(ldebug) print *,'varname is: ',trim(varname)
if(ldebug) print *,'***********************************************************************************************'

if(present(arr_shape)) then
  Ndims = size(arr_shape)
else
  Ndims = size(shape(arr))
endif
allocate(dims(Ndims))

if(present(arr_shape)) then
  dims = arr_shape
else
  dims = shape(arr)
endif

allocate(dim_ids(Ndims))

if(ldebug) print *,'Dimensions('//itoa(int(Ndims, iintegers))//') of array', dims
if(ldebug) then
  ncerr = nf90_set_log_level(2); call nccheck(ncerr)
endif

inquire(file=fname, exist=file_exists)
if(ldebug) print *,'file '//trim(fname)//' exists :: ', file_exists

ierr = 0

if(.not.file_exists) then
  !ncerr = nf90_create(fname, NF90_NOCLOBBER, ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
  !ncerr = nf90_create(fname, NF90_NETCDF4, ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
  !ncerr = nf90_create(fname, NF90_64BIT_OFFSET , ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
  ncerr = nf90_create(trim(fname), NF90_HDF5, ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
  if(ldebug) print *,'opened file to append: ',trim(fname),' ::',ncerr
else
  ncerr = nf90_open(path=trim(fname), mode=NF90_WRITE, ncid=ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
  if(ldebug) print *,'opened file to write: ',trim(fname),' ::',ncerr
endif
if(ncerr.ne.NF90_NOERR) call CHKERR(int(ncerr, mpiint),'stopping because we couldnt write to file')
ierr = 0

ncerr = nf90_redef(ncid); call nccheck(ncerr); def_mode=.True.
if(ldebug) print *,'Set NetCDF File to def_mode',ncerr

do i=1,Ndims
  write(dimname,FMT='("dim",I0)') i

  !See if we can find the dimension
  ncerr = nf90_inq_dimid(ncid, trim(varname)//'_'//trim(dimname), dim_ids(i)) ; call nccheck(ncerr)
  if(ldebug) print *,'Check if dimension '//trim(varname)//'_'//trim(dimname)//' exists',ncerr
  if(ncerr.ne.NF90_NOERR) then !Otherwise create it
    ncerr = nf90_def_dim(ncid, trim(varname)//'_'//trim(dimname), dims(i), dim_ids(i) ); call nccheck(ncerr)
    if(ldebug) print *,'Create dimension :: ',i,dims(i),dim_ids(i),trim(varname)//'_'//trim(dimname),' :: ',ncerr
  endif
enddo

!Check if variable already exists
ncerr = nf90_inq_varid(ncid, varname, varid) ; call nccheck(ncerr)
var_exists = ncerr.eq.NF90_NOERR
if(ldebug) print *,'Check if variable exists',ncerr,'::',var_exists

if(.not.var_exists) then !Otherwise create it
  if(ldebug) print *,'Create variable '//itoa(varid)
  ncerr = nf90_def_var(ncid, varname, NF90_FLOAT, dim_ids, varid, deflate_level=deflate_lvl) ; call nccheck(ncerr)
  ierr = ierr + ncerr
  call CHKWARN(ncerr, 'Could not create new variable '//trim(varname)//' varid: '//itoa(varid))
endif


if(def_mode) then !if we created the file here, we need to leave define mode...
  ncerr = nf90_enddef(ncid); call nccheck(ncerr )
  if(ldebug) print *,'Putting netcdf file in datamode',ncerr
  call CHKWARN(ncerr, 'Could not put netcdf file in datamode')
endif

if(ldebug) print *,'writing data'
ncerr = nf90_put_var(ncid, varid, values=arr, &
  start=startp, count=countp, stride=stride, map=map)
if(ldebug) print *,'written data',ncerr

ncerr = nf90_close(ncid) ; call nccheck(ncerr); ierr = ierr + ncerr
if(ldebug) print *,'closed file',ncerr

call release_file_lock(flock_unit, ncerr); call CHKERR(ncerr, 'Could not release filelock')
