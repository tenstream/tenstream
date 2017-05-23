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
integer :: Ndims,ncid,i,varid

integer,allocatable :: dims(:),dim_ids(:)
character(len=10) :: dimname

logical :: file_exists, var_exists, def_mode

integer(iintegers),intent(out) :: ierr
lockfile=trim(groups(1))//'.lock'

fname=trim(groups(1))

do iwait=1,int(maxwait/waitinterval)
open(v,file=lockfile,status='new',err=99)
write(v,*) 'file is locked by process: ',get_pid_macro()

varname=trim(groups(2))
do i=3,size(groups)
varname=trim(varname)//'.'//trim(groups(i))
enddo

if(ldebug) print *,'***********************************************************************************************'
if(ldebug) print *,'writing to file: ',trim(fname)
if(ldebug) print *,'varname is: ',trim(varname)
if(ldebug) print *,'***********************************************************************************************'


Ndims = size(shape(arr))

allocate(dims   (Ndims))
dims = shape(arr)
allocate(dim_ids(Ndims))

inquire(file=fname, exist=file_exists)
if(ldebug) print *,'file ',trim(fname),'exists ::', file_exists

if(.not.file_exists) then
  ierr = nf90_create(fname, NF90_NOCLOBBER, ncid) ; call nccheck(ierr)
  !ierr = nf90_create(fname, NF90_HDF5, ncid) ; call nccheck(ierr)
  if(ldebug) print *,'opened file to append: ',trim(fname),' ::',ierr
else
  ierr = nf90_open(fname, NF90_WRITE, ncid) ; call nccheck(ierr)
  if(ldebug) print *,'opened file to write: ',trim(fname),' ::',ierr
endif
if(ierr.ne.NF90_NOERR) stop 'stopping because we couldnt write to file'
ierr = nf90_redef(ncid); call nccheck(ierr); def_mode=.True.
if(ldebug) print *,'Set NetCDF File to def_mode',ierr

do i=1,Ndims
  write(dimname,FMT='("dim",I0)') i

  !See if we can find the dimension
  ierr = nf90_inq_dimid(ncid, trim(varname)//'_'//trim(dimname), dim_ids(i)) ; call nccheck(ierr)
  if(ldebug) print *,'Check if dimension exists',ierr
  if(ierr.ne.NF90_NOERR) then !Otherwise create it
    ierr = nf90_def_dim(ncid, trim(varname)//'_'//trim(dimname), dims(i), dim_ids(i) ); call nccheck(ierr)
    if(ldebug) print *,'Create dimension :: ',trim(varname)//'_'//trim(dimname),' :: ',ierr
  endif
enddo

!Check if variable already exists
ierr = nf90_inq_varid(ncid, varname, varid) ; call nccheck(ierr)
var_exists = ierr.eq.NF90_NOERR
if(ldebug) print *,'Check if variable exists',ierr,'::',var_exists

if(.not.var_exists) then !Otherwise create it
  if(ldebug) print *,'Create variable'
  ierr = nf90_def_var(ncid, varname, NF90_FLOAT, dim_ids, varid,deflate_level=deflate_lvl) ; call nccheck(ierr)
endif

if(def_mode) then !if we created the file here, we need to leave define mode...
  ierr = nf90_enddef(ncid); call nccheck(ierr )
  if(ldebug) print *,'Putting netcdf file in datamode',ierr
endif

tmp=arr
ierr = nf90_put_var(ncid, varid, tmp) ; call nccheck(ierr)
if(ldebug) print *,'written data',ierr

ierr = nf90_close(ncid) ; call nccheck(ierr)
if(ldebug) print *,'closed file',ierr

close(v,status='delete')
return
99 continue
if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
call cpusleep(waitinterval)
enddo
stop 'Couldnt write file .. waited now for quite a while but we couldnt write to it :('


