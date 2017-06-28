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
integer, dimension(size(shape(arr))) :: dimIDs,dimsize

logical :: var_exists

integer :: ncerr
integer(iintegers),intent(out) :: ierr
ierr=0

if(ldebug) print *,'Looking for lock file for reading netcdf File',trim(groups(1))//'.lock'
lockfile=trim(groups(1))//'.lock'

fname=trim(groups(1))

do iwait=1,int(maxwait/waitinterval)
  open(v,file=lockfile,status='new',err=99)
  write(v,*) 'file is locked by process: ',get_pid_macro()

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
    ncerr = nf90_inquire_variable(ncid, varid, dimids = dimIDs) ; call nccheck(ncerr) ; ierr = ierr+ncerr

    do i=1,size(shape(arr))
      ncerr = nf90_inquire_dimension(ncid, dimIDs(i), len = dimsize(i)) ; call nccheck(ncerr) ;  ierr = ierr+ncerr
    enddo
  endif

