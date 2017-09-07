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

  ncerr = nf90_get_var(ncid, varid, tmp) ; call nccheck(ncerr) ; ierr = ierr+ncerr
  ncerr = nf90_close(ncid)               ; call nccheck(ncerr) ; ierr = ierr+ncerr

  100 continue
  close(v,status='delete')
  if(ldebug) print *,'Reading Netcdf file returned with ierr: ',ierr

  if(allocated(tmp) ) then
    arr=tmp
    deallocate(tmp)
  endif
  return

  99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
  call cpusleep(waitinterval)
enddo
stop 'Couldnt read file .. waited now for quite a while but we couldnt open it :('

