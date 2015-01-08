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

