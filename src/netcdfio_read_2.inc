  ncerr = nf90_get_var(ncid, varid, arr) ; call nccheck(ncerr) ; ierr = ierr+ncerr
  ncerr = nf90_close(ncid)               ; call nccheck(ncerr) ; ierr = ierr+ncerr

  close(v,status='delete')
  if(ldebug) print *,'Reading Netcdf file returned with ierr: ',ierr
  return
  99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
  call cpusleep(waitinterval)
enddo



