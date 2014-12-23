  ierr = nf90_get_var(ncid, varid, arr) ; call nccheck(ierr)
  ierr = nf90_close(ncid) ; call nccheck(ierr)

  close(v,status='delete')
  return
  99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
  call cpusleep(waitinterval)
enddo


