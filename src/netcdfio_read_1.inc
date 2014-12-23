character(len=*) :: groups(:)

character(len=300) :: fname,varname
integer :: i,ncid,varid
integer, dimension(size(shape(arr))) :: dimIDs,dimsize

logical :: var_exists

integer :: ncerr
integer(iintegers),intent(out) :: ierr
ierr=0

lockfile=trim(groups(1))//'.lock'

fname=trim(groups(1))

do iwait=1,int(maxwait/waitinterval)
  open(v,file=lockfile,status='new',err=99)
  write(v,*) 'file is locked by process: ',getpid()

  varname=trim(groups(2))
  do i=3,size(groups)
    varname=trim(varname)//'.'//trim(groups(i))
  enddo
  if(ldebug) print *,'reading from file: ',trim(fname),'variable: ',trim(varname),' ::Ndims',size(shape(arr))

  ncerr = nf90_open(trim(fname), NF90_NOWRITE, ncid) ; call nccheck(ncerr) ; ierr = ierr+ncerr 
  ncerr = nf90_inq_varid(ncid, trim(varname), varid) ; call nccheck(ncerr) ; ierr = ierr+ncerr

  var_exists = ncerr.eq.NF90_NOERR
  if(ldebug) print *,'Reading NetCDF :: varname: ',trim(varname), ' :: exists?',var_exists

  if(var_exists) then ! How big are the data dimensions?
    ncerr = nf90_inquire_variable(ncid, varid, dimids = dimIDs) ; call nccheck(ncerr) ; ierr = ierr+ncerr

    do i=1,size(shape(arr))
      ncerr = nf90_inquire_dimension(ncid, dimIDs(i), len = dimsize(i)) ; call nccheck(ncerr) ;  ierr = ierr+ncerr
    enddo
  endif

