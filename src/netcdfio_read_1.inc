character(len=*) :: groups(:)

character(len=300) :: fname,varname
integer :: i,ncid,varid
integer, dimension(size(shape(arr))) :: dimIDs,dimsize

integer(iintegers),intent(out) :: ierr
lockfile=trim(groups(1))//'.lock'

fname=trim(groups(1))

do iwait=1,int(maxwait/waitinterval)
  open(v,file=lockfile,status='new',err=99)
  write(v,*) 'file is locked by process: ',getpid()

  varname=trim(groups(2))
  do i=3,size(groups)
    varname=trim(varname)//'.'//trim(groups(i))
  enddo
  print *,'reading from file: ',fname,'variable: ',varname,' ::Ndims',size(shape(arr))

  ierr = nf90_open(fname, NF90_NOWRITE, ncid) ; call nccheck(ierr)
  ierr = nf90_inq_varid(ncid, varname, varid) ; call nccheck(ierr) 

  ! How big are the data dimensions?
  ierr = nf90_inquire_variable(ncid, varid, dimids = dimIDs) ; call nccheck(ierr)

  do i=1,size(shape(arr))
    ierr = nf90_inquire_dimension(ncid, dimIDs(i), len = dimsize(i)) ; call nccheck(ierr)
  enddo


