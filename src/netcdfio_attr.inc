    !subroutine set_global_attribute_r32(fname, attr_name, attr)
      character(len=*), intent(in) :: fname, attr_name
      integer :: ncid, ierr
      ierr = nf90_open(trim(fname), nf90_write, ncid); call nccheck(ierr)
      ierr = nf90_redef(ncid); call nccheck(ierr)
      ierr = nf90_put_att (ncid, nf90_global, trim(attr_name), attr); call nccheck(ierr)
      ierr = nf90_enddef(ncid); call nccheck(ierr)
      ierr = nf90_close(ncid); call nccheck(ierr)
    !end subroutine

