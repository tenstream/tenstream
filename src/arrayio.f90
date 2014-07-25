module arrayIO
  use mpi
  use HDF5
  USE data_parameters, ONLY :   &
      ireals,    &
      iintegers
  implicit none

  private
  public :: h5write,h5load,write_bin,write_ascii,read_bin

  integer :: u=10,v=11
  real(ireals),parameter :: maxwait=600 !in seconds
  real(ireals),parameter :: waitinterval=.01 ! amount of cpu time to wait before trying anew in seconds 
  integer :: iwait
  character(310) :: lockfile

  interface h5write
    module procedure h5write_1d,h5write_2d,h5write_3d,h5write_4d,h5write_5d,h5write_7d
  end interface
  interface h5load
    module procedure h5load_1d,h5load_2d,h5load_3d,h5load_4d,h5load_5d,h5load_7d, h5load_1d_int,h5load_2d_int
  end interface
  interface write_bin
    module procedure wb1d,wb2d,wb3d,wb4d
  end interface
  interface write_ascii
    module procedure wa1d,wa2d,wa3d
  end interface
  interface read_bin
    module procedure rb1d,rb2d,rb3d,rb4d
  end interface

  contains 

    subroutine wa1d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:)
        character(10) :: suffbin='.txt',suffdescr='.txt.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        open (unit=v, file=trim(fname), status="replace")
        write(v, *) a
        close(v)
    end subroutine
    subroutine wa2d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:,:)
        character(10) :: suffbin='.txt',suffdescr='.txt.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        open (unit=v, file=trim(fname), status="replace")
        write(v, *) a
        close(v)
    end subroutine
    subroutine wa3d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:,:,:)
        character(10) :: suffbin='.txt',suffdescr='.txt.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        open (newunit=v, file=trim(fname), status="replace")
        write(v, *) a
        close(v)
    end subroutine

    subroutine wb1d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:)
        integer :: reclen
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        inquire(iolength=reclen) a
        open (unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        write(unit=v,rec=1) a
        close(unit=v)
    end subroutine
    subroutine rb1d(basename,o,ierr)
        character(100),intent(in) :: basename
        real(ireals),allocatable :: o(:)
        integer :: reclen,dims(1)
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        logical :: file_exists
        integer(iintegers),intent(out) :: ierr
        ierr=0
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)
        inquire(file=fname, exist=file_exists)
        if(.not.file_exists) then
          ierr=5
          print *,'file doesnt exist: ',trim(fname),' : ',trim(descr)
          return
        endif

        !        print *,'reading description from',trim(descr)
        open (newunit=u, file=trim(descr), status="old")
        read (u, *) dims
        close(u)
        !        print *,'dimensions of array are',dims
        allocate(o(dims(1)))

        inquire(iolength=reclen) o
        open(unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        read(unit=v,rec=1) o
        close(unit=v)
    end subroutine
    subroutine wb2d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:,:)
        integer :: reclen
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        inquire(iolength=reclen) a
        open (unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        write(unit=v,rec=1) a
        close(unit=v)
    end subroutine
    subroutine rb2d(basename,o,ierr)
        character(100),intent(in) :: basename
        real(ireals),allocatable :: o(:,:)
        integer :: reclen,dims(2)
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        logical :: file_exists
        integer(iintegers),intent(out) :: ierr
        ierr=0
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)
        inquire(file=fname, exist=file_exists)
        if(.not.file_exists) then
          ierr=5
          print *,'file doesnt exist: ',trim(fname),' : ',trim(descr)
          return
        endif

        !        print *,'reading description from',trim(descr)
        open (newunit=u, file=trim(descr), status="old")
        read (u, *) dims
        close(u)
        !        print *,'dimensions of array are',dims
        allocate(o(dims(1),dims(2)))

        inquire(iolength=reclen) o
        open(unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        read(unit=v,rec=1) o
        close(unit=v)
    end subroutine
    subroutine wb3d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:,:,:)
        integer :: reclen
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        inquire(iolength=reclen) a
        open (unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        write(unit=v,rec=1) a
        close(unit=v)
    end subroutine
    subroutine rb3d(basename,o,ierr)
        character(100),intent(in) :: basename
        real(ireals),allocatable :: o(:,:,:)
        integer :: reclen,dims(3)
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        logical :: file_exists
        integer(iintegers),intent(out) :: ierr
        ierr=0
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)
        inquire(file=fname, exist=file_exists)
        if(.not.file_exists) then
          ierr=5
          print *,'file doesnt exist: ',trim(fname),' : ',trim(descr)
          return
        endif

        !        print *,'reading description from',trim(descr)
        open (newunit=u, file=trim(descr), status="old")
        read (u, *) dims
        close(u)
        !        print *,'dimensions of array are',dims
        allocate(o(dims(1),dims(2),dims(3)))

        inquire(iolength=reclen) o
        open(unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        read(unit=v,rec=1) o
        close(unit=v)
    end subroutine
    subroutine wb4d(basename,a)
        character(100),intent(in) :: basename
        real(ireals),intent(in) :: a(:,:,:,:)
        integer :: reclen
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)

        print *,'writing data to:  ',trim(fname),'   ; shape to:   ',trim(descr)

        open (newunit=u, file=trim(descr), status="replace")
        write(u, *) shape(a)
        close(u)

        inquire(iolength=reclen) a
        open (unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        write(unit=v,rec=1) a
        close(unit=v)
    end subroutine
    subroutine rb4d(basename,o,ierr)
        character(100),intent(in) :: basename
        real(ireals),allocatable :: o(:,:,:,:)
        integer :: reclen,dims(4)
        character(10) :: suffbin='.bin',suffdescr='.bin.descr'
        character(100) :: fname,descr
        logical :: file_exists
        integer(iintegers),intent(out) :: ierr
        ierr=0
        fname = trim(basename) // trim(suffbin)
        descr = trim(basename) // trim(suffdescr)
        inquire(file=fname, exist=file_exists)
        if(.not.file_exists) then
          ierr=5
          print *,'file doesnt exist: ',trim(fname),' : ',trim(descr)
          return
        endif

        !        print *,'reading description from',trim(descr)
        open (newunit=u, file=trim(descr), status="old")
        read (u, *) dims
        close(u)
        !        print *,'dimensions of array are',dims
        allocate(o(dims(1),dims(2),dims(3),dims(4)))

        inquire(iolength=reclen) o
        open(unit=v, file=trim(fname), form='unformatted',access='direct',recl=reclen)
        read(unit=v,rec=1) o
        close(unit=v)
    end subroutine

    subroutine h5write_open(groups,dims,id,dataset,dataspace,ierr)
        character(len=*) :: groups(:)
        integer(HSIZE_T),intent(in) :: dims(:)
        integer(HID_T) :: id(:),dataset,dataspace
        integer(iintegers),intent(out) :: ierr

        integer(HSIZE_T) :: chunk(size(dims) )

        integer :: k,lastid,hferr
        logical :: compression,file_exists,link_exists
        integer(HID_T) :: dcpl

        ierr=0 ; lastid = ubound(id,1)
        call h5open_f(hferr); ierr=ierr+hferr
        if(ierr.ne.0) return
        chunk=max(1, dims/2)

!        print *,'Hdf5 chunks are:',chunk
        if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
            &  and scnd is at least hdf5 root /, and third is data name'
        !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))

        compression=check_for_zip_filter()

        inquire(file=trim(groups(1)), exist=file_exists)
        if(.not.file_exists) then
          !                print *,'File does not yet exists, creating new one'
          call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
        else   
          !                print *,'File does exist, opening it RW'
          call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
        endif

        !        print *,'Creating Compression Filter'
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
        call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
        call h5pset_chunk_f(dcpl, size(dims), chunk, hferr) ; ierr=ierr+hferr

        do k=2,lastid
          call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
          if(link_exists) then
            !                        print *,'opening group',trim(groups(k))
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          else
            !                        print *,'Create Groups',trim(groups(k))
            call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          endif
        enddo

        !        print *,'Creating DataSpace'
        call h5screate_simple_f(size(dims), dims, dataspace, hferr) ; ierr=ierr+hferr
        if(hferr.ne.0) return

        !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
        call h5lexists_f(id(lastid), trim(groups(lastid+1)), link_exists, hferr) ; ierr=ierr+hferr
        if(link_exists) then
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
        else
          call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
        endif

        if(ierr.ne.0) return
    end subroutine
    subroutine h5write_close(id, dataset,dataspace,ierr)
        integer(HID_T) :: id(:),dataset,dataspace
        integer(iintegers),intent(out) :: ierr
        integer :: k,hferr

          !        print *,'Closing handles'
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=ubound(id,1),2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr

          !        print *,'Data is now:',sum(arr)/size(arr)

          call h5close_f(hferr)  ; ierr=ierr+hferr
    end subroutine
    subroutine h5write_1d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 1 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_2d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 2 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_3d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 3 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_4d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 4 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_5d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:,:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 5 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_6d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:,:,:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 6 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5write_7d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),intent(in) :: arr(:,:,:,:,:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer :: hferr
        integer(HSIZE_T) :: dims( 7 )

        integer(HID_T) :: id(size(groups)-1),dataset,dataspace

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()

          dims = shape(arr)
          call h5write_open(groups,dims,id,dataset,dataspace,ierr)

          !        print *,'Write to dataset'
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr

          call h5write_close(id, dataset,dataspace,ierr)

          close(v,status='delete')
          return
          99 if(rand() .gt. .999) write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
!    subroutine h5write_2d(groups,arr,ierr)
!        character(len=*) :: groups(:)
!        real(ireals),intent(in) :: arr(:,:)
!        integer(iintegers),intent(out) :: ierr
!
!        integer,parameter :: rank=2
!        integer :: k,lastid,hferr
!        integer(HSIZE_T) :: dims(rank),chunk(rank)
!        logical :: compression,file_exists,link_exists
!
!        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
!        integer(HID_T) :: dcpl
!
!        lockfile=trim(groups(1))//'.lock'
!        do iwait=1,int(maxwait/waitinterval)
!          open(v,file=lockfile,status='new',err=99)
!          write(v,*) 'file is locked by process: ',getpid()
!
!          call h5open_f(hferr); ierr=ierr+hferr
!
!          ierr=0 ; lastid = ubound(id,1)
!          dims=[ubound(arr,1),ubound(arr,2)]
!          chunk=[dims(1),dims(2)]/2
!          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
!              &  and scnd is at least hdf5 root /, and third is data name'
!          !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))
!
!          compression=check_for_zip_filter()
!
!          inquire(file=trim(groups(1)), exist=file_exists)
!          if(.not.file_exists) then
!            !                print *,'File does not yet exists, creating new one'
!            call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
!          else   
!            !                print *,'File does exist, opening it RW'
!            call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
!          endif
!
!          !        print *,'Creating Compression Filter'
!          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
!          call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
!          call h5pset_chunk_f(dcpl, rank, chunk, hferr) ; ierr=ierr+hferr
!
!          do k=2,lastid
!            call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
!            if(link_exists) then
!              !                        print *,'opening group',trim(groups(k))
!              call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            else
!              !                        print *,'Create Groups',trim(groups(k))
!              call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            endif
!          enddo
!
!          !        print *,'Creating DataSpace'
!          call h5screate_simple_f(rank, dims, dataspace, hferr) ; ierr=ierr+hferr
!          if(hferr.ne.0) return
!
!          !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
!          call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
!
!          if(ierr.ne.0) return
!          !        print *,'Write to dataset'
!          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr
!
!          !        print *,'Closing handles'
!          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
!          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
!          do k=lastid,2,-1
!            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
!          enddo
!          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
!
!          !        print *,'Data is now:',sum(arr)/size(arr)
!
!          call h5close_f(hferr)  ; ierr=ierr+hferr
!          close(v,status='delete')
!          return
!          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
!          call cpusleep(waitinterval)
!        enddo
!    end subroutine
!    subroutine h5write_3d(groups,arr,ierr)
!        character(len=*) :: groups(:)
!        real(ireals),intent(in) :: arr(:,:,:)
!        integer(iintegers),intent(out) :: ierr
!
!        integer,parameter :: rank=3
!        integer :: k,lastid,hferr
!        integer(HSIZE_T) :: dims(rank),chunk(rank)
!        logical :: compression,file_exists,link_exists
!
!        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
!        integer(HID_T) :: dcpl
!
!        lockfile=trim(groups(1))//'.lock'
!        do iwait=1,int(maxwait/waitinterval)
!          open(v,file=lockfile,status='new',err=99)
!          write(v,*) 'file is locked by process: ',getpid()
!
!          call h5open_f(hferr); ierr=ierr+hferr
!
!          ierr=0 ; lastid = ubound(id,1)
!          dims=[ubound(arr,1),ubound(arr,2),ubound(arr,3)]
!          chunk=[dims(1),dims(2),dims(3)]/2
!          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
!              &  and scnd is at least hdf5 root /, and third is data name'
!          !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))
!
!          compression=check_for_zip_filter()
!
!          inquire(file=trim(groups(1)), exist=file_exists)
!          if(.not.file_exists) then
!            !                print *,'File does not yet exists, creating new one'
!            call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
!          else   
!            !                print *,'File does exist, opening it RW'
!            call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
!          endif
!
!          !        print *,'Creating Compression Filter'
!          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
!          call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
!          call h5pset_chunk_f(dcpl, rank, chunk, hferr) ; ierr=ierr+hferr
!
!          do k=2,lastid
!            call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
!            if(link_exists) then
!              !                        print *,'opening group',trim(groups(k))
!              call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            else
!              !                        print *,'Create Groups',trim(groups(k))
!              call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            endif
!          enddo
!
!          !        print *,'Creating DataSpace'
!          call h5screate_simple_f(rank, dims, dataspace, hferr) ; ierr=ierr+hferr
!          if(hferr.ne.0) return
!
!          !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
!          call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
!
!          if(ierr.ne.0) return
!          !        print *,'Write to dataset'
!          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr
!
!          !        print *,'Closing handles'
!          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
!          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
!          do k=lastid,2,-1
!            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
!          enddo
!          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
!
!          !        print *,'Data is now:',sum(arr)/size(arr)
!
!          call h5close_f(hferr)  ; ierr=ierr+hferr
!          close(v,status='delete')
!          return
!          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
!          call cpusleep(waitinterval)
!        enddo
!    end subroutine
!    subroutine h5write_4d(groups,arr,ierr)
!        character(len=*) :: groups(:)
!        real(ireals),intent(in) :: arr(:,:,:,:)
!        integer(iintegers),intent(out) :: ierr
!
!        integer,parameter :: rank=4
!        integer :: k,lastid,hferr
!        integer(HSIZE_T) :: dims(rank),chunk(rank)
!        logical :: compression,file_exists,link_exists
!
!        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
!        integer(HID_T) :: dcpl
!
!        lockfile=trim(groups(1))//'.lock'
!        do iwait=1,int(maxwait/waitinterval)
!          open(v,file=lockfile,status='new',err=99)
!          write(v,*) 'file is locked by process: ',getpid()
!
!          call h5open_f(hferr); ierr=ierr+hferr
!
!          ierr=0 ; lastid = ubound(id,1)
!          dims=[ubound(arr,1),ubound(arr,2),ubound(arr,3),ubound(arr,4)]
!          chunk=[dims(1),dims(2),dims(3),dims(4)]/2
!          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
!              &  and scnd is at least hdf5 root /, and third is data name'
!          !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))
!
!          compression=check_for_zip_filter()
!
!          inquire(file=trim(groups(1)), exist=file_exists)
!          if(.not.file_exists) then
!            !                print *,'File does not yet exists, creating new one'
!            call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
!          else   
!            !                print *,'File does exist, opening it RW'
!            call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
!          endif
!
!          !        print *,'Creating Compression Filter'
!          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
!          call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
!          call h5pset_chunk_f(dcpl, rank, chunk, hferr) ; ierr=ierr+hferr
!
!          do k=2,lastid
!            call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
!            if(link_exists) then
!              !                        print *,'opening group',trim(groups(k))
!              call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            else
!              !                        print *,'Create Groups',trim(groups(k))
!              call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            endif
!          enddo
!
!          !        print *,'Creating DataSpace'
!          call h5screate_simple_f(rank, dims, dataspace, hferr) ; ierr=ierr+hferr
!          if(hferr.ne.0) return
!
!          !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
!          call h5lexists_f(id(lastid), trim(groups(lastid+1)), link_exists, hferr) ; ierr=ierr+hferr
!          if(link_exists) then
!            call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
!          else
!            call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
!          endif
!
!          if(ierr.ne.0) return
!          !        print *,'Write to dataset'
!          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr
!
!          !        print *,'Closing handles'
!          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
!          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
!          do k=lastid,2,-1
!            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
!          enddo
!          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
!
!          !        print *,'Data is now:',sum(arr)/size(arr)
!
!          call h5close_f(hferr)  ; ierr=ierr+hferr
!          close(v,status='delete')
!          return
!          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
!          call cpusleep(waitinterval)
!        enddo
!    end subroutine
!    subroutine h5write_5d(groups,arr,ierr)
!        character(len=*) :: groups(:)
!        real(ireals),intent(in) :: arr(:,:,:,:,:)
!        integer(iintegers),intent(out) :: ierr
!
!        integer,parameter :: rank=5
!        integer :: k,lastid,hferr
!        integer(HSIZE_T) :: dims(rank),chunk(rank)
!        logical :: compression,file_exists,link_exists
!
!        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
!        integer(HID_T) :: dcpl
!
!        lockfile=trim(groups(1))//'.lock'
!        do iwait=1,int(maxwait/waitinterval)
!          open(v,file=lockfile,status='new',err=99)
!          write(v,*) 'file is locked by process: ',getpid()
!
!          call h5open_f(hferr); ierr=ierr+hferr
!
!          ierr=0 ; lastid = ubound(id,1)
!          dims=[ubound(arr,1),ubound(arr,2),ubound(arr,3),ubound(arr,4),ubound(arr,5)]
!          chunk=max(1_HSIZE_T,[dims(1),dims(2),dims(3),dims(4),dims(5)]/2)
!          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
!              &  and scnd is at least hdf5 root /, and third is data name'
!          !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))
!
!          compression=check_for_zip_filter()
!
!          inquire(file=trim(groups(1)), exist=file_exists)
!          if(.not.file_exists) then
!            !                print *,'File does not yet exists, creating new one'
!            call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
!          else   
!            !                print *,'File does exist, opening it RW'
!            call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
!          endif
!
!          !        print *,'Creating Compression Filter'
!          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
!          call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
!          call h5pset_chunk_f(dcpl, rank, chunk, hferr) ; ierr=ierr+hferr
!
!          do k=2,lastid
!            call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
!            if(link_exists) then
!              !                        print *,'opening group',trim(groups(k))
!              call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            else
!              !                        print *,'Create Groups',trim(groups(k))
!              call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            endif
!          enddo
!
!          !        print *,'Creating DataSpace'
!          call h5screate_simple_f(rank, dims, dataspace, hferr) ; ierr=ierr+hferr
!          if(hferr.ne.0) return
!
!          !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
!          call h5lexists_f(id(lastid), trim(groups(lastid+1)), link_exists, hferr) ; ierr=ierr+hferr
!          if(link_exists) then
!            call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
!          else
!            call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
!          endif
!
!          if(ierr.ne.0) return
!          !        print *,'Write to dataset'
!          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr
!
!          !        print *,'Closing handles'
!          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
!          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
!          do k=lastid,2,-1
!            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
!          enddo
!          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
!
!          !        print *,'Data is now:',sum(arr)/size(arr)
!
!          call h5close_f(hferr)  ; ierr=ierr+hferr
!          close(v,status='delete')
!          return
!          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
!          call cpusleep(waitinterval)
!        enddo
!    end subroutine
!    subroutine h5write_7d(groups,arr,ierr)
!        character(len=*) :: groups(:)
!        real(ireals),intent(in) :: arr(:,:,:,:,:,:,:)
!        integer(iintegers),intent(out) :: ierr
!
!        integer,parameter :: rank=7
!        integer :: k,lastid,hferr
!        integer(HSIZE_T) :: dims(rank),chunk(rank)
!        logical :: compression,file_exists,link_exists
!
!        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
!        integer(HID_T) :: dcpl
!
!        lockfile=trim(groups(1))//'.lock'
!        do iwait=1,int(maxwait/waitinterval)
!          open(v,file=lockfile,status='new',err=99)
!          write(v,*) 'file is locked by process: ',getpid()
!
!          call h5open_f(hferr); ierr=ierr+hferr
!
!          ierr=0 ; lastid = ubound(id,1)
!          dims = shape(arr) ![ubound(arr,1),ubound(arr,2),ubound(arr,3),ubound(arr,4),ubound(arr,5),ubound(arr,6),ubound(arr,7)]
!          chunk = max(1, [dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)]/2 )
!          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
!              &  and scnd is at least hdf5 root /, and third is data name'
!          !        print *,'writing hdf5 file ',trim(groups(1)), 'with groups: ',groups(2:ubound(groups,1))
!
!          compression=check_for_zip_filter()
!
!          inquire(file=trim(groups(1)), exist=file_exists)
!          if(.not.file_exists) then
!            !                print *,'File does not yet exists, creating new one'
!            call h5fcreate_f(trim(groups(1)), H5F_ACC_TRUNC_F, id(1), hferr) ; ierr=ierr+hferr
!          else   
!            !                print *,'File does exist, opening it RW'
!            call h5fopen_f(trim(groups(1)), H5F_ACC_RDWR_F, id(1), hferr) ; ierr=ierr+hferr
!          endif
!
!          !        print *,'Creating Compression Filter'
!          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hferr) ; ierr=ierr+hferr
!          call h5pset_deflate_f(dcpl, 9, hferr) ; ierr=ierr+hferr
!          call h5pset_chunk_f(dcpl, rank, chunk, hferr) ; ierr=ierr+hferr
!
!          do k=2,lastid
!            call h5lexists_f(id(k-1), trim(groups(k)), link_exists, hferr) ; ierr=ierr+hferr
!            if(link_exists) then
!              !                        print *,'opening group',trim(groups(k))
!              call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            else
!              !                        print *,'Create Groups',trim(groups(k))
!              call h5gcreate_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
!            endif
!          enddo
!
!          !        print *,'Creating DataSpace'
!          call h5screate_simple_f(rank, dims, dataspace, hferr) ; ierr=ierr+hferr
!          if(hferr.ne.0) return
!
!          !        print *,'Create Dataset id',id(lastid),'name ',trim(groups(lastid+1))
!          call h5dcreate_f(id(lastid), trim(groups(lastid+1)), H5T_NATIVE_DOUBLE, dataspace, dataset, hferr, dcpl) ; ierr=ierr+hferr
!
!          if(ierr.ne.0) return
!          !        print *,'Write to dataset'
!          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dble(arr), dims, hferr) ; ierr=ierr+hferr
!
!          !        print *,'Closing handles'
!          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr
!          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
!          do k=lastid,2,-1
!            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
!          enddo
!          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
!
!          !        print *,'Data is now:',sum(arr)/size(arr)
!
!          call h5close_f(hferr)  ; ierr=ierr+hferr
!          close(v,status='delete')
!          return
!          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
!          call cpusleep(waitinterval)
!        enddo
!    end subroutine

    subroutine h5load_1d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable :: arr(:)
        double precision,allocatable :: dparr(:)
        integer(iintegers),intent(out) :: ierr

        integer(HSIZE_T) :: dims(1),maxdims(1)
        logical :: file_exists,link_exists,compression
        integer :: k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete') 

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          !        do k=1,size(groups)
          !          print *,'opening hdf5 file groups(',k,') ',trim(groups(k))
          !        enddo

          compression=check_for_zip_filter()

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          !Get dataset creation properties list
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1)))
          allocate(dparr(dims(1)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)

        enddo
    end subroutine
    subroutine h5load_2d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable :: arr(:,:)
        double precision,allocatable :: dparr(:,:)
        integer(iintegers),intent(out) :: ierr

        integer(HSIZE_T) :: dims(2),maxdims(2)
        logical :: file_exists,link_exists
        integer :: k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete') 

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          !        print *,'opening hdf5 file ',trim(groups(1))

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          !Get dataset creation properties list
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2)))
          allocate(dparr(dims(1),dims(2)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5load_3d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable :: arr(:,:,:)
        double precision,allocatable :: dparr(:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer(HSIZE_T) :: dims(3),maxdims(3)
        logical :: file_exists,link_exists
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          !        print *,'opening 3d hdf5 file ',trim(groups(1))

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2),dims(3)))
          allocate(dparr(dims(1),dims(2),dims(3)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'loaded 3d data, with shape',shape(arr)
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5load_4d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable :: arr(:,:,:,:)
        double precision,allocatable :: dparr(:,:,:,:)
        integer(HSIZE_T) :: dims(4),maxdims(4)
        integer(iintegers),intent(out) :: ierr
        logical :: file_exists,link_exists,compression
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          print *,'opening 4d hdf5 file ',trim(groups(1))

          compression=check_for_zip_filter()

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2),dims(3),dims(4)))
          allocate(dparr(dims(1),dims(2),dims(3),dims(4)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          !        print *,'loaded 4d data, with shape',shape(arr)
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5load_5d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable :: arr(:,:,:,:,:)
        double precision,allocatable :: dparr(:,:,:,:,:)
        integer(iintegers),intent(out) :: ierr

        integer(HSIZE_T) :: dims(5),maxdims(5)
        logical :: file_exists,link_exists,compression
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          compression=check_for_zip_filter()

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2),dims(3),dims(4),dims(5)))
          allocate(dparr(dims(1),dims(2),dims(3),dims(4),dims(5)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5load_7d(groups,arr,ierr)
        character(len=*) :: groups(:)
        real(ireals),allocatable,dimension(:,:,:,:,:,:,:) :: arr
        double precision,allocatable,dimension(:,:,:,:,:,:,:) :: dparr
        integer(HSIZE_T) :: dims(7),maxdims(7)
        integer(iintegers),intent(out) :: ierr
        logical :: file_exists,link_exists
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
          allocate(dparr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))

          call h5dread_f(dataset, H5T_NATIVE_DOUBLE, dparr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr
          arr=real(dparr); deallocate(dparr)

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          !        print *,'loaded 4d data, with shape',shape(arr)
          call h5close_f(hferr)
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine

    subroutine h5load_1d_int(groups,out_arr,ierr)
        character(len=*) :: groups(:)
        integer(iintegers),allocatable :: out_arr(:)
        integer(iintegers),intent(out) :: ierr

        integer,allocatable :: arr(:)
        integer(HSIZE_T) :: dims(1),maxdims(1)
        logical :: file_exists,link_exists
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          !        print *,'opening hdf5 file ',trim(groups(1))

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          !Get dataset creation properties list
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1)))

          !        print *,'read',ierr
          call h5dread_f(dataset, H5T_NATIVE_INTEGER, arr, dims, hferr) ; ierr=ierr+hferr

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          call h5close_f(hferr)
          allocate(out_arr(ubound(arr,1)))
          out_arr = arr
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine h5load_2d_int(groups,out_arr,ierr)
        character(len=*) :: groups(:)
        integer(iintegers),allocatable :: out_arr(:,:)
        integer(iintegers),intent(out) :: ierr

        integer,allocatable :: arr(:,:)
        integer(HSIZE_T) :: dims(2),maxdims(2)
        logical :: file_exists,link_exists
        integer k,lastid,hferr,rank
        integer(HID_T) :: id(size(groups)-1),dataset,dataspace
        character(300) :: name

        lockfile=trim(groups(1))//'.lock'
        do iwait=1,int(maxwait/waitinterval)
          open(v,file=lockfile,status='new',err=99)
          write(v,*) 'file is locked by process: ',getpid()
          close(v,status='delete')

          call h5open_f(hferr)
          ierr=0 ; lastid = ubound(id,1)
          if(size(groups).lt.3) print *,'ARGHHH :: need at least 3 group entries, first is filename &
              &   and scnd is at least hdf5 root /, and third is data name'
          !        print *,'opening hdf5 file ',trim(groups(1))

          inquire(file=trim(groups(1)), exist=file_exists)
          if(.not.file_exists) then
            ierr=-5
            print *,'file doesnt exist: ',trim(groups(1))
            close(v,status='delete') ; return
          endif

          call h5fopen_f(trim(groups(1)), H5F_ACC_RDONLY_F, id(1), hferr) ; ierr=ierr+hferr
          !First check if everything is here:
          name = ''
          do k=2,lastid+1
            name = trim(name)//'/'//trim(groups(k))
            call h5lexists_f(id(1), trim(name), link_exists, hferr) ; ierr=ierr+hferr
            !                print *,'checked if ',trim(name),' exists: ',link_exists,ierr
            if (.not.link_exists) then
              ierr=-6
              call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
              close(v,status='delete') ; return
            endif
          enddo

          do k=2,lastid
            call h5gopen_f(id(k-1), trim(groups(k)), id(k), hferr) ; ierr=ierr+hferr
          enddo
          !Get dataset creation properties list
          call h5dopen_f(id(lastid), trim(groups(lastid+1)), dataset, hferr) ; ierr=ierr+hferr
          !        print *,'groups all open',ierr

          call h5dget_space_f(dataset, dataspace, hferr) ; ierr=ierr+hferr
          call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, rank) ; ierr=ierr+min(0,rank)
          call h5sclose_f(dataspace,hferr) ; ierr=ierr+hferr

          allocate(arr(dims(1),dims(2)))

          call h5dread_f(dataset, H5T_NATIVE_INTEGER, arr, dims, hferr) ; ierr=ierr+hferr
          !        print *,'read',ierr

          call h5dclose_f(dataset,hferr) ; ierr=ierr+hferr
          do k=lastid,2,-1
            call h5gclose_f(id(k),hferr) ; ierr=ierr+hferr
          enddo
          call h5fclose_f(id(1),hferr) ; ierr=ierr+hferr
          !        print *,'closed',ierr

          !        print *,'Data average is now:',sum(arr)/size(arr)
          call h5close_f(hferr)

          allocate(out_arr(ubound(arr,1),ubound(arr,2) ))
          out_arr = arr
          return
          99 if(rand() .gt. .999)     write(*,*) 'lock file already exists: ',lockfile
          call cpusleep(waitinterval)
        enddo
    end subroutine
    subroutine cpusleep(sec)
        real(ireals) :: sec
        real(ireals) :: t1,t2,dither
        call random_number(dither)
        dither = ( -1.0_ireals + 2*dither ) *.1*sec !10% dithering on time
        call cpu_time(t1)
        do
          call cpu_time(t2)
          if(t2-t1.gt.dither) return
        enddo
    end subroutine
    function check_for_zip_filter()
        logical :: check_for_zip_filter
        integer :: hferr

        CALL h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, check_for_zip_filter, hferr)
        if (.not.check_for_zip_filter) then
          write(*,'("szip compression filter not available. :( -- I might not be able to load some data... reconfigure the HDF5 Installation",/)')
          !            call exit()
        else
          !           print *,'Compression filter is ready'
        endif
    end function

  end module

  !program main
  !      use arrayIO
  !      implicit none
  !
  !      integer :: r(10,2,6)
  !      integer,allocatable :: o(:,:,:)
  !
  !      character(100) :: fname='test3d'
  !
  !
  !      integer i,j
  !
  !      do i=1,2
  !        do j=1,10
  !          r(j,i,:) = 10*(i-1)+j
  !        enddo
  !      enddo
  !
  !      call write_bin(fname,r)
  !      r = 0
  !      call read_bin(fname,o)
  !
  !      print *,o
  !
  !end program


