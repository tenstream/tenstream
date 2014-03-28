module data_parameters
      integer,parameter :: &
      ireals=selected_real_kind(15,200), &
      iintegers=selected_int_kind(12)
      real(ireals),parameter :: pi=3.141592653589793
      real(ireals),parameter :: zero=0._ireals,one=1._ireals
      integer(iintegers) ,parameter :: dir_streams=8, diff_streams=10
      integer(iintegers) ,parameter :: i0=0,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9,i10=10,i11=11

end module
