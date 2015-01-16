# ARCH Linux
set(CMAKE_C_COMPILER "mpcc_r")
set(CMAKE_CXX_COMPILER "mpCC_r")
set(CMAKE_Fortran_COMPILER "mpxlf2003_r")
set(Fortran_COMPILER_WRAPPER mpxlf2003_r)

set(USER_Fortran_FLAGS "-qextname -qfree=F90 -qwarn64 -qnosave -qinitauto=FFF00000 -qflag=w:e -qsuffix=cpp=f90 ")
set(USER_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions -qinitauto=ff -qsigtrap")
set(USER_Fortran_FLAGS_DEBUG "-O0 -qfullpath -C -g -qflttrp=enable:inexact:invalid:nanq:overflow:zerodivide -qsigtrap -qinitauto")

set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.2.1.1/include")
set(NETCDF_LIB_1       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.8/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.8/lib/libhdf5.a")
set(SZIP_LIB           "/sw/aix53/szip-2.1/lib/libsz.a")
include_directories(${NETCDF_INCLUDE_DIR})
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z pessl dl )

