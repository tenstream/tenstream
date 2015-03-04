# ARCH Linux
set(CMAKE_C_COMPILER "mpcc_r")
set(CMAKE_CXX_COMPILER "mpCC_r")
set(CMAKE_Fortran_COMPILER "mpxlf2003_r")
set(Fortran_COMPILER_WRAPPER mpxlf2003_r)

set(USER_Fortran_FLAGS "-qextname -qfree=F90 -qwarn64 -qnosave -qflag=w:e -qsuffix=cpp=f90 -qxlf90=autodealloc -q64 -qinitauto=FFF00000 -qflttrap=en:ov:zero:inv:imp ")
set(USER_Fortran_FLAGS_RELEASE "-O4 -qtune=pwr6 -qarch=pwr6 -qhot -qstrict=none:exceptions -qinitauto=ff -qsigtrap")
set(USER_Fortran_FLAGS_DEBUG "-O0 -qfullpath -C -g9 -qinitauto=FFF00000 -qflttrp=enable:inexact:invalid:nanq:overflow:zerodivide -qsigtrap -qinitauto")

set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.2.1.1/include")
set(NETCDF_LIB_1       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.8/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.8/lib/libhdf5.a")
set(SZIP_LIB           "/sw/aix53/szip-2.1/lib/libsz.a")
set(LAPACK_LIB         "/sw/aix53/lapack-3.2.0/lib/liblapack.a")
include_directories(${NETCDF_INCLUDE_DIR})
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} ${LAPACK_LIB} essl blas mass z dl )
message(STATUS "USERLIBS: ${LIBS}")

include_directories("/sw/aix53/valgrind-3.3.0/include")
