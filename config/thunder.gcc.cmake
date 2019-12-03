# Thunder configuration with custom petsc and hdf5 install in Home
# module setting:
#   module switch openmpi openmpi/1.6.5-static-gcc48
#   module switch gcc gcc/4.8.2
#

set(CMAKE_C_COMPILER       "mpicc")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS               " -cpp --std=c99 ")
set(USER_Fortran_FLAGS         " -cpp -fbacktrace -finit-real=nan -ffree-line-length-none -no-wrap-margin ")
set(USER_Fortran_FLAGS_RELEASE " -funroll-all-loops -O3 -march=native -mtune=native ")
set(USER_Fortran_FLAGS_DEBUG   " -W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow ")

set(NETCDF_DIR      "/sw/squeeze-x64/netcdf-4.2-static")
set(NETCDF_DIR_F90  "/sw/squeeze-x64/netcdf_fortran-4.2-static-gcc48")

set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
set(LIBS ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
