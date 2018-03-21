# Default GCC
#
# Config script for the linux environment at the Meteorological Institute Munich 08.Jan.2015
#
# and use this config file with `cmake <tenstream_root_dir> -DSYST:STRING=lmu_mim`

set(CMAKE_C_COMPILER   "mpiicc")
set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")

set(USER_C_FLAGS       "-W -std=c99 ")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g -sox -no-wrap-margin ")
set(USER_Fortran_FLAGS_RELEASE "-O3 -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-warn all -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")

set(NETCDF_DIR      "$ENV{NETCDF}")
set(NETCDF_DIR_F90  "$ENV{NETCDF}")
#set(HDF5_DIR        "$ENV{HDF5ROOT}")
