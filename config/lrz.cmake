# LRZ Linux Cluster
# before, load modules:
#
# module load ccomp/intel
# module load fortran/intel
# module load gcc
# module load netcdf

set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xHOST -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")


set(NETCDF_DIR "$ENV{NETCDF_BASE}")

