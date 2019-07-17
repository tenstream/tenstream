# LRZ Linux Cluster
# before, load modules:
#
# module unload mkl mpi.intel intel
#
# module load intel
# module load mpi.intel
# module load mkl
# module load netcdf/4.6.1-intel-impi-hdf5v1.10-parallel
# module load netcdf-fortran/4.4.4-intel-impi-hdf5v1.10

set(CMAKE_C_COMPILER   "mpiicc")
set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(USER_MPIEXEC "srun")

set(NETCDF_DIR "$ENV{NETCDF_BASE}")
set(NETCDF_DIR_F90 "$ENV{NETCDF_FORTRAN_BASE}")
