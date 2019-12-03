# LRZ Linux Cluster
# before, load modules:
#
# module purge
#
# module load lrz/nodev cmake git
# module load intel/19.1
# module load mpi.intel/2019
# module load mkl/2019_s

set(CMAKE_Fortran_COMPILER   "mpiifort")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g")
set(USER_Fortran_FLAGS_RELEASE "-march=native -O3 -ftz -fp-model fast=2 -no-prec-div -no-prec-sqrt -fast-transcendentals")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(USER_MPIEXEC "mpiexec")
set(PETSC_SKIP_TEST_RUNS True)
