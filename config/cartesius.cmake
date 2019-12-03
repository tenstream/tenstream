# Cartesius supercomputer at SARA Amsterdam
# modules loaded:

# module unload compilerwrappers
# module load hdf5 netcdf cmake python git
# module load blas mkl

# Petsc installed with:
# module unload compilerwrappers
# base_opt="--with-fortran --with-fortran-interfaces --with-shared-libraries=1 --with-blas-lapack-dir=$SURFSARA_MKL_ROOT"
# intel_compilers="--with-cc=$(which mpiicc) --with-fc=$(which mpiifort) --with-cxx=$(which mpiicpc)"
# export PETSC_ARCH='prod_icc'
# ./configure $intel_compilers $base_opt $batch_compile_opts \
# '--with-debugging=0' COPTFLAGS='-mkl' FOPTFLAGS='-mkl' && make

set(NETCDF_DIR      "$ENV{SURFSARA_NETCDF_ROOT}")
set(NETCDF_DIR_F90  "$ENV{SURFSARA_NETCDF_ROOT}")
set(BLA_PREFER_PKGCONFIG "$ENV{SURFSARA_MKL_ROOT}")

set(USER_C_FLAGS "-nowarn -std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend-source -g -mkl ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xCORE-AVX2 -fp-model source -fno-omit-frame-pointer")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(CMAKE_C_COMPILER   "mpiicc")
set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")
