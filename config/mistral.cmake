# Mistral supercomputer at DKRZ Hamburg
# modules loaded:

# module add intel/16.0.2 mxm/3.3.3002 fca/2.5.2393 bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
#
# export OMPI_MCA_pml=cm         # sets the point-to-point management layer
# export OMPI_MCA_mtl=mxm        # sets the matching transport layer (MPI-2 one-sided comm.)
# export MXM_RDMA_PORTS=mlx5_0:1
# export OMPI_MCA_coll=^ghc             # disable BULLs GHC algorithm for collectives
# export OMPI_MCA_coll_fca_priority=95
# export OMPI_MCA_coll_fca_enable=1
#
# export MPIMODE=bullxmpi
# export MYmpif90=mpif90
# export MYmpicc=mpicc
# export MYmpicxx=mpicxx
#
# # NETCDF LIBPATH
# export NETCDFCROOT=/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-$MPIMODE-intel14/
# export NETCDFFROOT=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-$MPIMODE-intel14/
# export HDF5ROOT=/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-$MPIMODE-intel14/
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFFROOT/lib:$NETCDFCROOT/lib:
# export PATH=$PATH:$NETCDFFROOT/bin:$NETCDFCROOT/bin
# 
# module load cmake

set(NETCDF_DIR      "$ENV{NETCDFCROOT}")
set(NETCDF_DIR_F90  "$ENV{NETCDFFROOT}")
set(HDF5_DIR        "$ENV{HDF5ROOT}")

# Petsc installed with:
# export PETSC_ARCH=fast-$MPIMODE
#./configure                           \
#  --with-cc=$(which  $MYmpicc)        \
#  --with-fc=$(which  $MYmpif90)       \
#  --with-cxx=$(which $MYmpicxx)       \
#  --with-shared-libraries=0           \
#  --with-netcdf-dir=$NETCDFCROOT      \
#  --with-hdf5-dir=$HDF5ROOT           \
#  --with-blas-lapack-dir=$MKL         \
#  --with-fortran                      \
#  --with-fortran-interfaces           \
#  --with-cmake=$(which cmake)         \
#  --with-precision=single             \
#  --download-metis                    \
#  --download-parmetis                 \
#  --with-debugging=0                  \
#  --with-valgrind-dir=/usr            \
#  COPTFLAGS='-O3 -xAVX -mkl' \
#  FOPTFLAGS='-O3 -xAVX -mkl' \
#  && make all test install

set(CMAKE_C_COMPILER       "$ENV{MYmpicc}")
set(CMAKE_CXX_COMPILER     "$ENV{MYmpicxx}")
set(CMAKE_Fortran_COMPILER "$ENV{MYmpif90}")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g -sox -no-wrap-margin ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -xCORE-AVX2 -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")
