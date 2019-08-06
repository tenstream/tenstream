#!/bin/bash

. $WORKDIR/.profile

mkdir -p $WORKDIR/TenstreamLUT_cpy
cd $WORKDIR/TenstreamLUT_cpy
cp -v $WORKDIR/TenstreamLUT/*.mmap4 $WORKDIR/TenstreamLUT_cpy
echo "-lut_basename $(pwd)/LUT" >> $HOME/.petscrc

cd $WORKDIR/tenstream || exit
rm -rf $WORKDIR/tenstream/build
mkdir build && cd build && cmake .. -DSYST=$SYST -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DNETCDF_DIR=$NETCDF_DIR -DCTEST_MPIRUN_FLAGS="--allow-run-as-root" || exit
make -j all || exit
make check || exit
