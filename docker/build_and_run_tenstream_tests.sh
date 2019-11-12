#!/bin/bash

. $WORKDIR/.profile

mkdir -p $WORKDIR/TenstreamLUT_cpy
cd $WORKDIR/TenstreamLUT_cpy
cp -v -rns $WORKDIR/TenstreamLUT/*.mmap4 $WORKDIR/TenstreamLUT_cpy
echo "-lut_basename $(pwd)/LUT" >> $HOME/.petscrc

cd $WORKDIR/tenstream || exit 1
rm -rf $WORKDIR/tenstream/build
mkdir build && cd build && cmake .. -DSYST=$SYST -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DNETCDF_DIR=$NETCDF_DIR -DCTEST_MPIRUN_FLAGS="--allow-run-as-root" || exit
make -j all || exit 2

make check | tee ctest_log

FAILED_TESTS=$(grep '...\*\*\*Failed' ctest_log | awk {'print $4'})
for failed_test in $FAILED_TESTS
do
  echo "Re-Running verbosely failed test $failed_test"
  ctest -V -R $failed_test
done
[ -z "$FAILED_TESTS" ] && (echo "Tests seem OK") || (echo "Some Tests failed"; exit 3)
