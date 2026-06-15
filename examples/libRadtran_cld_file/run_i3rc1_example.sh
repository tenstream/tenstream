#!/bin/bash
BIN=bin/ex_uvspec_cld_file
make -j -C $(dirname $BIN)/.. $(basename $BIN) || exit 1

RUN="mpirun -np 8"
TENSTREAM_SRC=$HOME/tenstream
CLDFILE=${TENSTREAM_SRC}/examples/libRadtran_cld_file/i3rc_cloudfile.nc
ATMFILE=${TENSTREAM_SRC}/examples/libRadtran_cld_file/afglus_100m.dat
OUTDIR=${SCRATCH:-.}

BASE_OPT="-thermal no -specint ecckd $@"

OUT=${OUTDIR}/example_i3rc1.nc
OPT=$BASE_OPT
if [ ! -e $OUT ]; then
  $RUN $BIN -wc $CLDFILE -atm_filename $ATMFILE -out $OUT $OPT
else
  echo "Skipping simulation because output already exists: $OUT"
fi

## With RayLi MonteCarlo solver
#NUMTHREADS=10
#RUN="mpirun -np 1 --cpus-per-rank $NUMTHREADS"
#OUT=$OUTDIR/example_i3rc1_rayli.nc
#OPT="$BASE_OPT -solver rayli -rayli_cyclic_bc -rayli_nthreads $NUMTHREADS"
#[ ! -e $OUT ] && $RUN $BIN -wc $CLDFILE -atm_filename $ATMFILE -out $OUT $OPT
