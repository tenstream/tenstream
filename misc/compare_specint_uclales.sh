#!/bin/bash

set -euo pipefail

BIN=ex_uclales_cld_file
CLD=/project/meteo/public/jakub/phillip1_0d_25m_srfc/phillip1_0d_25m_srfc.merged.nc

make -j -C $TENSTREAMROOT $BIN

MPIEXEC=${MPIEXEC:-srun -n 64 --mem=120G --time=08:00:00}

for specint in rrtmg repwvl; do
for solver in 2str 3_10; do
  ID=$(basename $CLD .nc).$specint.$solver
  OUT=$SCRATCH/compare_specint_uclales/$ID.nc
  LOG=$SCRATCH/compare_specint_uclales/$ID.log
  mkdir -p $(dirname $OUT)
  ${MPIEXEC} \
    $TENSTREAMROOT/bin/$BIN \
    -cld $CLD \
    -specint $specint \
    -solver $solver \
    -tstart 800 \
    -tend 820 \
    -log_view \
    -out $OUT \
    &> >(tee $LOG)
done
done
