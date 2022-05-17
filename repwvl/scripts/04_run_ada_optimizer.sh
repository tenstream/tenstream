#!/bin/bash

set -euo pipefail

. config.sh

MPIEXEC=${MPIEXEC:-srun -n1 -c 20 -C $MODULES_MARCH -p met-cl,met-ws,cluster --mem=40G}
SOLVER=basinhopper
VMODE=1

OPTIFILE=$WORK/opti.$(basename $REPWVL_LW)

if [ ! -e $OPTIFILE ]; then
  ${MPIEXEC} python $ADA/optipy/opti.py \
    --Nsamples $REPWVL_SAMPLES_LW \
    --maxiter 1000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_LW \
    $OPTIFILE
fi

if [ ! -e $REPWVL_LW ]; then
  python $ADA/optipy/merge2repwvl.py \
    --thermal \
    --verbose \
    $LBL_LW \
    $OPTIFILE \
    $REPWVL_LW
fi


OPTIFILE=$WORK/opti.$(basename $REPWVL_SW)

if [ ! -e $OPTIFILE ]; then
  ${MPIEXEC} python $ADA/optipy/opti.py \
    --Nsamples $REPWVL_SAMPLES_SW \
    --maxiter 1000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_SW \
    $OPTIFILE
fi

if [ ! -e $REPWVL_SW ]; then
  python $ADA/optipy/merge2repwvl.py \
    --verbose \
    --crsNO2 $ADA/lkp_annotate/crs_NO2_UBremen_cf.nc \
    --crsO3  $ADA/lkp_annotate/crs_O3_UBremen_cf.nc \
    --solar  $ADA/lkp_annotate/solar_spectrum.nc \
    $LBL_SW \
    $OPTIFILE \
    $REPWVL_SW
fi
