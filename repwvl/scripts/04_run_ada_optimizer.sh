#!/bin/bash

set -euo pipefail

. config.sh

MPIEXEC=${MPIEXEC:-srun -n1 -c 40 --mem=60G --time=12:00:00}
SOLVER=basinhopper
VMODE=1

OPTIFILE=$WORK/opti.$(basename $REPWVL_LW)

if [ ! -e $OPTIFILE ]; then
  ${MPIEXEC} python $ADA/optipy/opti.py \
    --Nsamples $REPWVL_SAMPLES_LW \
    --maxiter 10000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_LW \
    $OPTIFILE
  rm -f $REPWVL_LW # invalidate annotated file
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
    --maxiter 10000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_SW \
    $OPTIFILE
  rm -f $REPWVL_SW # invalidate annotated file
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
