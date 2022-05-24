#!/bin/bash

set -euo pipefail

. config.sh

MPIEXEC=${MPIEXEC:-srun -n1 -c 40 --mem=140G --time=24:00:00}
SOLVER=basinhopper
VMODE=1

OPTIFILE=$WORK/opti.$(basename $REPWVL_LW)

if [ ! -e $OPTIFILE ]; then
  echo "Running optimizer for $OPTIFILE"
  ${MPIEXEC} python $ADA/pyADA/pyADA.py \
    --Nsamples $REPWVL_SAMPLES_LW \
    --maxiter 1000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_LW \
    $OPTIFILE
  rm -f $REPWVL_LW # invalidate annotated file
else
  echo "Skipping optimizer for $OPTIFILE"
fi

if [ ! -e $REPWVL_LW ]; then
  echo "Running merge for $REPWVL_LW"
  python $ADA/pyADA/merge2repwvl.py \
    --thermal \
    --verbose \
    $LBL_LW \
    $OPTIFILE \
    $REPWVL_LW
else
  echo "Skipping merge for $REPWVL_LW"
fi


OPTIFILE=$WORK/opti.$(basename $REPWVL_SW)

if [ ! -e $OPTIFILE ]; then
  echo "Running optimizer for $OPTIFILE"
  ${MPIEXEC} python $ADA/pyADA/pyAD.py \
    --Nsamples $REPWVL_SAMPLES_SW \
    --maxiter 1000 \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_SW \
    $OPTIFILE
  rm -f $REPWVL_SW # invalidate annotated file
else
  echo "Skipping optimizer for $OPTIFILE"
fi

if [ ! -e $REPWVL_SW ]; then
  echo "Running merge for $REPWVL_SW"
  python $ADA/pyADA/merge2repwvl.py \
    --verbose \
    --crsNO2 $ADA/lkp_annotate/crs_NO2_UBremen_cf.nc \
    --crsO3  $ADA/lkp_annotate/crs_O3_UBremen_cf.nc \
    --solar  $ADA/lkp_annotate/solar_spectrum.nc \
    $LBL_SW \
    $OPTIFILE \
    $REPWVL_SW
  echo "Skipping merge for $REPWVL_SW"
fi
