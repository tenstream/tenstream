#!/bin/bash

set -euo pipefail

. config.sh

MAXTIME=24 # hours

MPIEXEC=${MPIEXEC:-srun -n1 -c 1 --mem=200G --time=${MAXTIME}:00:00}
SOLVER=basinhopper
VMODE=1

OPTIFILE=$WORK/opti.$(basename $REPWVL_LW)
LOG=${OPTIFILE}.stdout
rm -f $LOG

PIDS=()
(
if [ ! -e $OPTIFILE ]; then
  echo "Running optimizer for $OPTIFILE" &> >(tee -a $LOG)
  ${MPIEXEC} python $ADA/pyADA/pyADA.py \
    --Nsamples $REPWVL_SAMPLES_LW \
    --maxiter 1000 \
    --maxtime $((${MAXTIME}*3600-900)) \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_LW \
    $OPTIFILE \
    &> >(tee -a $LOG)
  rm -f $REPWVL_LW # invalidate annotated file
else
  echo "Skipping optimizer for $OPTIFILE" &> >(tee -a $LOG)
fi

if [ ! -e $REPWVL_LW ]; then
  echo "Running merge for $REPWVL_LW"
  python $ADA/pyADA/merge2repwvl.py \
    --thermal \
    --verbose \
    $LBL_LW \
    $OPTIFILE \
    $REPWVL_LW \
    &> >(tee -a $LOG)
else
  echo "Skipping merge for $REPWVL_LW" &> >(tee -a $LOG)
fi
) &
pids[0]=$!


OPTIFILE=$WORK/opti.$(basename $REPWVL_SW)
LOG=${OPTIFILE}.stdout

(
if [ ! -e $OPTIFILE ]; then
  echo "Running optimizer for $OPTIFILE" &> >(tee -a $LOG)
  ${MPIEXEC} python $ADA/pyADA/pyADA.py \
    --Nsamples $REPWVL_SAMPLES_SW \
    --maxiter 1000 \
    --maxtime $((${MAXTIME}*3600-900)) \
    --verbose \
    --solver $SOLVER \
    --var_mode $VMODE \
    $TRAIN_SW \
    $OPTIFILE \
    &> >(tee -a $LOG)
  rm -f $REPWVL_SW # invalidate annotated file
else
  echo "Skipping optimizer for $OPTIFILE" &> >(tee -a $LOG)
fi

if [ ! -e $REPWVL_SW ]; then
  echo "Running merge for $REPWVL_SW" &> >(tee -a $LOG)
  python $ADA/pyADA/merge2repwvl.py \
    --verbose \
    --crsNO2 $ADA/lkp_annotate/crs_NO2_UBremen_cf.nc \
    --crsO3  $ADA/lkp_annotate/crs_O3_UBremen_cf.nc \
    --solar  $ADA/lkp_annotate/solar_spectrum.nc \
    $LBL_SW \
    $OPTIFILE \
    $REPWVL_SW \
    &> >(tee -a $LOG)
  echo "Skipping merge for $REPWVL_SW" &> >(tee -a $LOG)
fi
) &
pids[1]=$!

trap 'kill ${pids}' EXIT
for pid in ${pids[*]}; do echo "waiting for job $pid to finish"; wait $pid; done
