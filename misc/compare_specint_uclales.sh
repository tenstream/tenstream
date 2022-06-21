#!/bin/bash

set -euo pipefail

BIN=ex_uclales_cld_file
CLD=/project/meteo/public/jakub/phillip1_0d_25m_srfc/phillip1_0d_25m_srfc.merged.nc

make -j -C $TENSTREAMROOT $BIN

MPIEXEC=${MPIEXEC:-srun -n 64 --mem=120G --time=08:00:00}

TSTART=800
TEND=810
TINC=1

DT=10

for specint in repwvl rrtmg
do
for solver in 2str 3_10 3_10_incomplete
do
  ID=$(basename $CLD .nc).$specint.$solver
  OUT=$SCRATCH/compare_specint_uclales/$ID.nc
  LOG=$SCRATCH/compare_specint_uclales/$ID.log

  OPT=""
  if [[ $solver == *"incomplete"* ]]; then
    OPT="$OPT \
      -solver 3_10 \
      -solar_dir_explicit \
      -solar_dir_ksp_view \
      -solar_dir_ksp_max_it 1 \
      -solar_dir_ksp_skip_residual \
      -solar_dir_ksp_ignore_max_it $(($TSTART * $DT)) \
      -solar_diff_explicit \
      -solar_diff_ksp_view \
      -solar_diff_ksp_max_it 1 \
      -solar_diff_ksp_skip_residual \
      -solar_diff_ksp_ignore_max_it $(($TSTART * $DT)) \
      -thermal_diff_explicit \
      -thermal_diff_ksp_view \
      -thermal_diff_ksp_max_it 1 \
      -thermal_diff_ksp_skip_residual \
      -thermal_diff_ksp_ignore_max_it $(($TSTART * $DT)) \
      -accept_incomplete_solve \
      -absorption_by_coeff_divergence \
    "
  fi

  mkdir -p $(dirname $OUT)

  if [ ! -e $OUT ]; then
  ${MPIEXEC} \
    $TENSTREAMROOT/bin/$BIN \
    -cld $CLD \
    -specint $specint \
    -solver $solver \
    -tstart $TSTART \
    -tend $TEND \
    -tinc $TINC \
    -log_view \
    -out $OUT \
    $OPT \
    &> >(tee $LOG)
  else
    echo "Skipping $ID because output $OUT already exists"
  fi
done
done
