#!/bin/bash
set -euo pipefail

echo "Using TENSTREAMROOT: $TENSTREAMROOT"
make -j -C $TENSTREAMROOT ex_pprts_specint_lw_sw

LATEST_DIFFUSE=$(ls -t S.*.nc| head -1)
echo "Using DIFFUSE model: $LATEST_DIFFUSE"

function run () {
  LOG=$1
  OPT=$2
  ANN=$3

  if [ ! -z "$ANN" ]; then
    ln -sf $(readlink -f $ANN) /tmp/ANN_dir2dir_ANN_.tau31.w020.aspect_zx23.g6.phi19.theta19.c3.nc
    ln -sf $(readlink -f $ANN) /tmp/ANN_dir2diff_ANN_.tau31.w020.aspect_zx23.g6.phi19.theta19.c3_10.nc
    ln -sf $(readlink -f $ANN) /tmp/ANN_diff2diff_ANN_.tau31.w020.aspect_zx23.g6.c10.nc
  fi

  BASEOPT="-solar no -Nz 32 -specint ecckd -pprts_view_geometry -solver 3_10 -ann_basename /tmp/ANN"

  if [ -e $LOG ]; then
    echo "Skipping $LOG"
  else
    echo "Running $LOG :: $OPT"
    $TENSTREAMROOT/bin/ex_pprts_specint_lw_sw $BASEOPT $OPT &> >(tee $LOG) || echo "failed"
  fi
}

LOGBASE="log"
run "${LOGBASE}.ann"  "-pprts_use_ANN" "$LATEST_DIFFUSE"
run "${LOGBASE}.3_10" "" ""
run "${LOGBASE}.mc"   "-solver mcdmda" ""
for ann in S.*.nc; do
  run "${LOGBASE}.ann.$ann"  "-pprts_use_ANN" "$ann"
done

for f in ${LOGBASE}*; do
  echo "** $f"
  ack "\d+ edir" $f|head -n 4||echo "no edir"
  echo "..."
  ack "\d+ edir" $f|tail -n 4||echo "no edir"
done
