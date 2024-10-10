#!/bin/bash
set -euo pipefail

echo "Using TENSTREAMROOT: $TENSTREAMROOT"

LATEST_DIFFUSE=$(ls -t S.*.nc| head -1)
echo "Using DIFFUSE model: $LATEST_DIFFUSE"

LOGBASE="log"
BIN=$TENSTREAMROOT/bin/ex_pprts_specint_lw_sw
BASEOPT="-solar no -Nz 140 -dx 100 -specint ecckd -pprts_view_geometry -solver 3_10 -ann_basename $SCRATCH/tmp/ANN"
SOPT=""

LOGBASE="li3rc"
CLD=$HOME/tenstream/examples/libRadtran_cld_file/i3rc_cloudfile.nc
BIN=$TENSTREAMROOT/bin/ex_uvspec_cld_file
MCOPT="-mcdmda_queue_size 1000000 -mcdmda_photons_per_px 200 -mcdmda_info -mcdmda_info_advance -mcdmda_batch_size 10000"
BASEOPT="-wc $CLD \
  -specint ecckd \
  -phi 200 -theta 40 \
  -solar no \
  -pprts_view_geometry \
  -solver 3_10 \
  -thermal_diff_explicit no \
  -thermal_diff_ksp_monitor \
  -ann_basename $SCRATCH/tmp/ANN \
  -log_view \
  $MCOPT\
  "
SOPT="srun -N 1 -n 64 --mem=128G --time=24:00:00 -p cluster"

make -j -C $TENSTREAMROOT $(basename $BIN)

function run () {
  LOG=$1
  OPT=$2
  ANN=$3

  if [ ! -z "$ANN" ]; then
    ln -sf $(readlink -f $ANN) $SCRATCH/tmp/ANN_dir2dir_ANN_.tau31.w020.aspect_zx23.g6.phi19.theta19.c3.nc
    ln -sf $(readlink -f $ANN) $SCRATCH/tmp/ANN_dir2diff_ANN_.tau31.w020.aspect_zx23.g6.phi19.theta19.c3_10.nc
    ln -sf $(readlink -f $ANN) $SCRATCH/tmp/ANN_diff2diff_ANN_.tau31.w020.aspect_zx23.g6.c10.nc
  fi


  if [ -e $LOG ]; then
    echo "Skipping $LOG"
  else
    echo "Running $LOG :: $OPT"
    $SOPT $BIN $BASEOPT $OPT -out data.$LOG &> >(tee $LOG) || echo "failed"
  fi
}

cat > cmp.py << EOF
import xarray as xr
import sys
import numpy as np
M = xr.open_dataset('data.${LOGBASE}.mc.merged', engine='netcdf4') #.mean(['ny','nx'])

klev = int((M.zlev<1e4).argmax())
klay = int((M.zlay<1e4).argmax())

D = xr.open_dataset(sys.argv[1], engine='netcdf4') #.mean(['ny','nx'])
D, M = [ _.isel(zlev=slice(klev,None)).isel(zlay=slice(klay,None)) for _ in (D, M) ]

rrmse = lambda a,b: float(np.sqrt(np.mean((a-b)**2)) / np.mean(b) * 100)
bias = lambda a,b: float((np.mean(a) / np.mean(b) - 1) * 100)
err_edn, bias_edn = rrmse(D.edn, M.edn), bias(D.edn, M.edn)
err_eup, bias_eup = rrmse(D.eup, M.eup), bias(D.eup, M.eup)
err_abso, bias_abso = rrmse(D.abso, M.abso), bias(D.abso, M.abso)
print(f"{sys.argv[1]:40s} :: edn {err_edn:5.2f}% ({bias_edn:5.2f}%) eup {err_eup:5.2f}% ({bias_eup:5.2f}%) abso {err_abso:5.2f}% ({bias_abso:5.2f}%)")
EOF

run "${LOGBASE}.2str"   "-solver 2str -schwarzschild no" ""
run "${LOGBASE}.disort" "-solver disort" ""
for i in {01..09}; do
  run "${LOGBASE}.mc${i}"     "-solver mcdmda" ""
done

#rm -f data.li3rc.mc.merged
if [ ! -e data.li3rc.mc.merged ]; then
  python -c "import xarray as xr; xr.open_mfdataset('data.li3rc.mc*',engine='netcdf4', concat_dim='iter', combine='nested').mean('iter').to_netcdf('data.li3rc.mc.merged')"
fi

run "${LOGBASE}.3_10"   "" ""
run "${LOGBASE}.ann"    "-pprts_use_ANN" "$LATEST_DIFFUSE"
for ann in S.*.nc; do
  run "${LOGBASE}.ann.$ann"  "-pprts_use_ANN" "$ann"
done

for f in ${LOGBASE}*; do
  echo "** $f"
  ack "\d+ edir" $f|head -n 4||echo "no edir"
  echo "..."
  ack "\d+ edir" $f|tail -n 4||echo "no edir"
done

for f in ${LOGBASE}*; do
  python cmp.py data.$f
done
