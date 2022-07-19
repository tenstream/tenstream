#!/bin/bash

set -euo pipefail

. config.sh

CMP=$WORK/compare_repwvl_to_lbl.py
cat > $CMP <<EOF
#!/usr/bin/env python3
import xarray as xr
import numpy as np
import sys

tc = dict(
    HEADER = '\033[95m',
    OKBLUE = '\033[94m',
    OKCYAN = '\033[96m',
    OKGREEN = '\033[92m',
    WARNING = '\033[93m',
    FAIL = '\033[91m',
    ENDC = '\033[0m',
    BOLD = '\033[1m',
    UNDERLINE = '\033[4m',
    )

rmse = lambda a, b: np.sqrt(np.mean((a - b) ** 2))
rel_rmse = lambda a, b: rmse(a, b) / np.abs(np.mean(b)) * 100
bias = lambda a, b: np.mean(a) - np.mean(b)
rel_bias = lambda a, b: (np.mean(a) / np.mean(b) - 1.0) * 100

print(f"Comparing {tc['HEADER']} {sys.argv[2]:>100} {tc['ENDC']}\nto benchmark {tc['HEADER']} {sys.argv[1]:>100} {tc['ENDC']}")
lbl = xr.open_dataset(sys.argv[1])
rep = xr.open_dataset(sys.argv[2])

print("Mean errors in total atmosphere")
for k in lbl.keys():
    err, rel, b, rb = [ _(rep[k].data, lbl[k].data) for _ in (rmse, rel_rmse, bias, rel_bias) ]
    print(f"{k:10s} lbl {lbl[k].mean().data:12.4g} rep {rep[k].mean().data:12.4f} bias {b:12.6f} ({rb:8.2f} %)")

print("Absorption errors along vertical layers")
for k in lbl.nlay:
    err, rel, b, rb = [ _(rep['abso'].isel(nlay=k).data, lbl['abso'].isel(nlay=k).data) for _ in (rmse, rel_rmse, bias, rel_bias) ]
    print(f"abso lay {k.data:3} lbl {lbl['abso'].isel(nlay=k).mean().data:12.4g} rep {rep['abso'].isel(nlay=k).mean().data:12.4f} bias {b:12.6f} ({rb:8.2f} %)")

print("Surface Irrdiance errors")
for v in ('edir', 'edn', 'eup'):
    if v in lbl.keys():
        err, rel, b, rb = [ _(rep[v].isel(nlev=-1).data, lbl[v].isel(nlev=-1).data) for _ in (rmse, rel_rmse, bias, rel_bias) ]
        print(f"{v:6s} srfce lbl {lbl[v].isel(nlev=-1).mean().data:12.4g} rep {rep[v].isel(nlev=-1).mean().data:12.4f} bias {b:12.6f} ({rb:8.2f} %)")
EOF


BIN=$TENSTREAM_ROOT/bin/ex_pprts_repwvl_lw_sw
if [ ! -e $BIN ]; then
  echo "Building $BIN"
  make -j -C $TENSTREAM_ROOT $(basename $BIN)
fi


MIE=$MIE_REPWVL
if [ ! -e $MIE_REPWVL ]; then
  MIE=$MIE_GENERAL
fi

FU=$FU_REPWVL
if [ ! -e $FU_REPWVL ]; then
  FU=$FU_GENERAL
fi

LOG_PREFIX=$WORK/log.$(basename $BIN).$REPWVL_SAMPLES_SW.$REPWVL_SAMPLES_LW

############
BASEOPT="-solver 2str -lwc .1 -iwc .1"
SOLAR="yes"
THERMAL="no"
LOG=${LOG_PREFIX}_SW${SOLAR}_LW${THERMAL}
if [ ! -e $LOG.lbl.nc ]; then
$BIN \
  $BASEOPT \
  -solar   $SOLAR   \
  -thermal $THERMAL \
  -repwvl_data_thermal $ADA/repwvl/rdata.lbl.thermal.nc  \
  -repwvl_data_solar   $ADA/repwvl/rdata.lbl.solar.nc    \
  -mie_wc $MIE_GENERAL \
  -fu_ice $FU_GENERAL \
  -out $LOG.lbl.nc \
  | tee $LOG.lbl
else
  echo "Skipping example computation because output $LOG.lbl.nc exists"
fi

if [ ! -e $LOG.repwvl.nc ]; then
$BIN \
  $BASEOPT \
  -solar   $SOLAR   \
  -thermal $THERMAL \
  -repwvl_data_thermal $REPWVL_LW \
  -repwvl_data_solar $REPWVL_SW \
  -mie_wc $MIE \
  -fu_ice $FU \
  -out $LOG.repwvl.nc \
  | tee $LOG.repwvl
else
  echo "Skipping example computation because output $LOG.repwvl.nc exists"
fi

python $CMP $LOG.lbl.nc $LOG.repwvl.nc

#############
SOLAR="no"
THERMAL="yes"
LOG=${LOG_PREFIX}_SW${SOLAR}_LW${THERMAL}
if [ ! -e $LOG.lbl.nc ]; then
$BIN \
  -solver 2str \
  -solar   $SOLAR   \
  -thermal $THERMAL \
  -repwvl_data_thermal $ADA/repwvl/rdata.lbl.thermal.nc  \
  -repwvl_data_solar   $ADA/repwvl/rdata.lbl.solar.nc    \
  -mie_wc $MIE_GENERAL \
  -fu_ice $FU_GENERAL \
  -out $LOG.lbl.nc \
  | tee $LOG.lbl
else
  echo "Skipping example computation because output $LOG.lbl.nc exists"
fi

if [ ! -e $LOG.repwvl.nc ]; then
$BIN \
  -solver 2str \
  -solar   $SOLAR   \
  -thermal $THERMAL \
  -repwvl_data_thermal $REPWVL_LW \
  -repwvl_data_solar $REPWVL_SW \
  -mie_wc $MIE \
  -fu_ice $FU \
  -out $LOG.repwvl.nc \
  | tee $LOG.repwvl
else
  echo "Skipping example computation because output $LOG.repwvl.nc exists"
fi

python $CMP $LOG.lbl.nc $LOG.repwvl.nc
