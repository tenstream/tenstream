#!/bin/bash

set -euo pipefail

. config.sh

if [ ! -e $TENSTREAM_ROOT/bin/compute_repwvl_training_data ]; then
  echo "Building compute_repwvl_training_data in TenStream build dir"
  make -j -C $TENSTREAM_ROOT compute_repwvl_training_data
fi

MPIEXEC=${MPIEXEC:-srun -n32 ${MODULES_MARCH:+-C} $MODULES_MARCH -p met-cl,met-ws,cluster --mem=40G --time=01:00:00}

if [ ! -e $TRAIN_LW ]; then
  ${MPIEXEC} \
    $TENSTREAM_ROOT/bin/compute_repwvl_training_data \
    -thermal \
    -atm $ADA/Atmosphere_Data/garand_profiles.nc \
    -bg_atm $PROJECT_ROOT/../examples/pprts/afglus_100m.dat \
    -repwvl_data_thermal $ADA/repwvl/rdata.lbl.thermal.nc \
    -repwvl_data_solar   $ADA/repwvl/rdata.lbl.solar.nc \
    -mie_wc $MIE_GENERAL \
    -fu_ice $FU_GENERAL \
    -out $TRAIN_LW
else
  echo "Skipping compute_repwvl_training_data thermal because <$TRAIN_LW> already exists"
fi

if [ ! -e $TRAIN_SW ]; then
  ${MPIEXEC} \
    $TENSTREAM_ROOT/bin/compute_repwvl_training_data \
    -solar \
    -atm $ADA/Atmosphere_Data/garand_profiles.nc \
    -bg_atm $PROJECT_ROOT/../examples/pprts/afglus_100m.dat \
    -repwvl_data_thermal $ADA/repwvl/rdata.lbl.thermal.nc \
    -repwvl_data_solar $ADA/repwvl/rdata.lbl.solar.nc \
    -mie_wc $MIE_GENERAL \
    -fu_ice $FU_GENERAL \
    -out $TRAIN_SW
else
  echo "Skipping compute_repwvl_training_data solar because <$TRAIN_SW> already exists"
fi
