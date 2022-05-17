#!/bin/bash

set -euo pipefail

. config.sh

ADA_URL=git@gitlab.lrz.de:libradtran/ada.git
ADA_BRANCH=FJ/optipy

if [ ! -e $ADA ]; then
  git clone -b $ADA_BRANCH $ADA_URL $ADA
fi

sed -i "s/NUMWVL=.*/NUMWVL=$NWVL_LBL_SW/" $ADA/config.solar
sed -i "s/NUMWVL=.*/NUMWVL=$NWVL_LBL_LW/" $ADA/config.thermal

cd $ADA
./01_prep_ARTS.sh
./02_get_HITRAN.sh

for cfg in config.solar config.thermal
do
  bash 03_run_ARTS.sh $cfg
  bash 07_run_repwvl_test.sh $cfg
done

if [ ! -e $LBL_SW ]; then
  echo "Expected to find ARTS line by line solar output at <$LBL_SW> but its not there"
  exit 1
fi
if [ ! -e $LBL_LW ]; then
  echo "Expected to find ARTS line by line thermal output at <$LBL_LW> but its not there"
  exit 1
fi
