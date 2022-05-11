#!/bin/bash

set -euo pipefail

LIBRADTRAN=$1

python gen_fu_ice_tables.py $LIBRADTRAN -v fu.ice.general.nc
echo ""
python gen_fu_ice_tables.py $LIBRADTRAN -v --repwvl_file repwvl_thermal.lut fu.ice.repwvl.nc
echo ""
python gen_fu_ice_tables.py $LIBRADTRAN -v --repwvl_file repwvl_solar.lut fu.ice.repwvl.nc
