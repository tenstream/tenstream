#!/bin/bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" # path to this file
PROJECT_ROOT="$(readlink -f $SCRIPTDIR/../)" # root of the tenstream repwvl dir

SCRIPTS=$PROJECT_ROOT/scripts # where repwvl scripts lie
DATA=$PROJECT_ROOT/data       # where tenstream relevant data will end up
WORK=$SCRATCH/repwvl_work     # where intermediate data is stored and can be deleted if project is finished

ADA=$WORK/ada                 # git clone of the ADA repo

NWVL_LBL_SW=2000 # number of line by line wavelengths (shortwave)
NWVL_LBL_LW=2000 # number of line by line wavelengths (longwave)

LBL_SW=$ADA/ARTSrun_solar/arts_line_by_line.nc
LBL_LW=$ADA/ARTSrun_thermal/arts_line_by_line.nc

TENSTREAM_ROOT=$TENSTREAMROOT # where the tenstream build is located

TRAIN_SW=$WORK/lbl_train_data_sw.nc
TRAIN_LW=$WORK/lbl_train_data_lw.nc

LIBRADTRAN_ROOT=$LIBRADTRANROOT # path to a libRadtran install (needed for Mie tool and fu ice data)

# General optprop files/settings:
MIE_GENERAL=$DATA/mie.wc.table.nc
FU_GENERAL=$DATA/fu.ice.general.nc


REPWVL_SAMPLES_SW=20
REPWVL_SAMPLES_LW=20

REPWVL_SW=$DATA/repwvl_sw_${REPWVL_SAMPLES_SW}.nc
REPWVL_LW=$DATA/repwvl_lw_${REPWVL_SAMPLES_LW}.nc

MIE_REPWVL=$DATA/mie.wc.repwvl_${REPWVL_SAMPLES_SW}_${REPWVL_SAMPLES_LW}.nc
FU_REPWVL=$DATA/fu.ice.repwvl_${REPWVL_SAMPLES_SW}_${REPWVL_SAMPLES_LW}.nc
