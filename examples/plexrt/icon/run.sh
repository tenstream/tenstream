#!/bin/bash
# Run the ICON example (ex_plex_rrtmg_icon) on the grid and data in this directory.
# Usage: BUILD=<tenstream build dir> [TYPE=...] [GRID=...] [DATA=...] [ATM=...] [OUT=...] [SRUN=...] ./run.sh
# TYPE is one of: rrtmg twostream disort plexrt twostreamvsrayli plexrtvsrayli rayli

[ "x$TYPE" == 'x' ] && TYPE="plexrt"

WDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TENSTREAM=$(cd "$WDIR/../../.." && pwd)
[ "x$BUILD" == 'x' ] && BUILD="$TENSTREAM/build"

[ "x$GRID" == 'x' ] && GRID="$WDIR/grid.ifc_ham_55km-diam_0626m.nc"
[ "x$DATA" == 'x' ] && DATA="$WDIR/icon_input.nc"
[ "x$ATM"  == 'x' ] && ATM="$WDIR/afglus.dat"

BIN=$BUILD/bin/ex_plex_rrtmg_icon
make -j -C $BUILD ex_plex_rrtmg_icon || exit

[ "x$OUT" == 'x' ] && OUT=$WDIR/out/out_${TYPE}
mkdir -p $(dirname $OUT)
cd $WDIR

#DEBUG="$DEBUG -start_in_debugger"
#DEBUG="$DEBUG -show_plex_coordinates"
#DEBUG="$DEBUG -show_migration_sf"

BASEOPT="-qv_data_string hus -lwc_data_string clw -qnc_data_string qnc -iwc_data_string cli -qni_data_string qni -thermal no"
#SOLVER="-N_first_bands_only 1"
SOLVER="$SOLVER -twostr_ratio 2"
NP=10
NC=1

if [ "x$TYPE" == "xdisort" ]; then
  SOLVER="$SOLVER -disort_only -disort_delta_scale"
fi
if [ "x$TYPE" == "xrrtmg" ]; then
  SOLVER="$SOLVER -rrtmg_only"
fi
if [ "x$TYPE" == "xtwostream" ]; then
  SOLVER="$SOLVER -twostr_only"
fi
if [ "x$TYPE" == "xtwostreamvsrayli" ]; then
  SOLVER="$SOLVER -twostr_only -plexrt_vacuum_domain_boundary"
fi
if [ "x$TYPE" == "xrayli" ]; then
  NP=1
  NC=10
  SOLVER="$SOLVER -plexrt_use_rayli -rayli_photons 100000 -plexrt_vacuum_domain_boundary"
fi

if [ "x$TYPE" == "xplexrtvsrayli" ]; then
  SOLVER="$SOLVER -plexrt_vacuum_domain_boundary"
fi

rm -f $OUT.*

MPIOPT="mpirun -n $NP -wdir $WDIR"
# on a cluster, e.g.: SRUN="salloc -N 1 -n $NP -c $NC --mem=30G --time=08:00:00"
$SRUN bash -c "$MPIOPT $BIN -grid $GRID -data $DATA -atm $ATM -out $OUT.h5 $BASEOPT $SOLVER $DEBUG | tee $OUT.log"
[ -e $OUT.h5 ] && petsc_gen_xdmf.py $OUT.h5
