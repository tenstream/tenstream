PP="debug_gcc"
PP="prod_single_gcc"

#GRID=$1
#DATA=$2
#ATM=$3
TYPE="rrtmg"
TYPE="twostream"
TYPE="disort"
TYPE="plexrt"
TYPE="twostreamvsrayli"
TYPE="plexrtvsrayli"
TYPE="rayli"

TENSTREAM=$HOME/tenstream/
WDIR=$TENSTREAM/misc/plex_from_icon/

[ "x$GRID" == 'x' ] && GRID="$WDIR/grid.ifc_ham_55km-diam_0626m.nc"
[ "x$DATA" == 'x' ] && DATA="$WDIR/icon_input.nc"
#[ "x$DATA" == 'x' ] && DATA="$WDIR/1246_1m_DOM01_ML_20130505T130000Z.nc"
[ "x$ATM"  == 'x' ] && ATM="$WDIR/afglus.dat"

module purge
module load Modules petsc/$PP

BIN=$TENSTREAM/$PP/bin/ex_plex_ex4
#make -C $PETSC_DIR && \
make -j -C $TENSTREAM/$PP ex_plex_ex4 || exit

OUT=$HOME/scratch/plex_from_icon/out_${TYPE}
mkdir -p $(dirname $OUT)
cd $WDIR

#DEBUG="$DEBUG -start_in_debugger"
#DEBUG="$DEBUG -show_plex_coordinates"
#DEBUG="$DEBUG -show_migration_sf"

BASEOPT="-qv_data_string hus -lwc_data_string clw -qnc_data_string qnc -iwc_data_string cli -qni_data_string qni -thermal no"
#SOLVER="-N_first_bands_only 1"
SOLVER="$SOLVER -twostr_ratio 2"
NP=10

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
  SOLVER="$SOLVER -plexrt_use_rayli -rayli_photons 100000 -plexrt_vacuum_domain_boundary"
fi

if [ "x$TYPE" == "xplexrtvsrayli" ]; then
  SOLVER="$SOLVER -plexrt_vacuum_domain_boundary"
fi

rm -f $OUT.*

MPIOPT="mpirun -np $NP -wdir $WDIR"
SRUN="salloc -p ws -C GPU -n $NP --mem=30G bash -c "
$SRUN "$MPIOPT $BIN -grid $GRID -data $DATA -atm_filename $ATM -out $OUT.h5 $BASEOPT $SOLVER $DEBUG | tee $OUT.log"
[ -e $OUT.h5 ] && petsc_gen_xdmf.py $OUT.h5
