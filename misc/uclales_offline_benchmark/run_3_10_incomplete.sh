#!/bin/bash
. config.sh

INC=${1:-1}
ITER=${2:-1}
IDENT=3_10_incomplete_${INC}_${ITER}

OPT="\
  -solver 3_10 \
  -solar_diff_explicit \
  -solar_diff_ksp_max_it $ITER \
  -solar_diff_ksp_complete_initial_run \
  -solar_diff_ksp_converged_reason \
  -solar_diff_ksp_monitor \
  -thermal_diff_explicit \
  -thermal_diff_ksp_max_it $ITER \
  -thermal_diff_ksp_complete_initial_run \
  -thermal_diff_ksp_converged_reason \
  -thermal_diff_ksp_monitor \
  -accept_incomplete_solve \
  -absorption_by_coeff_divergence \
  -repwvl_pprts_atm_view no \
  "

MPIEXEC="${MPIEXEC:-srun -J uclales-$IDENT -n8 -c1 --mem=8G --time=08:00:00}"

pids=()
jobnr=0
for SPEC in sw lw
do

  if [ x$SPEC == "xsw" ]; then
    S="-thermal no"
  else
    S="-solar no"
  fi

  OUT=${OUTPREFIX}/$IDENT/$IDENT.$SPEC.nc
  mkdir -p $(dirname $OUT)
  if [ ! -e $OUT ]; then
    sleep ${SLEEP:-0}
    ($MPIEXEC $TENSTREAMROOT/bin/$BIN $BASEOPT $OPT $S -tstart $TSTART -tend $TEND -tinc $INC -out $OUT &> >(tee $OUT.log)) &
    pids[${jobnr}]=$!
    jobnr=$(($jobnr +1))
  else
    echo "Skipping $OUT"
  fi
done

trap 'for pid in ${pids[*]}; do echo "kill $pid"; kill ${pid}; done' EXIT
for pid in ${pids[*]}; do echo "waiting for job $pid to finish"; wait $pid; done
echo "Finished all jobs"
