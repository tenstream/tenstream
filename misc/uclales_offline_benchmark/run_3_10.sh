#!/bin/bash
. config.sh
IDENT=3_10

OPT="\
  -solver 3_10 \
  #-thermal_diff_ksp_atol 1e-60 \
  #-solar_diff_ksp_atol 1e-60 \
  #-solar_dir_ksp_atol 1e-60 \
  #-thermal_diff_ksp_converged_reason \
  #-solar_diff_ksp_converged_reason \
  #-solar_dir_ksp_converged_reason \
  "

MPIEXEC="${MPIEXEC:-srun -J uclales-$IDENT -n8 -c1 --mem=8G --time=08:00:00}"

pids=()
jobnr=0

#for T in $(seq $TSTART 1 $TEND); do
  for SPEC in sw lw; do

    if [ x$SPEC == "xsw" ]; then
      S="-thermal no"
    else
      S="-solar no"
    fi

    OUT=${OUTPREFIX}/$IDENT/$IDENT.$SPEC.nc
    mkdir -p $(dirname $OUT)
    if [ ! -e $OUT ]; then
      sleep ${SLEEP:-0}
      ($MPIEXEC $TENSTREAMROOT/bin/$BIN $BASEOPT $OPT $S -tstart $TSTART -tend $TEND -out $OUT &> >(tee $OUT.log)) &
      pids[${jobnr}]=$!
      jobnr=$(($jobnr +1))
    else
      echo "Skipping $OUT"
    fi
  done
#done

trap 'for pid in ${pids[*]}; do echo "kill $pid"; kill ${pid}; done' EXIT
for pid in ${pids[*]}; do echo "waiting for job ${pid} to finish"; wait ${pid}; done
echo "Finished all jobs"
