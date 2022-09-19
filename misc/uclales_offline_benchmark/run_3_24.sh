#!/bin/bash
. config.sh
IDENT=3_24

OPT="\
  -solver 3_24 \
  "

MPIEXEC="${MPIEXEC:-srun -J uclales-$IDENT -n8 -c1 --mem=32G --time=08:00:00}"

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
