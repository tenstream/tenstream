#!/bin/bash
. config.sh

ITER=${1:-}
IDENT=rayli${ITER}

OPT="\
  -solver rayli \
  -rayli_cyclic_bc \
  "

MPIEXEC="srun -J uclales-$IDENT -n1 -c32 --mem=180G --time=48:00:00"

pids=()
jobnr=0

for T in $(seq $TSTART 1 $TEND); do
  for SPEC in sw lw; do

    if [ x$SPEC == "xsw" ]; then
      S="-thermal no -solar yes"
    else
      S="-solar no -thermal yes"
    fi

    OUT=${OUTPREFIX}/$IDENT/$IDENT.$T.$SPEC.nc
    mkdir -p $(dirname $OUT)
    if [ ! -e $OUT ]; then
      sleep ${SLEEP:-0}
      ($MPIEXEC $TENSTREAMROOT/bin/$BIN $BASEOPT $OPT $S -tstart $T -tend $T -out $OUT &> >(tee $OUT.log)) &
      pids[${jobnr}]=$!
      jobnr=$(($jobnr +1))
    else
      echo "Skipping $OUT"
    fi
  done
done

trap 'for pid in ${pids[*]}; do echo "kill $pid"; kill ${pid}; done' EXIT
for pid in ${pids[*]}; do echo "waiting for job $pid to finish"; wait $pid; done
echo "Finished all jobs"

bash merge_rayli.sh $ITER
