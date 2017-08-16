#!/bin/bash

TENSTREAM_DIR=$HOME/tenstream/

BIN=$TENSTREAM_DIR/build/bin/ex_rrtm_lw_sw

NPROC=$(cat /proc/cpuinfo |grep processor|wc -l)

SCALING_LOG="scaling.log"
SCALING_OUTPUT="scaling.result"

REPETITIONS=10
MPI_OPT="--bind-to core --map-by socket"

echo "Running tenstream rrtmg benchmark with $NPROC processors"
rm -f $SCALING_OUTPUT
for (( c=1; c<=$NPROC; c++ ))
do
  for (( r=1; r<=$REPETITIONS; r++ ))
  do
    echo "Running benchmark for NP==$c with mpi options $MPI_OPT"
    mpirun -np $c $MPI_OPT --wdir $TENSTREAM_DIR/examples/rrtm_lw_sw/ $BIN -log_summary | tee $SCALING_LOG
    TOTAL=$(grep "Time (sec):" $SCALING_LOG| awk '{print $3}')
    DIFFSOLVE=$(grep "calc_ediff:" $SCALING_LOG| awk '{print $3}')
    
    echo $c $TOTAL $DIFFSOLVE >> $SCALING_OUTPUT
  done
done

cat > scaling.py <<EOF
from pylab import *
res = np.loadtxt("$SCALING_OUTPUT")

# Filter for repetitions, i.e. search for minimum solve times
data = []
for ncpu in np.unique(res[:,0]):
  min_total = np.min(res[res[:,0]==ncpu][:,1])
  min_diff_solve = np.min(res[res[:,0]==ncpu][:,2])
  data.append([ncpu, min_total, min_diff_solve])

data = np.array(data)

plt.plot(data[:,0], data[:,1], label='Total Time')
plt.plot(data[:,0], data[:,2], label='Diffuse Solver')
xlabel('Number of Cores')
ylabel('Time to solution')
title('weak scaling, approx. {} dof per process'.format(3*3*10*164))
legend()
plt.savefig('scaling.png')
EOF

python scaling.py
