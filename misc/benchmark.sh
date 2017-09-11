#!/bin/bash

TENSTREAM_DIR=$HOME/tenstream/

BIN=$TENSTREAM_DIR/build/bin/ex_rrtm_lw_sw

NPROC=$(cat /proc/cpuinfo |grep processor|wc -l)

SCALING_LOG="scaling.log"
SCALING_OUTPUT="scaling.result"
RUNS_OUTPUT="scaling.full_log"

REPETITIONS=3
MPI_OPT="--bind-to core --map-by numa"
EX_OPT="-log_view -dir_ksp_atol 1e-6 -diff_ksp_atol 1e-6 -dir_ksp_rtol 1e-5 -diff_ksp_rtol 1e-5 -dir_ksp_converged_reason -diff_ksp_converged_reason -Nx 3 -Ny 3 -Nz 50"
SOLVER_OPT="-dir_ksp_type fbcgs -dir_pc_type bjacobi -dir_sub_pc_type ilu -dir_sub_pc_factor_levels 1 -dir_sub_ksp_type richardson -diff_ksp_type fbcgs -diff_pc_type bjacobi -diff_sub_pc_type ilu -diff_sub_pc_factor_levels 1 -diff_sub_ksp_type richardson"
#SOLVER_OPT="-twostr_only"

echo "Running tenstream rrtmg benchmark with $NPROC processors"
rm -f $SCALING_OUTPUT $RUNS_OUTPUT
for (( c=1; c<=$NPROC; c++ ))
do
  for (( r=1; r<=$REPETITIONS; r++ ))
  do
    echo "Running benchmark for NP==$c with mpi options $MPI_OPT" | tee -a $RUNS_OUTPUT
    echo "mpirun -np $c $MPI_OPT --wdir $TENSTREAM_DIR/examples/rrtm_lw_sw/ $BIN $EX_OPT" | tee -a $RUNS_OUTPUT
    mpirun -np $c $MPI_OPT --wdir $TENSTREAM_DIR/examples/rrtm_lw_sw/ $BIN $EX_OPT $SOLVER_OPT | tee $SCALING_LOG | tee -a $RUNS_OUTPUT
    TOTAL=$(grep "Time (sec):" $SCALING_LOG| awk '{print $3}')
    DIFFSOLVE=$(grep "calc_ediff:" $SCALING_LOG| awk '{print $3}')
    DIRSOLVE=$(grep "calc_edir:" $SCALING_LOG| awk '{print $3}')

    echo $c $TOTAL $DIFFSOLVE $DIRSOLVE >> $SCALING_OUTPUT
  done
done

cat > scaling.py <<EOF
from pylab import *
res = np.loadtxt("$SCALING_OUTPUT")

# Filter for repetitions, i.e. search for minimum solve times
data = []
for ncpu in np.unique(res[:,0]):
  min_total = np.min(res[res[:,0]==ncpu][:,1:3])
  min_diff_solve = np.min(res[res[:,0]==ncpu][:,3:5])
  min_dir_solve = np.min(res[res[:,0]==ncpu][:,5:7])
  data.append([ncpu, min_total, min_diff_solve, min_dir_solve])

data = np.array(data)

plt.plot(data[:,0], data[:,1], label='Total Time')
plt.plot(data[:,0], data[:,2], label='Diffuse Solver')
plt.plot(data[:,0], data[:,3], label='Direct Solver')
plt.plot(data[:,0], data[:,1] - data[:,2] - data[:,3], label='Miscellaneous')
xlabel('Number of Cores')
ylabel('Time to solution [s]')
title('weak scaling, approx. {} dof per process'.format(6*6*85*100))
legend(loc='best')
plt.savefig('scaling.png')
EOF

python scaling.py
