#TESTIDENT="thermal_diff"
#LSOLAR=False
#LTHERMAL=True
#DX=100

TESTIDENT=$1
LSOLAR=$2
LTHERMAL=$3
DX=$4

OUT=out_${TESTIDENT}_${DX}
mkdir -p $OUT

RESULT="${TESTIDENT}_log.npy"
RUNLOG="$OUT/${TESTIDENT}_run.log"
rm -f $OUT/$RUNLOG #$OUT/$RESULT

cat > $OUT/concat_${TESTIDENT}_logs.py << EOF
import os, sys, glob, numpy as np
logfiles = glob.glob('solver_log_${TESTIDENT}_*')
logfiles.sort(key=os.path.getmtime)
result_file = "$RESULT"
print("Result file:", result_file)

result = list(np.load(result_file, allow_pickle=True)) if os.path.exists(result_file) else list();
for ilogfile, logfile in enumerate(logfiles):
    have_entry_already=False
    for r in result:
        if r['logfile'] == logfile:
            have_entry_already=True
            break
    if have_entry_already:
        print("Have entry already", logfile)
        continue
    print("Appending: {} ({}/{})".format(logfile, ilogfile, len(logfiles)))
    try:
        exec(open(logfile).read());
    except NameError as e:
        continue

    delkeys = []
    for k in Stages['plexrtsolve_stage'].keys():
        try:
            count = Stages['plexrtsolve_stage'][k][1].get('count')
            if count == 0:
                delkeys.append(k)
        except Exception as e:
            pass
        try:
            count = Stages['plexrtsolve_stage'][k][1][0].get('count')
            if count == 0:
                delkeys.append(k)
        except Exception as e:
            pass
    [Stages['plexrtsolve_stage'].pop(k) for k in delkeys ]
    Stages.pop('Main Stage')

    label = "\n         -".join(Stages['description'].split(' -'))
    print(label)
    Stages['logfile'] = logfile
    result.append(Stages);
    del Stages
    if ilogfile%10==0:
        np.save(result_file, result)
np.save(result_file, result)
EOF

cat > $OUT/plot_${TESTIDENT}.py << EOF
from pylab import *
results=list(np.load("$RESULT", allow_pickle=True))
ldir=$LSOLAR

stages = []
if ldir:
  stages.append('plexrtsolve_Mdir'  )
  stages.append('plexrtsetup_Mdir'  )
  stages.append('plexrtdir2dir'     )
  stages.append('plexrtdir2diff'    )

stages.append('plexrtsolve_Mdiff' )
stages.append('plexrtsetup_Mdiff' )
stages.append('plexrtdiff2diff'   )
stages.append('PCSetUp'           )
stages.append('PCApply'           )
stages.append('KSPSetUp'          )
stages.append('KSPSolve'          )
stages.append('summary'           )

for isolve, solve in enumerate(results):
  print("Description: {}".format(solve['description']))
  for stage in stages:
    Nranks = len(solve['plexrtsolve_stage'][stage])
    time = np.sum([ solve['plexrtsolve_stage'][stage][id]['time'] for id in range(Nranks) ])
    print("{:>30} : {:5.4f}".format(stage, time))
  print()

from mpldatacursor import datacursor
lines=[]

if ldir:
  figure(1)
  clf(); title('times for dir solves ${TESTIDENT}')
  stages = ['plexrtsolve_Mdir', 'plexrtsetup_Mdir']

  for isolve, solve in enumerate(results):
      label = "dir" + solve['logfile'] + "\n-".join(solve['description'].split(' -'))
      Nranks = len(solve['plexrtsolve_stage'][stages[0]])
      time = np.sum([ [ solve['plexrtsolve_stage'][stage][id]['time'] for stage in stages ] for id in range(Nranks) ])
      line, = plot(isolve, time, 'o', label=label)
      lines.append(line)

  datacursor(display='multiple', draggable=True, fontsize='x-small',bbox=dict(alpha=.9))
  grid(True, axis='y', which='both')

figure(2)
clf(); title('times for diff solves ${TESTIDENT}')
stages = ['plexrtsolve_Mdiff', 'plexrtsetup_Mdiff']

for isolve, solve in enumerate(results):
    label = "diff" + solve['logfile'] + "\n-".join(solve['description'].split(' -'))
    Nranks = len(solve['plexrtsolve_stage'][stages[0]])
    time = np.sum([ [ solve['plexrtsolve_stage'][stage][id]['time'] for stage in stages ] for id in range(Nranks) ])
    line, = plot(isolve, time, 'o', label=label)
    lines.append(line)

datacursor(display='multiple', draggable=True, fontsize='x-small',bbox=dict(alpha=.9))
grid(True, axis='y', which='both')

plt.show()
EOF


WDIR="$HOME/tenstream/build/bin/"
LOGGING="-solar_dir_ksp_converged_reason -solar_diff_ksp_converged_reason -thermal_diff_ksp_converged_reason -${TESTIDENT}_ksp_view"
BASEOPT="-out /tmp/o.h5 -atm_filename atm.dat -dx $DX -Nx 32 -Ny 39 -Nz 40 -N_first_bands_only 10 -solar $LSOLAR -thermal $LTHERMAL -solve_iterations 10 -solve_iterations_scale 1 "
BASEOPT="$BASEOPT -solar_dir_pc_type asm -solar_dir_pc_asm_overlap 0 -solar_dir_sub_ksp_type richardson -solar_dir_sub_ksp_max_it 0 -solar_diff_sub_pc_type sor"
BASEOPT="$BASEOPT -solar_diff_pc_type asm -solar_diff_pc_asm_overlap 0 -solar_diff_sub_ksp_type richardson -solar_diff_sub_ksp_max_it 0 -solar_diff_sub_pc_type sor"
BIN="mpirun --bind-to core bin/ex_plex_rrtmg_fish"
BIN="mpirun $WDIR/ex_plex_rrtmg_fish"

runsolve() {
  SOLVER=$1
  HASH=$(echo $SOLVER|md5sum|cut -f1 -d" ")
  LOGFILE="$OUT/solver_log_${TESTIDENT}_${HASH}.py"
  [[ $(wc -l <file.txt) -le 10 ]] && rm $LOGFILE
  if [ ! -e $LOGFILE ]
  then
    echo "Computing $SOLVER -> $LOGFILE" |tee -a $RUNLOG
    LOG="-log_view :${LOGFILE}:ascii_info_detail"
    CMD="$BIN $BASEOPT $LOGGING $LOG $SOLVER"
    echo "Running: $CMD"
    #salloc -p vis -C GPU --ntasks-per-core=1 -N 1 -n 10 --time=00:05:00 -w met-cl-vis01 bash -c "\
    #salloc -p ws -C GPU -N 2 -n 20 --exclusive --time=00:05:00 -w met-ws-970r18,met-ws-970r17 bash -c "\
    salloc -p vis -N 1 -n 64 --exclusive --time=00:02:00 bash -c "\
      $CMD >> $RUNLOG; \
      echo \"Stages['description'] = \\\"$SOLVER\\\"\" |tee -a $LOGFILE; \
    "
  else
    echo "Skipping $SOLVER -> $LOGFILE" |tee -a $RUNLOG
  fi
  #python3 concat_logs.py $LOGFILE "$SOLVER"
}

#ILU Run
for it in {0..3}
do
  SOLVER="\
    -${TESTIDENT}_pc_type bjacobi\
    -${TESTIDENT}_sub_pc_type ilu\
    -${TESTIDENT}_sub_pc_factor_levels $it\
    -${TESTIDENT}_sub_pc_factor_fill $((it+2))\
    "
  runsolve "$SOLVER"
done

#ASM Eisenstat Runs
for overlap in {0..2}
do
  for it in {0..2}
  do
    for richardson in "" "-${TESTIDENT}_ksp_richardson_self_scale"
    do
      SOLVER="\
        -${TESTIDENT}_pc_type asm\
        -${TESTIDENT}_pc_asm_overlap $overlap \
        -${TESTIDENT}_sub_ksp_type richardson \
        $richardson \
        -${TESTIDENT}_sub_ksp_max_it $it \
        -${TESTIDENT}_sub_pc_type eisenstat\
        "
      runsolve "$SOLVER"
    done
  done
done
#ASM SOR Runs
for overlap in {0..2}
do
  for it in {0..2}
  do
    for richardson in "" "-${TESTIDENT}_ksp_richardson_self_scale"
    do
      SOLVER="\
        -${TESTIDENT}_pc_type asm\
        -${TESTIDENT}_pc_asm_overlap $overlap \
        -${TESTIDENT}_sub_ksp_type richardson \
        -${TESTIDENT}_sub_ksp_max_it $it \
        -${TESTIDENT}_sub_pc_type sor\
        $richardson \
        "
      runsolve "$SOLVER"
    done
  done
done
for overlap in {0..2}
do
  for it in {0..2}
  do
      SOLVER="\
        -${TESTIDENT}_pc_type asm\
        -${TESTIDENT}_pc_asm_overlap $overlap \
        -${TESTIDENT}_sub_ksp_type preonly \
        -${TESTIDENT}_sub_pc_type sor\
        "
      runsolve "$SOLVER"
  done
done

#ASM ILU Runs
for overlap in {0..2}
do
  for it in {0..2}
  do
    for richardson in "" "-${TESTIDENT}_ksp_richardson_self_scale"
    do
      SOLVER="\
        -${TESTIDENT}_pc_type asm\
        -${TESTIDENT}_pc_asm_overlap $overlap \
        -${TESTIDENT}_sub_ksp_type richardson \
        -${TESTIDENT}_sub_ksp_max_it $it \
        -${TESTIDENT}_sub_pc_type ilu\
        $richardson \
        "
      runsolve "$SOLVER"
    done
  done
done

# GAMG Runs
for limit in 100000 100
do
  for square in 0 1 5 10
  do
    for it in {0..3}
    do
      for thresh in {2..-3..-1}
      do
        for levels_ksp in "preonly" "richardson" "cg" "bcgs"
        do
          for levels_pc in "bjacobi -${TESTIDENT}_mg_levels_sub_pc_type ilu" "sor"
          do
            for coarse_ksp in "gmres" "preonly" "richardson" "cg"
            do
              for coarse_pc in "sor" "bjacobi -${TESTIDENT}_mg_coarse_sub_pc_type ilu"
              do
                for coarse_it in 0 1 10 100
                do
                  for GAMG_OPT in "" "-${TESTIDENT}_pc_gamg_reuse_interpolation"
                  do
                    SOLVER="\
                      -${TESTIDENT}_pc_type gamg\
                      -${TESTIDENT}_pc_gamg_agg_nsmooths 0\
                      -${TESTIDENT}_pc_gamg_sym_graph true\
                      -${TESTIDENT}_pc_gamg_threshold 1e$thresh \
                      -${TESTIDENT}_pc_gamg_square_graph $square \
                      -${TESTIDENT}_pc_gamg_coarse_eq_limit $limit \
                      -${TESTIDENT}_mg_levels_ksp_type richardson\
                      -${TESTIDENT}_mg_levels_ksp_richardson_self_scale\
                      -${TESTIDENT}_mg_levels_ksp_max_it $it \
                      -${TESTIDENT}_mg_levels_pc_type $levels_pc \
                      -${TESTIDENT}_mg_coarse_ksp_type $coarse_ksp \
                      -${TESTIDENT}_mg_coarse_ksp_max_it $coarse_it \
                      -${TESTIDENT}_mg_coarse_pc_type $coarse_pc \
                      $GAMG_OPT \
                      "
                    runsolve "$SOLVER"
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

#python3 concat_${TESTIDENT}_logs.py && python3 plot_${TESTIDENT}.py
