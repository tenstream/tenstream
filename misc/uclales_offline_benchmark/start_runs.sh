#!/bin/bash
. config.sh

make -j -C $TENSTREAMROOT $BIN

SLEEP=3

bash run_3_10.sh &
bash run_3_10_incomplete.sh 1 2 &
bash run_3_10_incomplete.sh 1 4 &
bash run_3_10_incomplete.sh 1 8 &
bash run_3_10_incomplete.sh 6 2 &
bash run_3_10_incomplete.sh 6 4 &
bash run_3_10_incomplete.sh 6 8 &
bash run_3_10_incomplete.sh 6 12 &
bash run_2str.sh &
bash run_2str_rrtmg.sh &
bash run_rayli.sh 0 &
#bash run_rayli.sh 1 &

python plot_scaling.py scaling.pdf $WORK/full/{2,3}*/*.log
