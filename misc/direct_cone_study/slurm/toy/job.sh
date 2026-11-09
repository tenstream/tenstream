#!/bin/bash
#SBATCH --job-name=cone_toy
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --array=0-119
#SBATCH --output=/project/meteo/work/Fabian.Jakub/claude_cone_study/logs/toy_%a.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
D=/project/meteo/work/Fabian.Jakub/claude_cone_study
mkdir -p $D/toy/results
python3 -I $D/toy/case.py $SLURM_ARRAY_TASK_ID $D/toy/results
