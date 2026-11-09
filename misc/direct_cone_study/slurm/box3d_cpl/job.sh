#!/bin/bash
#SBATCH --job-name=cone_cpl
#SBATCH --partition=cluster
#SBATCH --exclude=th-cl-262h05n3,usm-cl-826bac03,usm-cl-1116cs01
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-20
#SBATCH --output=/project/meteo/work/Fabian.Jakub/claude_cone_study/logs/cpl_%a.out
source ~/.profile >/dev/null 2>&1; luvspec >/dev/null 2>&1
cd /project/meteo/work/Fabian.Jakub/claude_cone_study/box3d_cpl
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" runs.txt)
name=${line%% *}; opts=${line#* }
rm -f ${name}_*
export UCX_TLS=self,tcp OMPI_MCA_btl=self,tcp
/usr/bin/time -f "$name %e s" mpirun -np 32 /project/meteo/work/Fabian.Jakub/claude_cone_study/bin/ex_pprts_box_cld_cpl $opts -pprts_xdmf ${name}_ < /dev/null > $name.log 2>&1
tail -1 $name.log
