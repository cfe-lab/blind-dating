#!/bin/bash
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

INCR=50
END=1000

name=$SLURM_JOB_NAME
inst=$SLURM_ARRAY_TASK_ID

rep=$SLURM_ARRAY_TASK_ID
pat_id=SIM_${rep}-with_ref
pat_id_short=SIM_${rep}
src=~/blind-dating/src/
raxml=raxmlHPC-PTHREADS-AVX
cpus=$SLURM_CPUS_PER_TASK

source ${name} && \
new_inst=$(expr "$inst" + "$INCR") && \
echo "next: "$new_inst && \
if [[ $new_inst -le $END ]]; then
  echo sbatch --get-user-env --mail-type=FAIL --mail-user=bjones@cfenet.ubc.ca --partition=slow --array=$new_inst --job-name=$name wrapper.sh
  sbatch --get-user-env --mail-type=FAIL --mail-user=bjones@cfenet.ubc.ca --partition=slow --array=$new_inst --job-name=$name wrapper.sh
fi
