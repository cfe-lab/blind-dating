#!/bin/bash
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

if [ -z "$SLURM_WRAPPER_INCR" ]; then
SLURM_WRAPPER_INCR=1
fi

if [ -z "$SLURM_WRAPPER_END" ]; then
SLURM_WRAPPER_END=10
fi

name=$SLURM_JOB_NAME
inst=$SLURM_ARRAY_TASK_ID

rep=$SLURM_ARRAY_TASK_ID
pat_id=SIM_${rep}-with_ref
pat_id_short=SIM_${rep}
src=~/blind-dating/src/
raxml=raxmlHPC-PTHREADS-AVX
cpus=$SLURM_CPUS_PER_TASK

source ${name} && \
new_inst=$(expr "$inst" + "$SLURM_WRAPPER_INCR") && \
echo "next: "$new_inst && \
if [[ $new_inst -le $SLURM_WRAPPER_END ]]; then
  echo /opt/scyld/slurm/bin/sbatch --get-user-env --export="SLURM_WRAPPER_INCR=$SLURM_WRAPPER_INCR,SLURM_WRAPPER_END=$SLURM_WRAPPER_END" --mail-type=FAIL --mail-user=bjones@cfenet.ubc.ca --partition=slow --array=$new_inst --job-name=$name wrapper.sh
  /opt/scyld/slurm/bin/sbatch --get-user-env --export="SLURM_WRAPPER_INCR=$SLURM_WRAPPER_INCR,SLURM_WRAPPER_END=$SLURM_WRAPPER_END" --mail-type=FAIL --mail-user=bjones@cfenet.ubc.ca --partition=slow --array=$new_inst --job-name=$name wrapper.sh
fi
