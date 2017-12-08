#!/bin/bash
#SBATCH --cpus=4
#SBATCH --mem=2G

rep=$SLURM_ARRAY_TASK_ID
pat_id=SIM_${rep}-with_ref
pat_id_short=SIM_${rep}
src=~/blind-dating/src/
raxml=raxmlHPC-PTHREADS-AVX
cpus=$SLURM_CPUS_PER_TASK

echo "Building" && \
R --slave --silent -f ${src}/build.tree.R --args simulated/${pat_id_short}_TRUE-with_ref.fasta trees/${pat_id}.nwk ${raxml} ${cpus} ${rep} && \
echo "Rooting" && \
R --slave --silent -f ${src}/root.tree.R --args trees/${pat_id}.nwk info/${pat_id_short}.cens.csv trees.rooted/${pat_id}.ogr.nwk 1 0 && \
echo "Regressing" && \
R --slave --silent -f ${src}/regression.R --args --tree=trees.rooted/${pat_id}.ogr.nwk --info=info/${pat_id_short}.cens.csv --patid=${pat_id}
