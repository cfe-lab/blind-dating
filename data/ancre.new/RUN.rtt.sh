#!/bin/bash
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4


raxml=raxmlHPC-PTHREADS-AVX
pat_id=$SLURM_JOB_NAME
src=~/blind-dating/src
rep=1989
cpus=$SLURM_CPUS_PER_TASK

if [ ! -e trees ]; then mkdir trees; fi
if [ ! -e trees.rooted ]; then mkdir trees.rooted; fi
if [ ! -e stats ]; then mkdir stats; fi

#echo "Building Tree"
#R --slave --silent -f ${src}/build.tree.R --args aligned/${pat_id}.fasta trees/${pat_id}.nwk ${raxml} $cpus ${rep}
echo "Rooting"
Rscript ${src}/root.tree.R trees/${pat_id}.nwk info/${pat_id}.csv trees.rooted/${pat_id}.rtt.nwk 1 0
echo "Regressing"
Rscript ${src}/regression.R --tree=trees.rooted/${pat_id}.rtt.nwk --info=info/${pat_id}.csv --patid=${pat_id}
echo "Plotting"
Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id}
