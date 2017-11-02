#!/bin/bash
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

raxml=raxmlHPC-PTHREADS-AVX
pat_id=$SLURM_JOB_NAME
full_pat_id=${pat_id}-rms
src=~/blind-dating/src
rep=1989
cpus=$SLURM_CPUS_PER_TASK

if [ ! -e trees ]; then mkdir trees; fi
if [ ! -e trees.rooted ]; then mkdir trees.rooted; fi
if [ ! -e stats ]; then mkdir stats; fi

#echo "Building Tree"
#R --slave --silent -f ${src}/build.tree.R --args aligned/${pat_id}.fasta trees/${pat_id}.nwk ${raxml} $cpus ${rep}
echo "Rooting"
Rscript ${src}/root.tree.R trees/${pat_id}.nwk info/${pat_id}.csv trees.rooted/${pat_id_full}.rtt.nwk 1 0
echo "Regressing"
Rscript ${src}/regression.R --tree=trees.rooted/${pat_id_full}.rtt.nwk --info=info/${pat_id}.csv --patid=${pat_id_full}
#echo "Plotting"
#Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.ogr.nwk --patid=${pat_id} --real --distmax=0.27 --distby=.04 --yearstart=1993 --yearend=2018 --therapy="2006-08-01"
