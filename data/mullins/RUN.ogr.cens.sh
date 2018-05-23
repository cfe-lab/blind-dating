#!/bin/bash
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1


raxml=raxmlHPC-PTHREADS-AVX
pat_id_tree=${SLURM_JOB_NAME}-with_ref
pat_id_short=${SLURM_JOB_NAME}.cens
pat_id=${pat_id_short}-with_ref
src=~/blind-dating/src
rep=1989
cpus=$SLURM_CPUS_PER_TASK

echo $pat_id_tree
echo $pat_id_short
echo $pat_id
echo trees.rooted/${pat_id}.ogr.nwk

if [ ! -e trees ]; then mkdir trees; fi
if [ ! -e trees.rooted ]; then mkdir trees.rooted; fi
if [ ! -e stats ]; then mkdir stats; fi
if [ ! -e plots ]; then mkdir plots; fi

#echo "Building Tree"
#R --slave --silent -f ${src}/build.tree.R --args aligned/${pat_id}.fasta trees/${pat_id}.nwk ${raxml} $cpus ${rep}
#echo "Rooting"
#Rscript ${src}/root.tree.R trees/${pat_id_tree}.nwk info/${pat_id_short}.csv trees.rooted/${pat_id}.ogr.nwk 0 0
echo "Regressing"
Rscript ${src}/regression.R --tree=trees.rooted/${pat_id}.ogr.nwk --info=info/${pat_id_short}.csv --patid=${pat_id}
#echo "Plotting"
#Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.ogr.nwk --patid=${pat_id} 
