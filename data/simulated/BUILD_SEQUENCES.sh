#!/bin/bash
#SBATCH --mem=1G

rep=$SLURM_ARRAY_TASK_ID

Rscript ../../src/gen_guide_trees.R ${rep} ${rep} 1989 ref/HXB2.fasta
