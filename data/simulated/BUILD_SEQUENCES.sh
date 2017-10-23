#!/bin/bash
#SBATCH --mem=1G

rep=$SLURM_ARRAY_TASK_ID

Rscript ../../src/gen_guide_trees.R ${rep} ${rep} ${rep} ref/HXB2.fasta
#python ../../src/fix_control.py control.unfixed.txt control.txt
#~/bin/indelible
#mv SIM_${rep}.fas simulated/
#mv SIM_${rep}_TRUE.fas simulated/
