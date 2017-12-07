#!/bin/bash
#SBATCH --mem=1G

reps=$1

Rscript ../../src/gen_guide_trees.R head 1 1989 ref/HXB2.fasta

cat control/control.unfixed.head.txt > control/control.unfixed.txt

for i in $(seq 1 ${reps}); do head -n 1 control/control.unfixed.${i}.txt; done >> control/control.unfixed.txt
for i in $(seq 1 ${reps}); do head -n 2 control/control.unfixed.${i}.txt | tail -n 1; done >> control/control.unfixed.txt
echo '[EVOLVE]' >> control/control.unfixed.txt
for i in $(seq 1 ${reps}); do head -n 3 control/control.unfixed.${i}.txt | tail -n 1; done >> control/control.unfixed.txt

python ../../src/fix_control.py control/control.unfixed.txt control.txt
~/bin/indelible

Rscript ../../src/censor.simulated.R info/ 50 1989
