#!/bin/bash
#SBATCH --mem=1G

cat control.unfixed.head.txt > control.unfixed.txt

for i in $(seq 1 1000); do head -n 1 control.unfixed.${i}.txt; done >> control.unfixed.txt
for i in $(seq 1 1000); do head -n 2 control.unfixed.${i}.txt | tail -n 1; done >> control.unfixed.txt
echo '[EVOLVE]' >> control.unfixed.txt
for i in $(seq 1 1000); do head -n 3 control.unfixed.${i}.txt | tail -n 1; done >> control.unfixed.txt

python ../../src/fix_control.py control.unfixed.txt control.txt
~/bin/indelible
mv *.fas simulated/
