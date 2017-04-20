#!/bin/sh

# creates simulated reads for the current working directory

src_folder=${0%/1_simulate_sequences.sh}
percent_censored=$1
r_seed=$2
indelible_seed=$3

R --silent --slave -f $src_folder/gen_guide_trees.R --args ${percent_censored} ${r_seed} ${indelible_seed}
python $src_folder/fix_control.py control.unfixed.txt control.txt
indelible
mkdir simulated 2>/dev/null
mv *.fas simulated/

