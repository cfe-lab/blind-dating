#!/bin/sh

# creates simulated reads for the current working directory

src_folder=${0%/1_simulate_sequences.sh}
latent=$1
latency_rate=$2
r_seed=$3
indelible_seed=$4

r -f $src_folder/gen_guide_trees.R --args $1 $2 $3 $4
python $src_folder/fix_control.py control.unfixed.txt control.txt
indelible
mv *.fas simulated/

