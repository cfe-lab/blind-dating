#!/bin/sh

# creates simulated reads for the current working directory

src_folder=${0%/1_simulate_sequences.sh}
unlatent_rate=$1
latent_rate=$2
r_seed=$3
indelible_seed=$4

r -f $src_folder/gen_guide_trees.R --args $unlatent_rate $latent_rate $r_seed $indelible_seed
python $src_folder/fix_control.py control.unfixed.txt control.txt
indelible
mv *.fas simulated/

