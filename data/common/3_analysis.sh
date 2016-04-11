#!/bin/sh


src_folder=${0%/3_analysis.sh}
tree_dir=$1
data_type=$2 # 0 = no dna, 1 = mixed (latency unknown), 2 = mixed (latency known)
use_rtt=$3
r_seed=$4
suffix=$5
data_name=$6

r -f $src_folder/plot.R --args "$tree_dir" $data_type $use_rtt $r_seed
mv data.csv data${suffix}.csv
mv stats.csv stats${suffix}.csv
mv hist.pdf hist${suffix}.pdf
mv plot.pdf plot${suffix}.pdf

r -f $src_folder/do_analysis.R --args "${suffix}" "$data_name"