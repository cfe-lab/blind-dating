#!/bin/sh


src_folder=${0%/3_analysis.sh}
tree_dir=$1
data_type=$2 # 0 = no dna, 1 = mixed (latency unknown), 2 = mixed (latency known)
use_rtt=$3

r -f $src_folder/plot.R --args $tree_dir $data_type $use_rtt