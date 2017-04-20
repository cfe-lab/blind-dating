#!/bin/sh

# Just a dummy wrapper around the 
# build tree script to keep everything organised

src_folder=${0%/2_build_trees.sh}
seed=$1

mkdir trees trees.rooted 2> /dev/null

cat patient_list.txt | while read p; do
	pat_id=$(echo $p | cut  -f ' ' -d 1)
	fasta=$(echo $p | cut -f ' ' -d 2-)
	echo ${pat_id}
	R --slave --silent -f ${src_folder}/build.tree.R --args $fasta trees/${pat_id}.nwk 10 ${seed}
	R --slave --silent -f ${src_folder}/root.tree.R --args trees/${pat_id}.nwk info/patient_***REMOVED***.csv trees.rooted/${pat_id}.rtt.nwk 1
done
