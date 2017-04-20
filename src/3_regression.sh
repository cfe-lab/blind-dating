#!/bin/bash

src_folder=${0%/3_regression.sh}

$TYPE=$1	# 0 Training/Censored, 1 RNA/DNA

mkdir plots stats 2>dev.null

cat patient_list.txt | while read p; do
	pat_id=$(echo $p | cut  -f ' ' -d 1)
	echo ${pat_id}
	R --slave --silent -f ${src_folder}/regression.R --args trees.rooted/${pat_id}.rtt.nwk patient_***REMOVED***.csv ${pat_id}
	R --slave --silent -f ${src_folder}/plot.R --args trees.rooted/${pat_id}.rtt.nwk patient_***REMOVED***.csv ${pat_id} $TYPE
done
