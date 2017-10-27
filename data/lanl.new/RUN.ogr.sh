#!/bin/bash

pat_id=$1
src=~/git/blind-dating/src

echo "Rooting"
Rscript ${src}/root.tree.R trees/${pat_id}.ogr.nwk info/${pat_id}.csv trees.rooted/${pat_id}.ogr.nwk 0 0
echo "Regressing"
Rscript ${src}/regression.R --tree=trees.rooted/${pat_id}.ogr.nwk --info=info/${pat_id}.csv --patid=${pat_id}.ogr --real
#echo "Plotting"
#Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.ogr.nwk --patid=${pat_id} --real --distmax=0.27 --distby=.04 --yearstart=1993 --yearend=2018 --therapy="2006-08-01"
