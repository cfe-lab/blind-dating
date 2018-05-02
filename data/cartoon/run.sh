#!/bin/bash

Rscript ../../src/regression.R --tree=trees.rooted/cartoon.nwk --info=info/cartoon.csv --patid=cartoon
Rscript ../../src/plot.R --tree=trees.rooted/cartoon.nwk --patid=cartoon  --distmin=-0.01 --distmax=.11 --distby=0.02 --yearstart=-2.01 --yearend=8 --yearby=1 --therapy=4.4666666 --cartoon --histfreqby=1 --maxcoltime=4.166666666 --mincoltime=0 --histheight=1.1
Rscript vl.R
