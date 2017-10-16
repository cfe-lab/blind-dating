#!/bin/bash

Rscript ../../src/regression.R --tree=trees.rooted/cartoon.nwk --info=info/cartoon.csv --patid=cartoon --real
Rscript ../../src/plot.R --tree=trees.rooted/cartoon.nwk --patid=cartoon  --distmin=-0.01 --distmax=.11 --distby=0.02 --yearstart=-1.01 --yearend=4 --yearby=1 --therapy=2.1 --cartoon
Rscript vl.R
