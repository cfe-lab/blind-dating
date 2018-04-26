#!/bin/bash

src=../../src
pat_id=SIM_1

Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id} --mincoltime=0 --maxcoltime=2304.71091998375 --histfreqby=5
