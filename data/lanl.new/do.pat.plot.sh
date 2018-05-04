#!/bin/bash

src=../../src
pat_id=patient_821

Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id} --mincoltime=0 --maxcoltime=2373 --histfreqby=10 --yearstart=-100 --histheight=1.1 --dnashapescale=0.6
