#!/bin/bash

src=../../src
pat_id=patient_821

Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id} --mincoltime=0 --maxcoltime=2373 --histfreqby=10 --histheight=1.1 --dnashapescale=0.6 --yearstart=-100 --marklatent

pat_id=patient_13889

Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id} --mincoltime=0 --maxcoltime=4837 --histfreqby=10 --yearby=730.5 --histheight=1.1 --dnashapescale=0.6 --marklatent
