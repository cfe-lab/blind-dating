#!/bin/bash

src=../../src
pat_id=patient_13334.cens

Rscript ${src}/plot.R --tree=trees.rooted/${pat_id}.rtt.nwk --patid=${pat_id} --mincoltime=0 --maxcoltime=748 --histfreqby=10
