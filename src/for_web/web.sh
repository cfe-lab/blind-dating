#!/bin/bash

# Proviral dating pipeline for web server
#
# web.sh INPUT_FOLDER OUTPUT_FOLDER
# 
# INPUT_FOLDER: path to the input folder
# 
# OUTPUT_FOLDER: path to the output folder
#
# INPUT_FOLDER must contain:
# (1) tree.nwk, newick tree
# (2) info.csv, comma delimited file with the columns: FULLSEQID (tips names of the newick tree),
#   COLDATE (collection date of sequences), CENSORED (0 if training and 1 if censored)
# (3) runid.txt, file containing a run identifier
#
# OUTPUT_FOLDER will contain: 
#
# (1) rooted_tree.nwk, newick tree rooted by RTT
# (2) stats.csv, one row (plus header) comma delimted file with columns: RunID (run identifier),
#   dAIC (difference between null AIC and linear model AIC), EstimatedRootDate (estimated date
#   of the root as per the linear model), EstimatedRootDate95(Low|High) (95% confidence interval of 
#   estimated root date), EstimatedEvolutionaryRate (estimated evolutionary rate as per the linear
#   model)
# (3) data.csv, comma delimited file with columns: ID (sequence ID), EstimatedDate (estimated 
#   date of the sequence as per the linear model), EstimatedDate(Low|High) (95% confidence 
#   interval of date of the sequence)
# (4) divergence_versus_time.pdf, PDF file of divergence versus time plot of the sequences with 
#   linear regression and ancestral traces
# (5) divergence_versus_time.png, PNG version of (4)


# R source directory
### TODO: change to approriate folder ###
BDSRC=~/working/blind-dating/src/for_web

# Arguments
INPUT=$1
OUTPUT=$2

# Input/Output paths  
TREE=${INPUT}/tree.nwk
INFO=${INPUT}/info.csv
RUNID=$(cat ${INPUT}/runid.txt)

ROOTEDTREE=${OUTPUT}/rooted_tree.nwk
STATS=${OUTPUT}/stats.csv
DATA=${OUTPUT}/data.csv
PLOT=${OUTPUT}/divergence_versus_time

# Root tree
Rscript ${BDSRC}/root_and_regress.R --runid=${RUNID} --tree=${TREE} --info=${INFO} \
	--rootedtree=${ROOTEDTREE} --data=${DATA} --stats=${STATS} &&

# Plot
Rscript ${BDSRC}/plot_divergence_vs_time.R --rootedtree=${ROOTEDTREE}  --info=${INFO} \
	--stats=${STATS} --plotprefix=${PLOT}
	
	