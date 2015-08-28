#!/bin/sh

# Before this is run, make sure that all the "good" trees have
# been put into /trees.good/
r -f analyze.all.tree.R

# Old: Do the analysis on all the output trees
#for f in tree/*; do 
#	g=${f#tree/} ;
#	r -f analyze.tree.R ${f} 0 analysis/${g%.tre}.pdf
#	r -f analyze.tree.R ${f} 1 analysis/rtt_${g%.tre}.pdf
#done
