#!/bin/sh

# Do the analysis on all the output trees
for f in tree/*; do 
	g=${f#tree/} ;
	./analyze.tree.R ${f} 0 analysis/${g%.tre}.pdf
	./analyze.tree.R ${f} 1 analysis/rtt_${g%.tre}.pdf
done
