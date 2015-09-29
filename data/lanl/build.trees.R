library(ape)
source('../common/raxml.R')
source('../common/fasttree.R')

patients <- c(
	10137, 10138, 10586, 10769, 10770, 
	13889, 16616, 16617, 16618, 16619, 
	34375, 34382, 34391, 34393, 34396, 
	34397, 34399, 34400, 34405, 34408, 
	34410, 34411, 820, 821, 822, 824, 
	825)


for(p in patients) {
	dna <- read.FASTA(sprintf('aligned/patient_%d.fasta', p))
	tree <- fasttree(dna)
	# tree <- raxml(dna, N=100, parsimony.seed=10000, bootstrap.seed=1000, outgrp="REFERENCE")
	write.tree(tree, sprintf('tree/patient_%d.tre', p))
}