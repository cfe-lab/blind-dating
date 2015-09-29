library(ape)
source('../common/raxml.R')
source('../common/fasttree.R')

patients <- c(1, 19, 58, 128, 129, 130, 131, 132, 133, 427, 600)

for(p in patients) {
	dna <- read.FASTA(sprintf('aligned_tweak/patient_%d.fasta', p))
	tree <- fasttree(dna)
#	tree <- raxml(dna, N=1000, parsimony.seed=10000, bootstrap.seed=1000, outgrp="REFERENCE")
	write.tree(tree, sprintf('tree/patient_%d.tre', p))
}