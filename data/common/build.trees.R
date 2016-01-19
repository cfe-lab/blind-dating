library(ape)
source('../common/raxml.R')
source('../common/fasttree.R')

args <- commandArgs(trailingOnly = T)

prefix <- args[1]
use.raxml <- as.integer(args[2])

patients <- read.table('patients_list.txt', col.names = c("A", "B"))

apply(patients, 1, function (p) {
	dna <- read.FASTA(p["B"])
	if (use.raxml)
		tree <- raxml(dna, N=100, parsimony.seed=10000, bootstrap.seed=1000)
	else
		tree <- fasttree(dna)
	write.tree(tree, sprintf('trees/%s_%d.tre', prefix, as.integer(p["A"])))
})
