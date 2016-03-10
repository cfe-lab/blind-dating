library(ape)

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

source(file.path(source.dir, 'raxml.R'), chdir=T)
source(file.path(source.dir, 'fasttree.R'), chdir=T)

args <- commandArgs(trailingOnly = T)

prefix <- args[1]
use.raxml <- as.integer(args[2])

patients <- read.table('patients_list.txt', col.names = c("A", "B"))

apply(patients, 1, function (p) {
	print(p["A"])

	dna <- read.FASTA(p["B"])
	if (use.raxml)
		tree <- raxml(dna, N=100, parsimony.seed=10000, bootstrap.seed=1000)
	else
		tree <- fasttree(dna)
	write.tree(tree, sprintf('trees/%s_%d.tre', prefix, as.integer(p["A"])))
})
