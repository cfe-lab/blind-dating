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

args <- commandArgs(trailingOnly = T)

fasta.file = args[1]
tree.file = args[2]
exec = args[3]
threads = as.integer(args[4])
seed = as.integer(args[5])
tmp.dir <- if (length(args) > 5) args[6] else NULL

tmp <- !is.null(tmp.dir)

tree <- raxml(fasta.file, N=100, threads=threads, executable=exec, parsimony.seed=seed, bootstrap.seed=seed, clear=!tmp, tmp.dir=if(tmp) tmp.dir else tempdir())
write.tree(tree, tree.file)
