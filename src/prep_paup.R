library(ape)
library(seqinr)
library(treeio)
library(optparse)

args.all <- commandArgs(trailingOnly = F)
source.dir <- "../../src"
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
source(file.path(source.dir, 'read.info.R'), chdir=T)

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--nexus", type='character')
op <- add_option(op, "--censor", type='integer', default=0)
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
fasta.file <- args$fasta
nexus.file <- args$nexus
censor <- args$censor

tree <- ape::read.tree(tree.file)
f <- seqinr::read.fasta(fasta.file)
f <- f[match(tree$tip.label, names(f))]
info <- read.info(info.file, tree$tip.label)

keep <- info$CENSORED == censor
f <- f[keep]
info <- info[keep, ]
tree <- ape::drop.tip(tree, which(!keep))

data <- data.frame(Date=info$COLDATE)
rownames(data) <- tree$tip.label

tree.data <- treeio::tree.data(