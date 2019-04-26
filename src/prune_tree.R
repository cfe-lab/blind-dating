library(ape)
library(seqinr)
library(optparse)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--prunetree", type='character')
args <- parse_args(op)

fasta.file <- args$fasta
tree.file <- args$tree
prune.tree.file <- args$prunetree

tree <- read.tree(tree.file)
f <- read.fasta(fasta.file)

prune.tree <- keep.tip(tree, names(f))

write.tree(prune.tree, prune.tree.file)