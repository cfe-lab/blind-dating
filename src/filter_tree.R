library(optparse)
library(ape)
library(magrittr)
library(parallel)

op <- OptionParser()
op <- add_option(op, "--trees", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--outputtrees", type='character')
op <- add_option(op, "--threads", type='numeric', default=2)
op <- add_option(op, "--dir", type='logical', action='store_true', default=F)
op <- add_option(op, "--newick", type='logical', action='store_true', default=F)
op <- add_option(op, "--nexus", type='logical', action='store_true', default=F)
op <- add_option(op, "--treetype", type='character', default='nexus')
args <- parse_args(op)

tree.file <- args$trees
info.file <- args$info
output.tree.file <- args$outputtrees
tree.type <- args$treetype
threads <- args$threads

is.dir <- args$dir
is.newick <- args$newick
is.nexus <- args$nexus

if (is.dir) {
	tree.type <- 'dir'
} else if (is.newick) {
	tree.type <- 'newick'
} else if (is.nexus) {
	tree.type <- 'nexus'
}

if (tree.type == 'dir') {
	tree.files <- dir(tree.file) 
	trees <- paste(tree.file, tree.files, sep="/") %>%
		lapply(read.tree)
	
	names(trees) <- gsub(".nwk", "", tree.files)
} else if (tree.type == 'newick') {
	trees <- read.tree(tree.file)
	if (class(trees) == "phylo")
		trees <- list(trees)
} else {
	trees <- read.nexus(tree.file)
}

info <- read.csv(info.file, stringsAsFactors=F) %>%
	subset(
		CENSORED > 0 & DUPLICATE %in% trees[[1]]$tip.label
	)

tips <- split(info, info$DUPLICATE) %>%
	mclapply(
		function(x) {
			split(x, x$TYPE) %>%
				lapply(
					function(y) {
						y[1, ]
					}
				) %>%
				do.call(rbind, .) %$%
				FULLSEQID %>%
				paste0(":0", collapse=",") %>%
				paste0("(", ., ");") %>%
				read.tree(text=.)
		},
		mc.cores=threads
	)

dup.trees <- mclapply(
	trees,
	keep.tip,
	names(tips),
	mc.cores=threads
) %>%
	mclapply(
		function(dup.tree) {
			for (i in 1:length(tips)) {
				dup.tree <- bind.tree(
					dup.tree,
					tips[[i]],
					where=which(
						dup.tree$tip.label == names(tips)[i]
					)
				)
			}
			dup.tree
		},
		mc.cores=threads
	)

write.nexus(dup.trees, file=output.tree.file)