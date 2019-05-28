library(optparse)
library(ape)
library(magrittr)
library(parallel)

op <- OptionParser()
op <- add_option(op, "--trees", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--outputtrees", type='character')
op <- add_option(op, "--threads", type='numeric', default=2)
args <- parse_args(op)

tree.file <- args$trees
info.file <- args$info
output.tree.file <- args$outputtrees
threads <- args$threads

trees <- read.nexus(tree.file)

info <- read.csv(info.file) %>%
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