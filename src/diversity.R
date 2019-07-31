library(optparse)
library(ape)
library(tidyverse)
library(parallel)

op <- OptionParser()
op <- add_option(op, "--trees", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--output", type='character')
op <- add_option(op, "--threads", type='numeric', default=2)
args <- parse_args(op)

tree.file <- args$trees
info.file <- args$info
output.file <- args$output
threads <- args$threads

trees <- read.nexus(tree.file)

info <- read.csv(info.file, stringsAsFactors=F)
info <- info[match(trees[[1]]$tip.label, info$FULLSEQID), ]
info.s <- split(info, info$TYPE)

pd.data <- mclapply(
	1:length(trees),
	function(i) {
		tree <- trees[[i]]
		
		pd <- sapply(
			info.s,
			function(x) {
				subtree <- keep.tip(tree, x$FULLSEQID)
				
				sum(subtree$edge.length) / length(subtree$tip.label)
			}
		)
		
		data <- data.frame(Treeid=names(trees)[i], Type=names(info.s), PD=pd)
		
		write.csv(data, paste(output.file, names(trees)[i], "csv", sep="."), row.names=F)
		
		data
	},
	mc.cores=threads
) %>% 
	do.call(rbind, .)

write.csv(pd.data, output.file, row.names=F)