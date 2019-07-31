library(optparse)
library(ape)
library(tidyverse)
library(parallel)
library(magrittr)

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

ed.data <- mclapply(
	1:length(trees),
	function(i) {
		tree <- trees[[i]]
		
		tree$edge.length <- sapply(
			tree$edge[ ,1],
			. %>%
				equals(tree$edge[ ,1]) %>%
				sum
		)
		
		edge.length <- node.depth.edgelength(tree)
		
		ed <- sapply(
			info$FULLSEQID,
			function(x) {
				1 / edge.length[which(x == tree$tip.label)]
			}
		)
		
		ed.mean <- sapply(
			unique(info$TYPE),
			function(x) {
				mean(ed[x == info$TYPE])
			}
		)
		
		data <- data.frame(
			Treeid=names(trees)[i],
			Type=unique(info$TYPE),
			ED=ed.mean
		)
		
		write.csv(
			data,
			paste(output.file, names(trees)[i], "csv", sep="."),
			row.names=F)
		
		data
	},
	mc.cores=threads
) %>% 
	do.call(rbind, .)

write.csv(ed.data, output.file, row.names=F)