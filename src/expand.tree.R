library(optparse)
library(magrittr)
library(ape)

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--fulltree", type='character')
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
full.tree.file <- args$fulltree

tree <- read.tree(tree.file)
treeLine <- readLines(tree.file)
info <- read.csv(info.file)

treeLineRep <- treeLine
for (x in tree$tip.label) {
	info.dup <- subset(info, DUPLICATE == x)
	
	if (nrow(info.dup) > 1) {
		treeLineRep %<>% gsub(
			paste0(x, ":"),
			paste0(
				"(",
				paste0(
					info.dup$FULLSEQID,
					":0",
					collapse=","
				),
				"):"
			),
			.
		)
	}
}

writeLines(treeLineRep, full.tree.file)
