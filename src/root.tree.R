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

source(file.path(source.dir, 'read.info.R'), chdir=T)

args <- commandArgs(trailingOnly = T)

tree.file <- args[1]
info.file <- args[2]
rooted.tree.file <- args[3]
use.rtt <- as.integer(args[4])		# 0 = no, 1 = yes (only plasma), 2 = yes (all)
	
tree.read <- function(tr) {	
	tree <- read.tree(paste(tr, sep='/'))
		
	if (sum("REFERENCE" == tree$tip.label) > 0) {
		if (use.rtt)			
			drop.tip(tree, "REFERENCE")
		else
			drop.tip(root(tree, "REFERENCE"), "REFERENCE")
	}
	else
		tree
			
}

tree <- tree.read(tree.file)
info <- read.info(info.file, tree$tip.label)
plasma.dates <- as.numeric(as.Date(info$COLDATE))
tip.type <- info$CENSORED
	
if (use.rtt == 1)
	plasma.dates[tip.type == 1] <- NA
	
if (use.rtt)
	tree <- rtt(tree, plasma.dates, objective='rms', opt.tol=1e-8)

tree$node.label <- paste0("N.", 1:tree$Nnode)

write.tree(tree, rooted.tree.file)