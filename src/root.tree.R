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
use.date <- if (length(args) >= 5) as.integer(args[5]) else 1
method <- if (length(args) >= 6) args[6] else 'correlation'		# 'correlation', 'rms' or 'rsquared'
	
tree.read <- function(tr) {	
	tree <- read.tree(paste(tr, sep='/'))
		
	if (any("REFERENCE" == tree$tip.label)) {
		if (use.rtt > 0)			
			drop.tip(tree, "REFERENCE")
		else
			drop.tip(root(tree, "REFERENCE"), "REFERENCE")
	}
	else
		tree
			
}

tree <- tree.read(tree.file)

if (use.rtt > 0) {
	info <- read.info(info.file, tree$tip.label)
	plasma.dates <-  if (use.date == 1) as.numeric(as.Date(info$COLDATE)) else info$COLDATE
	tip.type <- info$CENSORED
}
	
if (use.rtt == 1)
	plasma.dates[tip.type != 0] <- NA
	
if (use.rtt)
	tree <- rtt(tree, plasma.dates, objective=method, opt.tol=1e-8)

tree$node.label <- paste0("N.", 1:tree$Nnode)

write.tree(tree, rooted.tree.file)
