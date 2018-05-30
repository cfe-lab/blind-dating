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

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--rootedtree", type='character')
op <- add_option(op, "--ogr", type='logical', action='store_true', default=F)
op <- add_option(op, "--useall", type='logical', action='store_true', default=F)
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--method", type='character', default='correlation')
op <- add_option(op, "--ogrname", type='character', default="REFERENCE")
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
rooted.tree.file <- args$rootedtree
use.rtt <- 	as.numeric(!args$ogr) * (1 + as.numeric(args$usall))	# 0 = no, 1 = yes (only plasma), 2 = yes (all)
use.date <- args$real
method <- args$method		# 'correlation', 'rms' or 'rsquared'
ogr.name <- args$ogrname
	
tree.read <- function(tr) {	
	tree <- read.tree(paste(tr, sep='/'))
		
	if (any(tree$tip.label == ogr.name)) {
		if (use.rtt > 0)			
			drop.tip(tree, ogr.name)
		else
			drop.tip(root(tree, ogr.name), ogr.name)
	}
	else
		tree
			
}

tree <- tree.read(tree.file)

if (use.rtt > 0) {
	info <- read.info(info.file, tree$tip.label)
	plasma.dates <-  if (use.date) as.numeric(as.Date(info$COLDATE)) else info$COLDATE
	tip.type <- info$CENSORED
}
	
if (use.rtt == 1)
	plasma.dates[tip.type != 0] <- NA
	
if (use.rtt)
	tree <- rtt(tree, plasma.dates, objective=method, opt.tol=1e-8)

tree$node.label <- paste0("N.", 1:tree$Nnode)

write.tree(tree, rooted.tree.file)
