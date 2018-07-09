library(ape)
library(optparse)

source("/opt/blind-dating/rtt.R", chdir=T)

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--rootedtree", type='character')
op <- add_option(op, "--ogr", type='logical', action='store_true', default=F)
op <- add_option(op, "--useall", type='logical', action='store_true', default=F)
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--method", type='character', default='correlation')
op <- add_option(op, "--ogrname", type='character', default="REFERENCE")
op <- add_option(op, "--usedups", type='logical', action='store_true', default=F)
op <- add_option(op, "--threads", type='numeric', default=1)
op <- add_option(op, "--freqweights", type='logical', action='store_true', default=F)
op <- add_option(op, "--settings", type='character', default=NA)
args <- parse_args(op)

settings.file <- args$settings
if (!is.na(settings.file)) {
	settings <- readLines(settings.file)
	settings.filter <- unlist(lapply(op@options, function(x) settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]))
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}

tree.file <- args$tree
info.file <- args$info
rooted.tree.file <- args$rootedtree
use.rtt <- 	as.numeric(!args$ogr) * (1 + as.numeric(args$useall))	# 0 = no, 1 = yes (only plasma), 2 = yes (all)
use.date <- args$real
method <- args$method		# 'correlation', 'rms' or 'rsquared'
ogr.name <- args$ogrname
use.dups <- args$usedups
threads <- args$threads
freq.weights <- args$freqweights
	
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
	info.all <- read.csv(info.file, stringsAsFactors=T)
	info <- info.all[match(tree$tip.label, info.all$FULLSEQID), ]
	tip.type <- info$CENSORED
	weights <- NA
	
	if (use.dups) {
		plasma.dates <- lapply(tree$tip.label, function(x) with(subset(info.all, DUPLICATE == x), if (use.date) as.numeric(as.Date(COLDATE)) else COLDATE))
		
		if (use.dups) {
			weights <- lapply(tree$tip.label, function(x) with(subset(info.all, DUPLICATE == x), COUNT))
		}
	} else {
		plasma.dates <- if (use.date) as.numeric(as.Date(info$COLDATE)) else info$COLDATE
		
		if (use.dups) {
			weights <- info$COUNT
		}
	}
}
	
if (use.rtt == 1) {
	plasma.dates[tip.type != 0] <- 0
	weights[tip.type != 0] <- 0
}
	
if (use.rtt)
	tree <- rtt(tree, plasma.dates, weights=weights, ncpu=threads, objective=method, opt.tol=1e-8)

tree$node.label <- paste0("N.", 1:tree$Nnode)

write.tree(tree, rooted.tree.file)
