library(ape)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")

get.val <- function(x, default) if (is.null(x)) default else x

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

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--rootedtree", type='character')
op <- add_option(op, "--ogr", type='logical', action='store_true')
op <- add_option(op, "--useall", type='logical', action='store_true')
op <- add_option(op, "--real", type='logical', action='store_true')
op <- add_option(op, "--method", type='character')
op <- add_option(op, "--ogrname", type='character')
op <- add_option(op, "--usedups", type='logical', action='store_true')
op <- add_option(op, "--threads", type='numeric')
op <- add_option(op, "--weight", type='character')
op <- add_option(op, "--alltraining", type='logical', action='store_true')
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
ogr <- get.val(args$ogr, F)
use.all <- get.val(args$useall, F)
use.date <- get.val(args$real, F)
method <- get.val(args$method, 'correlation')		# 'correlation', 'rms' or 'rsquared'
ogr.name <- get.val(args$ogrname, "REFERENCE")
use.dups <- get.val(args$usedups, F)
threads <- get.val(args$threads, 1)
weight <- get.val(args$weight, NA)
all.training <- args$alltraining

use.rtt <- as.numeric(!ogr) * (1 + as.numeric(use.all))	# 0 = no, 1 = yes (only plasma), 2 = yes (all)

tree <- tree.read(tree.file)

if (use.rtt > 0) {
	info.all <- read.csv(info.file, stringsAsFactors=F)
	info <- info.all[match(tree$tip.label, info.all$FULLSEQID), ]
	
	if (use.dups) {
		plasma.dates <- lapply(tree$tip.label, function(x) with(subset(info.all, DUPLICATE == x), if (use.date) as.numeric(as.Date(COLDATE)) else COLDATE))
		
		if (!is.na(weight)) {
			weights <- lapply(tree$tip.label, function(x) info.all[info.all$DUPLICATE == x, weight])
		} else {
			weights <- lapply(tree$tip.label, function(x) rep(1, sum(info.all$DUPLICATE == x)))
		}
	} else {
		plasma.dates <- if (use.date) as.numeric(as.Date(info$COLDATE)) else info$COLDATE
		
		if (!is.na(weight))
			weights <- info[, weight]
	}
}
	
if (use.rtt == 1) {
	if (use.dups) {
		weights <- lapply(1:length(tree$tip.label), function(i) {
			w <- weights[[i]]
			w[with(subset(info.all, DUPLICATE == tree$tip.label[i]), CENSORED > 0 | (!all.training & CENSORED != 0))] <- 0
			w
		})
	} else {
		tip.type <- info$CENSORED
		filter <- tip.type > 0 | (!all.training & tip.type != 0)
		if (is.na(weight))
			plasma.dates[filter] <- NA
		else
			weights[filter] <- 0
	}
}
	
if (use.rtt > 0) {
	if (!use.dups && is.na(weight)) {
		tree <- rtt(tree, plasma.dates, ncpu=threads, objective=method, opt.tol=1e-8)
	} else {
		source(file.path(bd.src, "rtt.R"))
		tree <- rtt(tree, plasma.dates, weights=weights, ncpu=threads, objective=method, opt.tol=1e-8)
	}
}

tree$node.label <- paste0("N.", 1:tree$Nnode)

write.tree(tree, rooted.tree.file)
