library(ape)
library(optparse)

source("~/git/node.dating/R/node.dating.R")

get.val <- function(x, default) if (is.null(x)) default else x

make.mean.row <- function(label, data, info, weight) {
	info.dup <- subset(info, DUPLICATE == label)
	data.dup <- subset(data, tip.label == label)
	
	if (data.dup$censored <= 0)
		info.dup <- subset(info.dup, CENSORED <= 0)
	
	data.dup$date[1] <- if (!is.na(weight)) {
		weighted.mean(info.dup$COLDATE, info.dup[weight, ])
	} else {
		mean(info.dup$COLDATE)
	}
	
	data.dup[1, ]
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--stats", type='character')
op <- add_option(op, "--output", type='character')
op <- add_option(op, "--plotgroups", type='logical', action='store_true')
op <- add_option(op, "--plotdups", type='logical', action='store_true')
op <- add_option(op, "--nsteps", type='numeric')
op <- add_option(op, "--weight", type='character')
op <- add_option(op, "--real", type='logical', action='store_true')
op <- add_option(op, "--freeroot", type='logical', action='store_true')
op <- add_option(op, "--settings", type='character', default=NA)
args <- parse_args(op)

settings.file <- args$settings
if (!is.na(settings.file)) {
	settings <- readLines(settings.file)
	settings.filter <- unlist(
		lapply(op@options, function(x)
			settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]
		)
	)
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}

tree.file <- args$tree
data.file <- args$data
info.file <- args$info
stats.file <- args$stats
output.file <- args$output
use.dates <- get.val(args$real, FALSE)
free.root <- get.val(args$freeroot, FALSE)

args$plotgroups <- get.val(args$plotgroups, FALSE)
args$plotdups <- get.val(args$plotdups, FALSE)
args$nsteps <- get.val(args$nsteps, 1000)
args$weight <- get.val(args$weight, NA)

tree <- read.tree(tree.file)
data <- read.csv(data.file)
info <- read.csv(info.file)
stats <- read.csv(stats.file)

colnames(data) <- c(
	"tip.label",
	"type",
	"censored",
	"date",
	"dist",
	"weight",
	"est.date",
	"date.diff"
)
data <- data[match(tree$tip.label, data$tip.label), ]

if (args$plotgroups) {
	info.groups <- subset(info, DUPLICATE %in% data$tip.label)
	info.groups <- info.groups[order(!info.groups$FULLSEQID %in% data$tip.label), ]
	info.split <- split(info.groups, info.groups$TYPE)
	info <- as.data.frame(do.call(rbind, lapply(info.split, function(x) {
		x[sapply(unique(x$DUPLICATE), function(y) which(x$DUPLICATE == y)[1]), ]
	})))
} else if (args$plotdups) {
	info <- subset(info, DUPLICATE %in% data$tip.label)
} else {
	info <- subset(info, FULLSEQID %in% tree$tip.label)
}

info$COLDATE <- as.numeric(as.Date(info$COLDATE))

data.mean <- do.call(
	rbind,
	lapply(tree$tip.label, make.mean.row, data, info, args$weight)
)

mu <- stats$Model.Slope

node.dates <- tryCatch(
	{
		if (free.root)
			stop("free root")
		estimate.dates(
			tree,
			c(data.mean$date, stats$Estimated.Root.Date, rep(NA, tree$Nnode - 1)),
			mu,
			node.mask=1:(length(tree$tip.label) + 1),
			lik.tol=0,
			nsteps=args$nsteps,
			show.steps=100,
			opt.tol=1e-16
		)
	},
	error=function(e) {
		cat("Error: ")
		message(e)
		cat("\n")
		estimate.dates(
			tree,
			data.mean$date,
			mu,
			lik.tol=0,
			nsteps=args$nsteps,
			show.steps=100,
			opt.tol=1e-16
		)
	}
)

write.csv(data.frame(dates=node.dates), output.file)