library(ape)
library(ggtree)
library(phylobase)
library(optparse)
source("~/git/node.dating/R/node.dating.R")

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "read.info.R"), chdir=T)

get.val <- function(x, default) if (is.null(x)) default else x

can.reach <- function(tree, node) {
	root <- length(tree$tip.label) + 1
	nodes <- node
	
	while (node != root) {
		node <- tree$edge[tree$edge[, 2] == node, 1]
		
		nodes <- c(nodes, node)
	}
	
	nodes
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--timetree", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--stats", type='character')
op <- add_option(op, "--real", type='logical', action='store_true')
op <- add_option(op, "--training", type='numeric')
op <- add_option(op, "--liktol", type='numeric')
op <- add_option(op, "--steps", type='numeric')
op <- add_option(op, "--verbose", type='logical', action='store_true')
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
time.tree.file <- args$timetree
data.file <- args$data
stats.file <- args$stats
real <- get.val(args$real, F)
training <- get.val(args$training, 0)
lik.tol <- get.val(args$lik.tol, 0)
nsteps <- get.val(args$steps, 1000)
verbose <- get.val(args$verbose, F)

tree <- read.tree(tree.file)
info <- read.info(info.file, tree$tip.label)

n.tips <- length(tree$tip.label)

dates <- as.numeric(
	if (real) {
		as.Date(info$COLDATE)
	} else {
		info$COLDATE
	}
)
training.tips <- which(info$CENSORED == training)

data <- data.frame(
	tip.label=tree$tip.label,
	type=info$TYPE,
	censored=info$CENSORED,
	date=info$COLDATE,
	dist=node.depth.edgelength(tree)[1:n.tips]
)

total.node.order <- get.node.order(tree)
reach <- unlist(lapply(training.tips, can.reach, tree=tree))
node.order <- c(
	total.node.order[total.node.order %in% reach],
	rev(total.node.order[! total.node.order %in% reach])
)

dates.for.nd <- data$date
dates.for.nd[-training] <- NA

model <- lm(date ~ dist, data=data, subset=censored==training)

mu <- estimate.mu(tree, dates.for.nd)

node.dates <- estimate.dates(
	tree,
	dates.for.nd,
	mu=mu,
	node.mask=training,
	node.order=node.order,
	max.date=max(data$date),
	nsteps=1000
)

time.tree <- tree
time.tree$edge.length <- apply(tree$edge, 1, function(x) node.dates[2] - node.dates[1])
write.tree(time.tree, time.tree.file)

data <- as.data.frame(
	cbind(
		data,
		est.date=new.dates[1:n.tips],
		date.diff=new.dates[1:n.tips] - data$date
	)
)
write.table(
	data,
	data.file,
	col.names=c(
		"ID",
		"Type",
		"Censored",
		"Collection Date",
		"Divergence",
		"Estimated Date",
		"Date Difference"
	),
	row.names=F,
	sep=","
)