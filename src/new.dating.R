library(ape)
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

concord <- function(x, y) {
	mu.x <- sum(x) / length(x)
	mu.y <- sum(y) / length(y)
	s.x <- sum((x - mu.x)^2) / length(x) 
	s.y <-  sum((y - mu.y)^2) / length(y)
	s.xy <- sum((x - mu.x) * (y - mu.y)) / length(y)
	
	2 * s.xy / (s.x + s.y + (mu.x - mu.y)^2)
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--timetree", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--stats", type='character')
op <- add_option(op, "--patid", type='character')
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
	settings.filter <- unlist(
		lapply(
			op@options,
			function(x)
				settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]
		)
	)
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}

tree.file <- args$tree
info.file <- args$info
time.tree.file <- args$timetree
data.file <- args$data
stats.file <- args$stats
pat.id <- args$patid
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
	date=dates,
	dist=node.depth.edgelength(tree)[1:n.tips],
	weight=1
)

total.node.order <- get.node.order(tree)
reach <- unlist(lapply(training.tips, can.reach, tree=tree))
node.order <- c(
	total.node.order[total.node.order %in% reach],
	rev(total.node.order[! total.node.order %in% reach])
)

dates.for.nd <- data$date
dates.for.nd[-training.tips] <- NA

g <- lm(dist ~ date, data, subset=censored == 0)
g.null <- lm(dist ~ date, data, subset=censored == 0)

mu <- coef(g)[[2]]
		
node.dates <- estimate.dates(
	tree,
	dates.for.nd,
	mu=mu,
	node.mask=training.tips,
	node.order=node.order,
	max.date=max(data$date),
	nsteps=nsteps,
	show.steps=if (verbose) 100 else 0,
	lik.tol=lik.tol
)
tree.lik <- logLik.phylo.strict(tree, node.dates, mu)

time.tree <- tree
time.tree$edge.length <- apply(
	tree$edge,
	1,
	function(x)
		node.dates[x[2]] - node.dates[x[1]]
)
write.tree(time.tree, time.tree.file)

data <- as.data.frame(
	cbind(
		data,
		est.date=node.dates[1:n.tips],
		date.diff=node.dates[1:n.tips] - data$date
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
		"Weight",
		"Estimated Date",
		"Date Difference"
	),
	row.names=F,
	sep=","
)

stats <- data.frame(
	pat=pat.id,
	samples.rna=sum(data$censored == training),
	samples.dna=sum(data$censored == 1),
	samples.total=n.tips,
	rna.time.points=length(unique(data[data$censored == training, 'date'])),
	dna.time.points=length(unique(data[data$censored == 1, 'date'])),
	total.time.points=length(unique(data$date)),
	min.rna.time.point=min(data[data$censored == training, 'date']),
	max.rna.time.point=max(data[data$censored == training, 'date']),
	min.dna.time.point=min(data[data$censored == 1, 'date']),
	max.dna.time.point=max(data[data$censored == 1, 'date']),
	min.time.point=min(data$date),
	max.time.point=max(data$date),
	Likelihood=tree.lik,
	LM.fit=as.numeric(AIC(g.null) - AIC(g) > 10),
	mu=mu,
	root.date=node.dates[n.tips + 1],
	RMSD=sqrt(sum(data$date.diff^2) / nrow(data)),
	MAE=sum(abs(data$date.diff)) / nrow(data),
	concord=concord(data$date, data$est.date),
	cens.RMSD=sqrt(
		sum(data$date.diff[data$censored == 1]^2) /
			sum(data$censored == 1)
	),
	cens.MAE=sum(abs(data$date.diff[data$censored == 1])) / 
		sum(data$censored == 1),
	cens.concord=with(
		subset(data, censored == 1),
		concord(date, est.date)
	)
)
stats.col.names <- c(
	"Patient",
	"Training Samples",
	"Censored Samples",
	"Total Samples",
	"Training Time Points",
	"Censored Time Points",
	"Total Time Points",
	"Minimum Training Time Point",
	"Maximum Training Time Point",
	"Minimum Censored Time Point",
	"Maximum Censored Time Point",
	"Minimum Time Point",
	"Maximum Time Point",
	"Likelihood",
	"LM Fit",
	"Model Slope",
	"Estimated Root Date",
	"Total RMSD",
	"Total MAE",
	"Total Concordance",
	"Censored RMSD",
	"Censored MAE",
	"Censored Concordance"
)
write.table(
	stats,
	stats.file,
	col.names=stats.col.names,
	row.names=F,
	sep=","
)