#!/usr/bin/Rscript

library(ape)
library(phylobase)
library(ggtree)
library(optparse)

get.date <- function(date) {
	as.character(as.Date(date, origin=as.Date("1970-01-01")))
}

get.child.edges <- function(tree, node) {
	e <- which(tree$edge[, 1] %in% node)
	if (length(e) > 0)
		c(e, get.child.edges(tree, tree$edge[e, 2]))
	else
		NULL
}

apply.axes <- function(p, flipped, scaled) {
	if (flipped) {
		p <- p + scale_x_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by), limits=c(dist.min, dist.max))
		if (scaled) {
			p <- p + y.scale +
				geom_abline(intercept=-stats[, "Model.Intercept"] / stats[, "Model.Slope"], slope=1 / stats[, "Model.Slope"], color="#0060b0b0", linetype=2) +
				geom_hline(yintercept=THERAPY_START, colour="#60600080", linetype=2)
			if (!is.na(pat.id2))
				p <- p + geom_abline(intercept=-stats.2[, "Model.Intercept"] / stats.2[, "Model.Slope"], slope=1 / stats.2[, "Model.Slope"], color="#003058b0", linetype=2)
		} else
			p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
	} else {
		p <- p + scale_y_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by), limits=c(dist.min, dist.max))
		if (scaled) {
			p <- p + x.scale +
				geom_abline(intercept=stats[, "Model.Intercept"], slope=stats[, "Model.Slope"], color="#0060b0b0", linetype=2) +
				geom_vline(xintercept=THERAPY_START, colour="#60600080", linetype=2)
			if (!is.na(pat.id2))
				p <- p + geom_abline(intercept=stats.2[, "Model.Intercept"], slope=stats.2[, "Model.Slope"], color="#003058b0", linetype=2)
		} else
			p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
	}
	p
}

apply.theme <- function(p, flipped=F, scaled=T) {
	apply.axes(p + theme_bw() +
		theme(
			legend.position=c(.01, .99),
			legend.justification=c(0, 1),
			legend.spacing=unit(5, 'points'),
			legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
			legend.background=element_blank(),
			legend.box.background=element_blank(),
			panel.grid.major=element_blank(),
		    panel.grid.minor=element_blank()
		) +
		scale_colour_manual(name="", breaks=type.break, labels=type.label, values=type.value, limits=type.break) +
		scale_shape_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(16, 18), limits=c(0, 1)) +
		scale_size_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(2.5, 3), limits=c(0, 1)) +
		guides(colour=type.guide),
		flipped, scaled)
}

LIK_TOL <- 1e-1

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--patid", type='character')
op <- add_option(op, "--patid2", type='character')
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--distmax", type='double')
op <- add_option(op, "--distmin", type='double', default=-.01)
op <- add_option(op, "--distby", type='double')
op <- add_option(op, "--yearstart", type='character')
op <- add_option(op, "--yearend", type='character')
op <- add_option(op, "--therapy", type='character', default="0000-01-01")
args <- parse_args(op)

tree.file <- args$tree
pat.id <- args$patid
pat.id2 <- if ("patid2" %in% names(args)) args$patid2 else NA
use.real <- args$real # 0 Training/Censored, 1 RNA/DNA
dist.min <- args$distmin
dist.max <- args$distmax
dist.by <- args$distby
if (use.real) {
	year.start <- args$yearstart
	year.end <- args$yearend
}
THERAPY_START <- as.numeric(as.Date(args$therapy))

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")
if (!is.na(pat.id2)) {
	data.2.file <- paste0("stats/", pat.id2, ".data.csv")
	stats.2.file <- paste0("stats/", pat.id2, ".stats.csv")
}
pdf.file <- paste0("plots/", pat.id, ".pdf")
pdf.disttree.file <- paste0("plots/", pat.id, ".disttree.pdf")
pdf.tree.file <- paste0("plots/", pat.id, ".tree.pdf")
pdf.hist.file <- paste0("plots/", pat.id, ".hist.pdf")

tree <- ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
stats <- read.csv(stats.file, stringsAsFactors=F)
if (!is.na(pat.id2)) {
	data.2 <- read.csv(data.2.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
	stats.2 <- read.csv(stats.2.file, stringsAsFactors=F)
		
	mu <- rep(stats[, "Model.Slope"], nrow(tree$edge))
	clade <- get.child.edges(tree, getMRCA(tree, data$tip.label[data$censored == -1]))
	mu[clade] <- stats.2[, "Model.Slope"]
	clade.tips <- data$edge[clade, 2]
	clade.tips <- clade.tips[clade.tips <= length(tree$tip.label)]
		
	data[clade.tips, ] <- data.2[clade.tips, ]
	data$censored[data$censored == -1] <- 0
	
} else
	mu <- stats[, "Model.Slope"]

node.dates <- estimate.dates(tree, data$date, mu, lik.tol=LIK_TOL, show.steps=1000, nsteps=0)

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type=rep("NODE", tree$Nnode), censored=rep(NA, tree$Nnode), date=rep(NA, tree$Nnode), dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], est.date=rep(NA, tree$Nnode), date.diff=rep(NA, tree$Nnode))), node.date=node.dates))

ptree <- phylo4d(tree, all.data=data.all)

if (use.real) {
	date.ticks <- as.character(seq(floor(as.numeric(year.start) / 5) * 5 + 5, as.numeric(year.end), by=5))
	
	x.scale <- scale_x_continuous(name="Year", breaks=as.numeric(as.Date(paste0(date.ticks, "-01-01"))), labels=date.ticks, limits=as.numeric(as.Date(paste0(c(year.start, year.end), "-01-01"))))
	y.scale <- scale_y_continuous(name="Year", breaks=as.numeric(as.Date(paste0(date.ticks, "-01-01"))), labels=date.ticks, limits=as.numeric(as.Date(paste0(c(year.start, year.end), "-01-01"))))
	type.guide <- guide_legend(override.aes=list(shape=15, size=5))
	if (length(unique(data$type)) == 4) {
		type.break <- c("PLASMA", "PBMC", "PBMC (cultured)", "WHOLE BLOOD")
		type.label <- c("Plasma RNA", "PBMC DNA", "Cultured PBMC DNA", "Whole Blood DNA")
		type.value <- c('black', 'red', 'cyan', 'darkgreen')
	} else {
		type.break <- c("PLASMA", "PBMC")
		type.label <- c("Plasma RNA", "PBMC DNA")
		type.value <- c('black', 'red')
	}
} else {
	x.scale <- scale_x_continuous(name="Days post simulation")
	y.scale <- scale_y_continuous(name="Days post simulation")
	type.guide <- 'legend'
	type.break <- c("Training", "Censored")
	type.label <- c("Training", "Censored")
	type.value <- c('black', 'red')
}

pdf(pdf.file)
apply.theme(ggplot(data) +
	geom_point(aes(x=date, y=dist, colour=type, shape=factor(censored), size=factor(censored))))
dev.off()

pdf(pdf.disttree.file)
apply.theme(ggtree(ptree, colour="#49494980", size=.5, ladderize=F) +
	geom_tiplab(colour="#49494980", angle=90, hjust=-.1, size=1.5) + 
	geom_tippoint(aes(colour=type, shape=factor(censored), size=factor(censored))) +
	coord_flip(), flipped=T, scaled=F)
dev.off()

pdf(pdf.tree.file)
apply.theme(ggtree(ptree, yscale="node.date", colour="#49494980", size=.3, ladderize=F) +
	geom_tippoint(aes(colour=type, shape=factor(censored), size=factor(censored))) +
	coord_flip(), flipped=T)
dev.off()

pdf(pdf.hist.file)
ggplot(subset(data, censored == 1)) +
	geom_histogram(aes(x=date.diff, fill=factor(type, levels=type.break)), bins=10) +
	scale_fill_manual(name="", breaks=type.break, labels=type.label, values=type.value, limits=type.break) +
	scale_x_continuous("Difference in Estimated and Collection Date (days)") +
	scale_y_continuous("Count") +
	theme_bw() +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)
dev.off()