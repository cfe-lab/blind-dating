#!/usr/bin/Rscript

library(ape)
library(phylobase)
library(ggtree)
library(optparse)

source("~/git/node.dating/src/node.dating.R")


THERAPY_COLOUR <- "#a0a0a010"
THERAPY_COLOUR2 <- "#c0c0c010"

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

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

get.parent.edges <- function(tree, node) {
	e <- which(tree$edge[, 2] %in% node)
	if (length(e) > 0)
		c(e, get.child.edges(tree, tree$edge[e, 1]))
	else
		NULL
}

apply.axes <- function(p, flipped, scaled) {
	if (scaled) {
		if (flipped) {
			p <- p + 
				y.scale +
				scale_x_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by)) + 
				coord_flip(xlim=c(dist.min, dist.max)) +
				geom_abline(intercept=-stats[, "Model.Intercept"] / stats[, "Model.Slope"], slope=1 / stats[, "Model.Slope"], colour="#0050b0b0", linetype=2)
#				geom_path(aes(y=date, x=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(y=date, x=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (use.real) {
				if (is.na(therapy2))
					p$layers <- c(geom_rect(ymin=as.numeric(THERAPY_START), ymax=Inf, xmin=-Inf, xmax=Inf, fill=THERAPY_COLOUR, linetype=2), p$layers)
				else {
					p$layers <- c(
						geom_rect(ymin=as.numeric(THERAPY_START), ymax=as.numeric(therapyend), xmin=-Inf, xmax=Inf, fill=if (therapyend == therapy2) THERAPY_COLOUR2 else THERAPY_COLOUR, linetype=0),
						geom_rect(ymin=as.numeric(therapy2), ymax=Inf, xmin=-Inf, xmax=Inf, fill=THERAPY_COLOUR, linetype=0),
						p$layers
					)
				}
			}
			if (!is.na(pat.id2))
				p <- p + geom_abline(intercept=-stats.2[, "Model.Intercept"] / stats.2[, "Model.Slope"], slope=1 / stats.2[, "Model.Slope"], colour="#003058b0", linetype=2)
#					geom_path(aes(y=date, x=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(y=date, x=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		} else {
			p <- p + 
				x.scale +
				scale_y_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by)) +
				coord_cartesian(ylim=c(dist.min, dist.max)) +
				geom_abline(intercept=stats[, "Model.Intercept"], slope=stats[, "Model.Slope"], colour="#0050b0b0", linetype=2)
#				geom_path(aes(x=date, y=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(x=date, y=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (use.real) {
				if (is.na(therapy2))
					p$layers <- c(geom_rect(xmin=as.numeric(THERAPY_START), xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR, linetype=0), p$layers)
				else {
					p$layers <- c(
						geom_rect(xmin=as.numeric(THERAPY_START), xmax=as.numeric(therapyend), ymin=-Inf, ymax=Inf, fill=if (therapyend == therapy2) THERAPY_COLOUR2 else THERAPY_COLOUR, linetype=0),
						geom_rect(xmin=as.numeric(therapy2), xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR, linetype=0),
						p$layers
					)
				}
			}
			if (!is.na(pat.id2))
				p <- p + geom_abline(intercept=stats.2[, "Model.Intercept"], slope=stats.2[, "Model.Slope"], colour="#003058b0", linetype=2)
#					geom_path(aes(x=date, y=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(x=date, y=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		}
	} else
		p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
	p
}

apply.theme <- function(p, flipped=F, scaled=T) {
	apply.axes(p + theme_bw() +
		theme(
			text=element_text(size=35),
			axis.text=element_text(size=30, colour='black'),
			legend.text=element_text(size=30),
			legend.position=0,
			legend.justification=c(0, 1),
			legend.spacing=unit(0, 'cm'),
			legend.margin=margin(0, 0, 0, 0, 'cm'),
			legend.key.size=unit(1.2, 'cm'),
			legend.key=element_rect(fill="#00000000", colour="#00000000"),
			legend.background=element_blank(),
			legend.box.background=element_blank(),
			panel.grid.major=element_blank(),
		    panel.grid.minor=element_blank()
		) +
		scale_shape_manual(name="", breaks=type.label, labels=type.label, values=type.value, limits=type.label) +
		scale_colour_manual(name="", breaks=colour.break, labels=colour.label, values=colour.value, limits=colour.break) +
#		scale_size_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(1.67, 2), limits=c(0, 1)) +
		guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1)),
		flipped, scaled)
}

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
op <- add_option(op, "--yearby", type='double', default=5)
op <- add_option(op, "--therapy", type='character', default="0000-01-01")
op <- add_option(op, "--therapy2", type='character', default=NA)
op <- add_option(op, "--therapyend", type='character', default=NA)
op <- add_option(op, "--liktol", type='numeric', default=1e-3)
op <- add_option(op, "--usedups", type='logical', action='store_true', default=F)
args <- parse_args(op)

tree.file <- args$tree
pat.id <- args$patid
pat.id2 <- if ("patid2" %in% names(args)) args$patid2 else NA
use.real <- args$real # 0 Training/Censored, 1 RNA/DNA
dist.min <- args$distmin
dist.max <- args$distmax
dist.by <- args$distby
LIK_TOL <- args$liktol
therapy2 <- args$therapy2
therapyend <- args$therapyend
if (use.real) {
	year.start <- args$yearstart
	year.end <- args$yearend
	THERAPY_START <- as.numeric(as.Date(args$therapy))
} else {
	year.start <- as.numeric(args$yearstart)
	year.end <- as.numeric(args$yearend)
}
if (!is.na(therapy2)) {
	therapy2 <- as.Date(therapy2)
}
if (is.na(therapyend)) {
	therapyend <- therapy2
} else {
	therapyend <- as.Date(therapyend)
}
year.by <- args$yearby

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")
regression.file <- paste0("stats/", pat.id, ".regression.rds")
if (!is.na(pat.id2)) {
	data.2.file <- paste0("stats/", pat.id2, ".data.csv")
	stats.2.file <- paste0("stats/", pat.id2, ".stats.csv")
	regression.2.file <- paste0("stats/", pat.id2, ".regression.rds")
}
pdf.file <- paste0("plots/", pat.id, ".pdf")
pdf.disttree.file <- paste0("plots/", pat.id, ".disttree.pdf")
pdf.tree.file <- paste0("plots/", pat.id, ".tree.pdf")
pdf.hist.file <- paste0("plots/", pat.id, ".hist.pdf")
pdf.histdate.file <- paste0("plots/", pat.id, ".histdate.pdf")

tree <- ape::ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
stats <- read.csv(stats.file, stringsAsFactors=F)
g <- readRDS(regression.file)

ci.dates <- if (use.real) {
	seq(as.numeric(as.Date(paste0(year.start, "-01-01"))), as.numeric(as.Date(paste0(year.end, "-01-01"))), length.out=1000)
} else {
		seq(as.numeric(year.start), as.numeric(year.end), length.out=1000)
}
p <- predict(g, data.frame(date=ci.dates), interval='confidence')
data.ci <- data.frame(date=ci.dates, ci.low=p[,2], ci.high=p[,3])

if (!is.na(pat.id2)) {
	data.2 <- read.csv(data.2.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
	stats.2 <- read.csv(stats.2.file, stringsAsFactors=F)
	g.2 <- readRDS(regression.2.file)
	
	p.2 <- predict(g.2, data.frame(date=ci.dates), interval='confidence')
	data.ci <- as.data.frame(cbind(data.ci, ci.low.2=p.2[,2], ci.high.2=p.2[,3]))
	
	mu <- rep(stats.2[, "Model.Slope"], nrow(tree$edge))
	first.edges <- unique(unlist(lapply(which(data$censored == 0), get.parent.edges, tree=tree)))
	mu[first.edges] <- stats[, "Model.Slope"]
	
#	mu <- rep(stats[, "Model.Slope"], nrow(tree$edge))
	clade <- get.child.edges(tree, getMRCA(tree, data$tip.label[data$censored == -1]))
#	mu[clade] <- stats.2[, "Model.Slope"]
	clade.tips <- tree$edge[clade, 2]
	clade.tips <- clade.tips[clade.tips <= length(tree$tip.label)]
		
	data.old <- data
	data[clade.tips, ] <- data.2[clade.tips, ]
	data[data.old$censored == -1, "censored"] <- -1
	data[data.old$censored == 0, ] <- data.old[data.old$censored == 0, ]
	write.table(data, paste0("stats/", pat.id, ".comb.data.csv"), col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Estimated Date", "Date Difference"), row.names=F, sep=",")
	
	stats$Minimum.Time.Point <- min(stats$Mimimum.Time.Point, stats.2$Minimum.Time.Point)
	stats$Minimum.Training.Time.Point <- min(stats$Mimimum.Training.Time.Point, stats.2$Minimum.Training.Time.Point)
	stats$Maximum.Time.Point <- min(stats$Maximum.Time.Point, stats.2$Maximum.Time.Point)
	stats$Maximum.Training.Time.Point <- min(stats$Maximum.Training.Time.Point, stats.2$Maximum.Training.Time.Point)
	stats$Training.Samples <- stats$Training.Samples + stats.2$Training.Samples
	stats$Censored.RMSD <- sqrt(sum(data[data$censored == 1, "date.diff"]^2)/sum(data$censored == 1))
	stats <- as.data.frame(cbind(stats, stats.2[, c("AIC", "null.AIC", "p.value", "Model.Intercept", "Model.Slope", "Training.RMSE")]))
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
		"AIC",
		"null AIC",
		"p-value",
		"Model Intercept",
		"Model Slope",
		"Model Error",
		"Estimated Root Date",
		"Training RMSE",
		"Censored RMSD",
		"Total RMSD",
		"Training MAE",
		"Censored MAE",
		"Total MAE",
		"Total Concordance",
		"AIC 2",
		"null AIC 2",
		"p-value 2",
		"Model Intercept 2",
		"Model Slope 2",
		"Training RMSE 2"
	)
	write.table(stats, paste0("stats/", pat.id, ".comb.stats.csv"), row.names=F, col.names=stats.col.names, sep=",")

	colour.break <- c(0, -1, 1)
	colour.label <- c("Training", "Training 2", "Censored")
	colour.value <- c('black', 'darkblue', 'red')
} else {
	mu <- stats[, "Model.Slope"]

	colour.break <- c(0, 1)
	colour.label <- c("Training", "Censored")
	colour.value <- c('black', 'red')
}

if (use.real) {
	date.ticks <- as.character(seq(floor(as.numeric(year.start) / year.by) * year.by - year.by, as.numeric(year.end), by=year.by))
	x.scale <- scale_x_continuous(name="Collection Year", breaks=as.numeric(as.Date(paste0(date.ticks, "-01-01"))), labels=date.ticks, limits=as.numeric(as.Date(paste0(c(year.start, year.end), "-01-01"))))
	y.scale <- scale_y_continuous(name="Collection Year", breaks=as.numeric(as.Date(paste0(date.ticks, "-01-01"))), labels=date.ticks, limits=as.numeric(as.Date(paste0(c(year.start, year.end), "-01-01"))))
	x.scale.hist <- scale_x_continuous(name="Collection Year", breaks=as.numeric(as.Date(paste0(seq(as.numeric(year.start), as.numeric(year.end)), "-01-01"))), labels=as.character(seq(as.numeric(year.start), as.numeric(year.end))))
	
	type.break <- c("PLASMA", "PBMC", "PBMC (cultured)", "WHOLE BLOOD", "PBMC (REACTIVE)", "PBMC (CULTURE)")
	type.label <- c("RNA", "DNA", "DNA", "DNA", "RNA", "DNA")
	type.value <- c(16, 5, 5, 5, 16, 5)
} else {
	date.ticks <- as.character(seq(floor(year.start / year.by) * year.by + year.by, year.end, by=year.by))
	x.scale <- scale_x_continuous(name="Days post simulation", breaks=date.ticks, limits=c(year.start, year.end))
	y.scale <- scale_y_continuous(name="Days post simulation", breaks=date.ticks, limits=c(year.start, year.end))
	x.scale.hist <- scale_x_continuous(name="Days post simulation")
	
	type.break <- c("Training", "Censored")
	type.label <- c("Training", "Censored")
	type.value <- c(16, 5)
}

data$type <- type.label[match(data$type, type.break)]

type.mask <- type.label %in% data$type
type.label <- unique(type.label[type.mask])
type.value <- unique(type.value[type.mask])

node.dates <- estimate.dates(tree, c(data$date, stats$Estimated.Root.Date, rep(NA, tree$Nnode - 1)), mu, node.mask=length(tree$tip.label) + 1, lik.tol=0, nsteps=1000, show.steps=100, opt.tol=1e-16)

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type="NODE", censored=NA, date=NA, dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], est.date=NA, date.diff=NA)), node.date=node.dates))

ptree <- phylo4d(tree, all.data=data.all)

pdf(pdf.file)
apply.theme(ggplot(data) +
	geom_point(aes(x=date, y=dist, colour=factor(censored), shape=type), size=6)) +
	theme(legend.position=c(0.02, .98))
dev.off()

pdf(pdf.disttree.file)
apply.theme(ggtree(ptree, colour="#49494980", size=.6, ladderize=T) +
#	geom_tiplab(colour="#49494980", angle=90, hjust=-.1, size=1) + 
	geom_tippoint(aes(colour=factor(censored), shape=type), size=6) +
	geom_treescale(width=0.02, offset=2), flipped=F, scaled=F)
dev.off()

pdf(pdf.tree.file)
apply.theme(ggtree(ptree, yscale="node.date", colour="#6d4d4180", size=.4, ladderize=F) +
	geom_tippoint(aes(colour=factor(censored), shape=type), size=6),
	flipped=T)
dev.off()

if (F) {
pdf(pdf.hist.file)
ggplot(subset(data, censored == 1)) +
	geom_histogram(aes(x=date.diff, fill=factor(type, levels=type.break)), bins=10) +
	scale_fill_manual(name="", breaks=type.break, labels=type.label, values=type.value, limits=type.break) +
	scale_x_continuous("Difference in Estimated and Collection Date (days)") +
	scale_y_continuous("Count") +
	theme_bw() +
	theme(
		text=element_text(size=10),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)
dev.off()

pdf(pdf.histdate.file)
ggplot(subset(data, censored == 1)) +
	geom_histogram(aes(x=est.date, fill=factor(type, levels=type.break)), bins=10) +
	scale_fill_manual(name="", breaks=type.break, labels=type.label, values=type.value, limits=type.break) +
	x.scale.hist +
	scale_y_continuous("Count") +
	theme_bw() +
	theme(
		text=element_text(size=10),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)
dev.off()
}
