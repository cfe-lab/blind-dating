#!/usr/bin/Rscript

library(ape)
library(phylobase)
library(ggplot2)
library(ggtree)
library(optparse)

source("~/git/node.dating/src/node.dating.R")

MIN_COL_TIME <- 9720
MAX_COL_TIME <- 13371

THERAPY_COLOUR <- "#111111"
THERAPY_LTY <- 3
THERAPY_COLOUR2 <- "#111111"
THERAPY_LTY2 <- 6

TRAINING_COLOUR <- 'black'
TRAINING_COLOUR2 <- 'darkblue'
CENSORED_COLOUR <- 'red'
CENSORED_COLOUR2 <- "#ff8000"

REGRESS_COLOUR <- "#808080b0"
REGRESS_LTY <- 2
REGRESS_COLOUR2 <- "#aaaaaab0"
REGRESS_LTY2 <- 4

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
				x.div + 
				coord_flip(xlim=c(dist.min, dist.max))				
#				geom_path(aes(y=date, x=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(y=date, x=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (!is.na(THERAPY_START)) {
				if (is.na(therapy2))
					p$layers <- c(geom_hline(yintercept=as.numeric(THERAPY_START), linetype=THERAPY_LTY, colour=THERAPY_COLOUR), p$layers)
				else {
					p$layers <- c(
						geom_hline(yintercept=as.numeric(THERAPY_START), linetype=THERAPY_LTY2, colour=THERAPY_COLOUR2),
						geom_hline(yintercept=as.numeric(therapy2), linetype=THERAPY_LTY, colour=THERAPY_COLOUR),
						p$layers
					)
				}
			}
			p$layers <- c(geom_abline(intercept=-stats[, "Model.Intercept"] / stats[, "Model.Slope"], slope=1 / stats[, "Model.Slope"], colour=REGRESS_COLOUR, linetype=REGRESS_LTY, size=regression.size), p$layers)
			if (!is.na(pat.id2))
				p$layers <- c(geom_abline(intercept=-stats.2[, "Model.Intercept"] / stats.2[, "Model.Slope"], slope=1 / stats.2[, "Model.Slope"], colour=REGRESS_COLOUR2, linetype=REGRESS_LTY2, size=regression.size), p$layers)
#					geom_path(aes(y=date, x=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(y=date, x=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		} else {
			p <- p + 
				x.scale +
				y.div +
				coord_cartesian(ylim=c(dist.min, dist.max))				
#				geom_path(aes(x=date, y=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(x=date, y=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (!is.na(THERAPY_START)) {
				if (is.na(therapy2))
					p$layers <- c(geom_vline(xintercept=as.numeric(THERAPY_START), linetype=THERAPY_LTY, colour=THERAPY_COLOUR), p$layers)
				else {
					p$layers <- c(
						geom_vline(xintercept=as.numeric(THERAPY_START), linetype=THERAPY_LTY2, colour=THERAPY_COLOUR2),
						geom_vline(xintercept=as.numeric(therapy2), linetype=THERAPY_LTY, colour=THERAPY_COLOUR),
						p$layers
					)
				}
			}
			p$layers <- c(geom_abline(intercept=stats[, "Model.Intercept"], slope=stats[, "Model.Slope"], colour=REGRESS_COLOUR, linetype=REGRESS_LTY, size=regression.size), p$layers)
			if (!is.na(pat.id2))
				p$layers <- c(geom_abline(intercept=stats.2[, "Model.Intercept"], slope=stats.2[, "Model.Slope"], colour=REGRESS_COLOUR2, linetype=REGRESS_LTY2, size=regression.size), p$layers)
#					geom_path(aes(x=date, y=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(x=date, y=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		}
	} else
		p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
	p
}

apply.theme <- function(p, flipped=F, scaled=T, type.label.=type.label, type.value.=type.value, colour.break.=colour.break, colour.label.=colour.label, colour.value.=colour.value) {
	apply.axes(p + theme_bw() +
		theme(
			text=element_text(size=35),
			axis.text=element_text(size=30, colour='black'),
			legend.text=element_text(size=30),
			legend.position=0,
			axis.text.x=element_text(angle=60, hjust=1),
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
		scale_shape_manual(name="", breaks=type.label., labels=type.label., values=type.value., limits=type.label.) +
		scale_colour_manual(name="", breaks=colour.break., labels=colour.label., values=colour.value., limits=colour.break.) +
#		scale_size_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(1.67, 2), limits=c(0, 1)) +
		guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1)),
		flipped, scaled)
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--patid", type='character')
op <- add_option(op, "--patid2", type='character')
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--distmax", type='double', default=NA)
op <- add_option(op, "--distmin", type='double', default=NA)
op <- add_option(op, "--distby", type='double', default=NA)
op <- add_option(op, "--yearstart", type='character', default=NA)
op <- add_option(op, "--yearend", type='character', default=NA)
op <- add_option(op, "--yearby", type='double', default=NA)
op <- add_option(op, "--therapy", type='character', default=NA)
op <- add_option(op, "--therapy2", type='character', default=NA)
op <- add_option(op, "--therapyend", type='character', default=NA)
op <- add_option(op, "--liktol", type='numeric', default=1e-3)
op <- add_option(op, "--usedups", type='logical', action='store_true', default=F)
op <- add_option(op, "--cartoon", type='logical', action='store_true', default=F)
op <- add_option(op, "--xtitle", type='character', default="Years since first collection")
op <- add_option(op, "--histbymonth", type='logical', action='store_true', default=F)
op <- add_option(op, "--histheight", type='numeric', default=1.7)
op <- add_option(op, "--mincoltime", type='numeric', default=MIN_COL_TIME)
op <- add_option(op, "--maxcoltime", type='numeric', default=MAX_COL_TIME)
op <- add_option(op, "--histfreqby", type='numeric', default=2)
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
cartoon <- args$cartoon
x.title <- args$xtitle
hist.by.month <- args$histbymonth
hist.height <- args$histheight
MIN_COL_TIME <- args$mincoltime
MAX_COL_TIME <- args$maxcoltime
hist.freq.by <- args$histfreqby
if (use.real) {
	year.start <- args$yearstart
	year.end <- args$yearend
	THERAPY_START <- as.numeric(as.Date(args$therapy))
	
	if (!is.na(therapy2)) {
		therapy2 <- as.Date(therapy2)
	}
	
	if (!is.na(therapyend)) {
		therapyend <- as.Date(therapyend)
	}
} else {
	year.start <- as.numeric(args$yearstart)
	year.end <- as.numeric(args$yearend)
	THERAPY_START <- as.numeric(args$therapy)
}
if (is.na(therapyend)) {
	therapyend <- therapy2
	if (is.na(therapyend))  {
		therapyend <- Inf
	}
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
pdf.colour.disttree.file <- paste0("plots/", pat.id, ".colour.disttree.pdf")
pdf.colour.tree.file <- paste0("plots/", pat.id, ".colour.tree.pdf")
pdf.tree.file <- paste0("plots/", pat.id, ".tree.pdf")
pdf.hist.file <- paste0("plots/", pat.id, ".hist.pdf")
pdf.histdate.file <- paste0("plots/", pat.id, ".histdate.pdf")

tree <- ape::ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
data <- data[match(tree$tip.label, data$tip.label), ]
rownames(data) <- NULL
stats <- read.csv(stats.file, stringsAsFactors=F)
g <- readRDS(regression.file)

if (!is.na(pat.id2)) {
	data.2 <- read.csv(data.2.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
	data.2 <- data.2[match(tree$tip.label, data.2$tip.label), ]
	rownames(data.2) <- NULL
	stats.2 <- read.csv(stats.2.file, stringsAsFactors=F)
	g.2 <- readRDS(regression.2.file)
	
	data.2$LM <- 2
	data$LM <- 1
	data.2$censored[data.2$censored == 1] <- 2
	
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
	write.table(data, paste0("stats/", pat.id, ".comb.data.csv"), col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Estimated Date", "Date Difference", "LM"), row.names=F, sep=",")
	data <- data[, -8]
	
	stats$Minimum.Time.Point <- min(stats$Mimimum.Time.Point, stats.2$Minimum.Time.Point)
	stats$Minimum.Training.Time.Point <- min(stats$Minimum.Training.Time.Point, stats.2$Minimum.Training.Time.Point)
	stats$Maximum.Time.Point <- max(stats$Maximum.Time.Point, stats.2$Maximum.Time.Point)
	stats$Maximum.Training.Time.Point <- max(stats$Maximum.Training.Time.Point, stats.2$Maximum.Training.Time.Point)
	stats$Training.Samples <- stats$Training.Samples + stats.2$Training.Samples
	stats$Censored.RMSD <- sqrt(sum(data[data$censored > 0, "date.diff"]^2)/sum(data$censored > 0))
	stats$Censored.MAE <- sum(abs(data[data$censored > 0, "date.diff"]))/sum(data$censored > 0)
	stats <- as.data.frame(cbind(stats, stats.2[, c("AIC", "null.AIC", "p.value", "Model.Intercept", "Model.Slope", "Rsquared", "Training.RMSE", "Training.MAE")]))
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
	"Rsquared",
	"Estimated Root Date",
	"ERD CI low",
	"ERD CI high",
	"Model Fit",
	"Training RMSE",
	"Censored RMSD",
	"Total RMSD",
	"Training MAE",
	"Censored MAE",
	"Total MAE",
	"Total Concordance",
	"Bin Test (p)",
	"Bin Test (mean)",
	"T-Test (p)",
	"T-Test (mean)",
		"AIC 2",
		"null AIC 2",
		"p-value 2",
		"Model Intercept 2",
		"Model Slope 2",
		"Rsquared 2",
		"Training RMSE 2",
		"Training MAE 2"
	)
	write.table(stats, paste0("stats/", pat.id, ".comb.stats.csv"), row.names=F, col.names=stats.col.names, sep=",")

	colour.break <- c(0, -1, 1, 2)
	colour.label <- c("Training", "Training 2", "Censored", "Censored 2")
	colour.value <- c(TRAINING_COLOUR, TRAINING_COLOUR2, CENSORED_COLOUR, CENSORED_COLOUR2)
} else {
	mu <- stats[, "Model.Slope"]

	colour.break <- c(0, 1)
	colour.label <- c("Training", "Censored")
	colour.value <- c(TRAINING_COLOUR, CENSORED_COLOUR)
}

node.dates <- tryCatch(estimate.dates(tree, c(data$date, stats$Estimated.Root.Date, rep(NA, tree$Nnode - 1)), mu, node.mask=length(tree$tip.label) + 1, lik.tol=0, nsteps=1000, show.steps=100, opt.tol=1e-16), error=function(e) estimate.dates(tree, data$date, mu, lik.tol=0, nsteps=1000, show.steps=100, opt.tol=1e-16))

if (is.na(dist.min))
	dist.min <- -0.01
if (is.na(dist.max))
	dist.max <- max(data$dist) * 1.01
if (is.na(dist.by))
	dist.by <- 10^(floor(log10((dist.max - dist.min))))
if (use.real) {
	node.years <- as.numeric(as.character(as.Date(node.dates, origin="1970-01-01"), format="%Y"))

	if (is.na(year.by))
		year.by <- 2
	if (is.na(year.start))
		year.start <- floor(min(node.years) / year.by) * year.by
	if (is.na(year.end))
		year.end <- floor(max(node.years) / year.by + 1) * year.by
} else {
	if (is.na(year.by))
		year.by <- 365.25
	if (is.na(year.start))
		year.start <- floor(min(node.dates) / year.by) * year.by
	if (is.na(year.end))
		year.end <- floor(max(node.dates) / year.by + 1) * year.by
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
	if (cartoon) {
		date.ticks <- seq(floor(year.start / year.by) * year.by, year.end, by=year.by)
		x.scale <- scale_x_continuous(name="Years since first collection", breaks=date.ticks, limits=c(year.start, year.end))
		y.scale <- scale_y_continuous(name="Years since first collection", breaks=date.ticks, limits=c(year.start, year.end))
	} else {
		date.ticks <- seq(floor(year.start / year.by) * year.by, year.end, by=year.by)
		x.scale <- scale_x_continuous(name=x.title, breaks=date.ticks, labels=as.integer(date.ticks / 365.25), limits=c(year.start, year.end))
		y.scale <- scale_y_continuous(name=x.title, breaks=date.ticks, labels=as.integer(date.ticks / 365.25), limits=c(year.start, year.end))
	}
	
	type.break <- c("PLASMA", "PBMC")
	type.label <- c("RNA", "DNA")
	type.value <- c(16, 5)
}


x.div <- scale_x_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by))
y.div <- scale_y_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by))

data$type <- type.label[match(data$type, type.break)]

type.mask <- type.label %in% data$type
type.label <- unique(type.label[type.mask])
type.value <- unique(type.value[type.mask])

data$my.colour <- data$date
data$my.colour[data$censored > 0] <- "censored"
my.colour.break <- unique(data$my.colour)
my.colour.filter <- suppressWarnings(!is.na(as.numeric(my.colour.break)))
my.colour.value <- rep('black', length(my.colour.break))
my.colour.scale <- (as.numeric(my.colour.break[my.colour.filter]) - MIN_COL_TIME) / (MAX_COL_TIME - MIN_COL_TIME)
my.colour.value[my.colour.filter] <- hsv(my.colour.scale * .75, 0.5, 0.5)

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type="NODE", censored=NA, date=NA, dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], est.date=NA, date.diff=NA, my.colour=NA)), node.date=node.dates))

ptree <- phylo4d(tree, all.data=data.all)

if (F) {
pdf(pdf.file)
apply.theme(ggplot(data) +
	geom_point(aes(x=date, y=dist, colour=factor(censored), shape=type), size=6)) +
	theme(legend.position=c(0.02, .98))
dev.off()
}

if (cartoon) {
	point.size <- 8
	dist.tree.size <- 1
	tree.size <- 1
	scale.offset <- 0.3
	regression.size <- 2
} else {
	point.size <- 6
	dist.tree.size <- .6
	tree.size <- .4
	scale.offset <- 3
	regression.size <- 1
}	

#pdf(pdf.disttree.file)
#apply.theme(ggtree(ptree, colour="#49494980", size=dist.tree.size, ladderize=T) +
#	geom_tippoint(aes(colour=factor(censored), shape=type), size=point.size) +
#	geom_treescale(width=0.02, fontsize=7, offset=scale.offset),
#	flipped=F, scaled=F)
#dev.off()

pdf(pdf.colour.disttree.file)
apply.theme(ggtree(ptree, colour="#49494980", size=dist.tree.size, ladderize=T) +
	geom_tippoint(aes(colour=my.colour, shape=type), size=point.size) +
	geom_treescale(width=0.02, fontsize=7, offset=scale.offset),
	flipped=F, scaled=F, colour.value=my.colour.value, colour.label=my.colour.break, colour.break=my.colour.break)
dev.off()

#pdf(pdf.tree.file)
#apply.theme(ggtree(ptree, yscale="node.date", colour="#6d4d4180", size=tree.size, ladderize=F) +
#	geom_tippoint(aes(colour=factor(censored), shape=type), size=point.size),
#	flipped=T)
#dev.off()

pdf(pdf.colour.tree.file)
apply.theme(ggtree(ptree, yscale="node.date", colour="#49494980", size=tree.size, ladderize=F) +
	geom_tippoint(aes(colour=my.colour, shape=type), size=point.size),
	flipped=T, scaled=T, colour.value=my.colour.value, colour.label=my.colour.break, colour.break=my.colour.break)
dev.off()

data.hist <- subset(data, censored > 0)

date.levels <- sort(unique(data.hist$date))

if (use.real) {
	if (hist.by.month) {
		m.month <- as.numeric(as.character(as.Date(min(c(data.hist$est.date, data$date)), origin="1970-01-01"), "%m"))
		m.year <- as.numeric(as.character(as.Date(min(c(data.hist$est.date, data$date)), origin="1970-01-01"), "%Y"))
		M.month <- as.numeric(as.character(as.Date(max(c(data.hist$est.date, subset(data, censored <= 0)$date)), origin="1970-01-01"), "%m")) + 1
		M.year <- as.numeric(as.character(as.Date(max(c(data.hist$est.date, subset(data, censored <= 0)$date)), origin="1970-01-01"), "%Y"))
		
		if (M.month == 13) {
			M.year <- M.year + 1
			M.month <- 1
		}
		
		breaks <- if (M.year == m.year) {
			as.numeric(as.Date(paste0(m.year, "-", seq(m.month, M.month), "-01")))
		} else if (M.year == m.year + 1) {
			as.numeric(as.Date(c(paste0(m.year, "-", seq(m.month, 12), "-01"), paste0(M.year, "-", seq(1, M.month), "-01"))))
		} else {
			as.numeric(as.Date(c(paste0(m.year, "-", seq(m.month, 12), "-01"), paste0(unlist(lapply((m.year + 1):(M.year - 1), rep, 12)), "-", 1:12, "01"), paste0(M.year, "-", seq(1, M.month), "-01"))))
		}
	}
	else {
		m <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(min(c(data.hist$est.date, data$date)), origin="1970-01-01")))
		M <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(max(c(data.hist$est.date), subset(data, censored <= 0)$date), origin="1970-01-01"))) + 1
		breaks <- as.numeric(as.Date(paste0(seq(m, M), "-01-01")))
	}
	m <- breaks[1]
	M <- breaks[length(breaks[1])]
	date.vals <-  if (length(date.levels) == 4) c("#999999", "#666666", "#333333", "#000000") else if (length(date.levels) == 2) c("#666666", "#000000") else 'black'
	date.labs <- as.character(as.Date(date.levels, origin="1970-01-01"), format="%b. %Y")
} else {
	m <- floor(min(c(data.hist$est.date, data$date)) / 365.25)
	M <- floor(max(c(data.hist$est.date, subset(data, censored <= 0)$date)) / 365.25) + 1
	breaks <- seq(m, M) * 365.25
	date.vals <- rep('black', length(date.levels))
	date.labs <- rep('LAB', length(date.levels))
}

H <- max(hist(data.hist$est.date, breaks=breaks)$counts)

if (use.real) {
	p <- ggplot(data.hist, aes(x=est.date, fill=factor(date, levels=date.levels)))
} else {
	p <- ggplot(data.hist, aes(x=est.date))
}

if (!is.na(THERAPY_START)) {
	if (is.na(therapy2)) {
		p <-  p + geom_segment(x=as.numeric(THERAPY_START), xend=as.numeric(THERAPY_START), y=-Inf, yend=H * 1.05, linetype=THERAPY_LTY, colour=THERAPY_COLOUR)
	} else {
		p <- p + geom_segment(x=as.numeric(THERAPY_START), xend=as.numeric(THERAPY_START), y=-Inf, yend=H * 1.05, linetype=THERAPY_LTY2, colour=THERAPY_COLOUR2) +
			geom_segment(x=as.numeric(therapy2), xend=as.numeric(therapy2), y=-Inf, yend=H * 1.05, linetype=THERAPY_LTY, colour=THERAPY_COLOUR)
	}
}

p <- p + geom_segment(x=min(data$date), xend=min(data$date), y=H * 9 / 8, yend=H, arrow=arrow(length = unit(0.25, "cm")))

if (use.real) {
	p <- p + geom_histogram(breaks=breaks) + scale_fill_manual(name="Collection Date", values=date.vals, labels=date.labs)
} else {
	p <- p + geom_histogram(breaks=breaks, fill='black')
}


if (use.real) {
	if (hist.by.month) {
		if (M.year == m.year)  {
			p <- p + scale_x_continuous(name="Estimated integration month", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%b"))
		} else {
			p + scale_x_continuous(name="Estimated integration month", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%b %Y"))
		}
	}
	else {
		p <- p + scale_x_continuous(name="Estimated integration year", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%Y"))
	}
} else {
	p <- p + scale_x_continuous(name=x.title, breaks=breaks, labels=seq(m, M))
}
	
p <- p +
	guides(fill=guide_legend(override.aes=list(size=8, colour="#00000000"), ncol=2)) +
	theme_bw() +
	theme(
		text=element_text(size=35),
		axis.text=element_text(size=30, colour='black'),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
		legend.text=element_text(size=30),
		legend.title=element_text(size=30),
		legend.position=c(.98, .98),
		legend.justification=c(1, 1),
		legend.spacing=unit(0, 'cm'),
		legend.margin=margin(0, 0, 0, 0, 'cm'),
		legend.key.size=unit(1.2, 'cm'),
		legend.background=element_blank(),
		legend.box.background=element_blank(),
		panel.grid.major=element_blank(),
	    panel.grid.minor=element_blank()
	)
	
if (!use.real || length(date.levels) == 1) {
	p <- p + scale_y_continuous(name="Frequency", breaks=seq(0, H, by=hist.freq.by), limits=c(0, H * hist.height)) +
		theme(legend.position='none')
} else {
	p <- p + scale_y_continuous(name="Frequency", breaks=seq(0, H, by=hist.freq.by), limits=c(0, H * hist.height))
}

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=5, colormodel='rgb')
pdf(pdf.hist.file)
p
dev.off()