#!/usr/bin/Rscript

library(ape)
library(ggtree)
library(phylobase)

TRAINING_COLOUR <- 'black'
TRAINING_COLOUR2 <- 'darkblue'
CENSORED_COLOUR <- 'red'
CENSORED_COLOUR2 <- "#ff8000"

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=10, colormodel='rgb')

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
				coord_flip(xlim=c(dist.min, dist.max)) +
				geom_abline(intercept=-stats[, "Model.Intercept"] / stats[, "Model.Slope"], slope=1 / stats[, "Model.Slope"], colour=REGRESS_COLOUR, linetype=2, size=regression.size)
#				geom_path(aes(y=date, x=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(y=date, x=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (!is.na(THERAPY_START)) {
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
				p <- p + geom_abline(intercept=-stats.2[, "Model.Intercept"] / stats.2[, "Model.Slope"], slope=1 / stats.2[, "Model.Slope"], colour=REGRESS_COLOUR2, linetype=2, size=regression.size)
#					geom_path(aes(y=date, x=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(y=date, x=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		} else {
			p <- p + 
				x.scale +
				y.div +
				coord_cartesian(ylim=c(dist.min, dist.max)) +
				geom_abline(intercept=stats[, "Model.Intercept"], slope=stats[, "Model.Slope"], colour=REGRESS_COLOUR, linetype=2, size=regression.size)
#				geom_path(aes(x=date, y=ci.high), data=data.ci, linetype=3, colour="#0060b080") +
#				geom_path(aes(x=date, y=ci.low), data=data.ci, linetype=3, colour="#0060b080")
			if (!is.na(THERAPY_START)) {
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
				p <- p + geom_abline(intercept=stats.2[, "Model.Intercept"], slope=stats.2[, "Model.Slope"], colour=REGRESS_COLOUR2, linetype=2, size=regression.size)
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
		scale_shape_manual(name="", breaks=type.label, labels=type.label, values=type.value, limits=type.label) +
		scale_colour_manual(name="", breaks=colour.break, labels=colour.label, values=colour.value, limits=colour.break) +
#		scale_size_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(1.67, 2), limits=c(0, 1)) +
		guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1)),
		flipped, scaled)
}

args <- commandArgs(T)

tree.file <- args[1]
info.file <- args[2]
plot.file <- args[3]
dot.spread <- as.numeric(args[4]) # 0.0045 for pat 1, 0.0012 for pat 2

tree <- read.tree(tree.file)
info <- read.csv(info.file, stringsAsFactors=F)

info <- subset(info, KEPTDUP == 1)
info$COLDATE <- as.numeric(as.Date(info$COLDATE))

ftree <- fortify(tree, colour="#49494980", size=dist.tree.size, ladderize=T)

tip.points <- apply(ftree[1:length(tree$tip.label), ], 1, function (x) {
	dups <- subset(info, DUPLICATE == x[["label"]])
	dups <- dups[order(with(dups, COLDATE + max(COLDATE) * as.numeric(DUPLICATE != FULLSEQID))), ]
	dups$x <- as.numeric(x[["x"]]) + (0:(nrow(dups) - 1)) * dot.spread
	dups$y <- as.numeric(x[["y"]])
		
	bars <- dups[which(dups$COLDATE[1:(nrow(dups) - 1)] != dups$COLDATE[2:nrow(dups)]), ]
	if (nrow(bars) > 0) {
		bars$x <- bars$x + dot.spread / 2
		bars$TYPE <- 'BAR'
		bars$CENSORED <- -100
		bars$FULLSEQID <- 'BAR'
		bars$DUPLICATE <- 'BAR'
				
		rbind(dups, bars)
	} else
		dups
})

if (-1 %in% info$CENSORED) {
	clade <- get.child.edges(tree, getMRCA(tree, which(unlist(lapply(tip.points, function(x) x[1, 'CENSORED'] == -1)))))
	clade.tips <- tree$edge[clade, 2]
	clade.tips <- clade.tips[clade.tips <= length(tree$tip.label)]
	
	tip.points[clade.tips] <- lapply(tip.points[clade.tips], function(x) {x$CENSORED[x$CENSORED == 1] <- 2; x})
}

tip.points <- do.call(rbind, tip.points)

type.break <- c("PLASMA", "PBMC", "PBMC (cultured)", "WHOLE BLOOD", "PBMC (REACTIVE)", "PBMC (CULTURE)", "BAR")
type.label <- c("RNA", "DNA", "DNA", "DNA", "RNA", "DNA", "BAR")
type.value <- c(16, 5, 5, 5, 16, 5, 124)
point.size <- 2.5
dist.tree.size <- .6
tree.size <- .4
scale.offset <- 3
regression.size <- 1
colour.break <- c(0, -1, 1, 2, -100)
colour.label <- c("Training", "Training 2", "Censored", "Censored 2", "BAR")
colour.value <- c(TRAINING_COLOUR, TRAINING_COLOUR2, CENSORED_COLOUR, CENSORED_COLOUR2, 'black')

tip.points$TYPE <- type.label[match(tip.points$TYPE, type.break)]

type.mask <- type.label %in% tip.points$TYPE
type.label <- unique(type.label[type.mask])
type.value <- unique(type.value[type.mask])

colour.mask <- colour.break %in% tip.points$CENSORED
colour.label <- colour.label[colour.mask]
colour.value <- colour.value[colour.mask]
colour.break <- colour.break[colour.mask]

pdf(plot.file)
apply.theme(ggtree(ftree, colour="#49494980", size=dist.tree.size) +
	geom_point(aes(x=x, y=y, colour=factor(CENSORED), shape=TYPE, alpha=DUPLICATE == FULLSEQID), data=tip.points, size=point.size) +
	geom_treescale(width=0.02, fontsize=7, offset=scale.offset) +
	scale_alpha_manual(name="", breaks=c(F, T), limits=c(F, T), values=c(0.5, 1)), flipped=F, scaled=F)
dev.off()
