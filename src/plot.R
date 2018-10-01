#!/usr/bin/Rscript

library(ape)
library(phylobase)
library(ggplot2)
library(ggtree)
library(optparse)
library(chemCal)

source("/opt/node.dating.R")

#MIN_COL_TIME <- 9720
#MAX_COL_TIME <- 13371


get.val <- function(x, default) if (is.null(x)) default else x

THERAPY_COLOUR <- "#ffffa0"
THERAPY_LTY <- 3
THERAPY_COLOUR2 <- "#ffcccc"
THERAPY_LTY2 <- 6
THERAPY_COLOUR3 <- "#ffffa0"
THERAPY_LTY3 <- 3

TRAINING_COLOUR <- 'black'
TRAINING_COLOUR2 <- 'darkblue'
CENSORED_COLOUR <- 'red'
CENSORED_COLOUR2 <- "#ff8000"

REGRESS_COLOUR <- "#606060b0"
REGRESS_LTY <- 2
REGRESS_COLOUR2 <- "#aaaaaab0"
REGRESS_LTY2 <- 4

colour.blind <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999")

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb', useDingbats=F)

#print.plot <- function(..., width=7, height=7) pdf(..., width=width, height=height)
print.plot <- function(file.name, p, width=7, height=7) ggsave(gsub("pdf", "png", file.name), p, device='png', width=width, height=height)

make.mean.row <- function(label, data) {
	data.dup <- data[data$duplicate == label, ]
	
	if (any(data.dup$censored <= 0))
		data.dup <- subset(data.dup, censored <= 0)
		
	data.dup$tip.label[1] <- label
	data.dup$date[1] <- weighted.mean(data.dup$date, data.dup$weight)
	
	data.dup[1, ]
}

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

add.therapy <- function(therapy.start, therapy.end, therapy.type, flipped=F, height=Inf) {
	if (T) {
		therapy.colour <- if (therapy.type == 1) THERAPY_COLOUR else THERAPY_COLOUR2	
	
		if (flipped) {
			geom_rect(ymin=therapy.start, ymax=therapy.end, xmin=-Inf, xmax=height, fill=therapy.colour, linetype=2)
		} else {
			geom_rect(xmin=therapy.start, xmax=therapy.end, ymin=-Inf, ymax=height, fill=therapy.colour, linetype=2)
		}
	} else {
		therapy.line <- if (therapy.type == 1) THERAPY_LTY else THERAPY_LTY2
		
		if (flipped) {
			geom_hline(yintercept=therapy.start, linetype=therapy.line)
		} else {
			geom_vline(xintercept=therapy.start, linetype=therapy.line)
		}
	}
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
			p$layers <- c(geom_abline(intercept=stats[, "Model.Intercept"], slope=stats[, "Model.Slope"], colour=REGRESS_COLOUR, linetype=REGRESS_LTY, size=regression.size), p$layers)
			if (!is.na(pat.id2))
				p$layers <- c(geom_abline(intercept=stats.2[, "Model.Intercept"], slope=stats.2[, "Model.Slope"], colour=REGRESS_COLOUR2, linetype=REGRESS_LTY2, size=regression.size), p$layers)
#					geom_path(aes(x=date, y=ci.high.2), data=data.ci, linetype=3, colour="#00305880") +
#					geom_path(aes(x=date, y=ci.low.2), data=data.ci, linetype=3, colour="#00305880")
		}
		
		if (!is.na(THERAPY_START)) {
			if (is.na(therapy2))
				p$layers <- c(add.therapy(as.numeric(THERAPY_START), Inf, 1, flipped), p$layers)
			else {
				p$layers <- c(
					add.therapy(as.numeric(THERAPY_START), Inf, 1, flipped),
					add.therapy(as.numeric(therapy2), as.numeric(THERAPY_START), 2, flipped),
					p$layers
				)
			}
			
			if (!is.na(therapy3))
				p$layers <- c(add.therapy(as.numeric(therapy3), as.numeric(therapy3end), 1, flipped), p$layers)
		}
	} else
		p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
	p
}

apply.theme <- function(p, flipped=F, scaled=T, type.break.=type.break, type.value.=type.value, colour.break.=colour.break, colour.value.=colour.value) {
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
		scale_shape_manual(name="", breaks=type.break., values=type.value., limits=type.break.) +
		scale_colour_manual(name="", breaks=colour.break., values=colour.value., limits=colour.break.) +
		guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1)),
		flipped, scaled)
}

plot.hist <- function(data, pdf.hist.file) {
	data.hist <- subset(data, censored > 0)
		
	if (use.rainbow) {
		date.levels <- sort(unique(data.hist$date))
	} else {
		type.filter <- my.colour.break %in% data.hist$type
		type.levels <- my.colour.break[type.filter]
		colour.vals <- my.colour.value[type.filter]
		colour.labs <- type.levels
	}
	
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
				as.numeric(as.Date(c(paste0(m.year, "-", seq(m.month, 12), "-01"), paste0(unlist(lapply((m.year + 1):(M.year - 1), rep, 12)), "-", 1:12, "-01"), paste0(M.year, "-", seq(1, M.month), "-01"))))
			}
		}
		else {
			m <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(min(c(data.hist$est.date, data$date)), origin="1970-01-01")))
			M <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(max(c(data.hist$est.date), subset(data, censored <= 0)$date), origin="1970-01-01"))) + 1
			breaks <- as.numeric(as.Date(paste0(seq(m, M), "-01-01")))
		}
		m <- breaks[1]
		M <- breaks[length(breaks[1])]
		if (use.rainbow) {
			colour.vals <- if (length(date.levels) == 4) c("#cccccc", "#888888", "#444444", "#000000") else if (length(date.levels) == 3) c("#aaaaaa", "#555555", "#000000") else if (length(date.levels) == 2) c("#999999", "#000000") else 'black'
			colour.labs <- as.character(as.Date(date.levels, origin="1970-01-01"), format="%b. %Y")
		}
	} else {
		if (cartoon) {
			m <- floor(min(c(data.hist$est.date, data$date)))
			M <- floor(max(c(data.hist$est.date, subset(data, censored <= 0)$date))) + 1
			breaks <- seq(m, M)
			if (use.rainbow) {
				colour.vals <- rep('black', length(date.levels))
				colour.labs <- rep('LAB', length(date.levels))
			}
		} else {
			m <- floor(min(c(data.hist$est.date, data$date)) / 365.25)
			M <- floor(max(c(data.hist$est.date, subset(data, censored <= 0)$date)) / 365.25) + 1
			breaks <- seq(m, M) * 365.25
			if (use.rainbow) {
				colour.vals <- rep('black', length(date.levels))
				colour.labs <- rep('LAB', length(date.levels))
			}
		} 
	}
	
	H <- max(hist(data.hist$est.date, breaks=breaks, plot=F)$counts)
	
	if (use.rainbow) {
		if (use.real) {
			p <- ggplot(data.hist, aes(x=est.date, fill=factor(date, levels=date.levels)))
		} else {
			p <- ggplot(data.hist, aes(x=est.date))
		}
	} else {
		p <- ggplot(data.hist, aes(x=est.date, fill=factor(type, levels=type.levels)))
	}
	
	if (!is.na(THERAPY_START)) {
		p <-  p + add.therapy(as.numeric(THERAPY_START), Inf, 1, height=H * 1.05)
		if (!is.na(therapy2))
			p <-  p + add.therapy(as.numeric(therapy2), as.numeric(THERAPY_START), 2, height=H * 1.05)
	}
	
	if (!is.na(therapy3))
		p <- p + add.therapy(as.numeric(therapy3), as.numeric(therapy3end), 1, height=H * 1.05)
	
	p <- p + geom_segment(x=min(data$date), xend=min(data$date), y=H * 9 / 8, yend=H, arrow=arrow(length = unit(0.25, "cm")))
	
	if (use.rainbow) {
		if (use.real) {
			p <- p + geom_histogram(breaks=breaks) + scale_fill_manual(name="Collection Date", values=colour.vals, labels=colour.labs)
		} else {
			p <- p + geom_histogram(breaks=breaks, fill='black')
		}
	} else {
		p <- p + geom_histogram(breaks=breaks) + scale_fill_manual(name="Type", values=colour.vals, labels=colour.labs)
	}
	
	
	if (use.real) {
		if (hist.by.month) {
			if (M.year == m.year)  {
				p <- p + scale_x_continuous(name="Estimated Integration Month", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%b"))
			} else {
				p <- p + scale_x_continuous(name="Estimated Integration Month", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%b"))
			}
		}
		else {
			p <- p + scale_x_continuous(name="Estimated Integration Year", breaks=breaks, labels=as.character(as.Date(breaks, origin="1970-01-01"), "%Y"))
		}
	} else {
		p <- p + scale_x_continuous(name="Estimated Integration Year", breaks=breaks, labels=seq(m, M))
	}
	
	p <- p +
		guides(fill=guide_legend(override.aes=list(size=8, colour="#00000000"), ncol=2)) +
		theme_bw() +
		theme(
			text=element_text(size=30),
			axis.text=element_text(size=25, colour='black'),
			axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
			legend.text=element_text(size=25),
			legend.title=element_text(size=25),
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
	
	p <- p + scale_y_continuous(name="Frequency", breaks=seq(0, H, by=hist.dup.freq.by), limits=c(0, H * hist.height))
	
	if (use.rainbow && !(use.real && length(date.levels) > 1)) {
		p <- p + theme(legend.position='none')
	}
	
	ggsave(pdf.hist.file, p, height=5)
}

plot.dupes <- function(info, pdf.dup.disttree.file, pdf.dup.tree.file, pdf.dup.hist.file, disttree.dupalpha=0.4, tree.dupalpha=0.1) {
	info$COLDATE <- if (use.real) as.numeric(as.Date(info$COLDATE, origin="1970-01-01")) else info$COLDATE
	info <- info[order(info$TYPE), ]
	info <- info[order(info$COLDATE), ]
	info <- info[order(info$DUPLICATE), ]
	
	fort <- fortify(ptree, colour="#49494980", size=dist.tree.size, ladderize=T)
	fort.dup <- fort[match(info$DUPLICATE, fort$tip.label), ]
	
	fort.dup$tip.label <- info$FULLSEQID
	fort.dup$duplicate <- info$DUPLICATE
	fort.dup$censored <- info$CENSORED
	fort.dup$dup.censored <- info[match(info$DUPLICATE, info$FULLSEQID), "CENSORED"]
	fort.dup$type <- info$TYPE
	fort.dup$my.colour <- NA
	fort.dup$date <- info$COLDATE
	fort.dup$x.shift <- fort.dup$x
	fort.dup$my.type <- with(fort.dup, toupper(paste0(type, censored > 0)))
	
	if (use.rainbow) {
		fort.dup$my.colour <- fort.dup$date
		fort.dup$my.colour[fort.dup$censored > 0] <- "censored"
		my.colour.break <- unique(fort.dup$my.colour)
		my.colour.filter <- suppressWarnings(!is.na(as.numeric(my.colour.break)))
		my.colour.value <- rep('black', length(my.colour.break))
		my.colour.scale <- (as.numeric(my.colour.break[my.colour.filter]) - MIN_COL_TIME) / (MAX_COL_TIME - MIN_COL_TIME)
		my.colour.value[my.colour.filter] <- hsv(my.colour.scale * .75, 0.5, 0.5)
	} else {
		fort.dup$my.colour <- toupper(fort.dup$type)
	}
	
	fort.dup$censored <- as.numeric(fort.dup$censored <= 0) * fort.dup$censored + 
		as.numeric(fort.dup$censored > 0 & fort.dup$dup.censored <= 0) * (1 - fort.dup$dup.censored) + 
		as.numeric(fort.dup$dup.censored > 0 & fort.dup$censored > 0) * (fort.dup$dup.censored)
	
	data.comb <- fort.dup
	fort.dup <- subset(fort.dup, tip.label != duplicate)
	
	for (i in 1:nrow(fort.dup)) {
		if (i == 1 || fort.dup$duplicate[i - 1] != fort.dup$duplicate[i]) {
			shift <- 1
		} else {
			shift <- shift + 1
		}
		
		fort.dup$x.shift[i] <- fort.dup$x[i] + shift * dup.shift
	}
	
	print.plot(pdf.dup.disttree.file, apply.theme(ggtree(ptree, colour="#49494980", size=dist.tree.size, ladderize=T) +
					  	geom_tippoint(aes(colour=my.colour, shape=my.type), size=point.size) +
					  	geom_treescale(width=0.02, fontsize=7, offset=scale.offset) + 
					  	geom_point(aes(colour=my.colour, shape=my.type, x=x.shift, y=y), data=fort.dup, alpha=disttree.dupalpha, size=point.size),
					  flipped=F, scaled=F, colour.value.=my.colour.value, colour.break.=my.colour.break), height=10)

	data.mean <- do.call(rbind, lapply(tree$tip.label, make.mean.row, data.comb))
	dup.node.dates <- tryCatch(estimate.dates(tree, c(data.mean$date, stats$Estimated.Root.Date, rep(NA, tree$Nnode - 1)), mu, node.mask=length(tree$tip.label) + 1, lik.tol=0, nsteps=nsteps, show.steps=100, opt.tol=1e-16), error=function(e) estimate.dates(tree, data.mean$date, mu, lik.tol=0, nsteps=nsteps, show.steps=100, opt.tol=1e-16))
	mean.ptree <- phylo4d(tree, all.data=as.data.frame(cbind(rbind(data.mean[, c("tip.label", "type", "censored", "date", "dist", "weight", "date.diff", "est.date", "my.colour", "my.type")], data.frame(tip.label=paste0("N.", 1:tree$Nnode), type="NODE", censored=NA, date=NA, dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], weight=NA, est.date=NA, date.diff=NA, my.colour=NA, my.type=NA)), node.date=dup.node.dates)))
	
	print.plot(pdf.dup.tree.file, apply.theme(ggtree(mean.ptree, yscale="node.date", colour="#49494980", size=tree.size, ladderize=F) +
					  	geom_point(aes(colour=my.colour, shape=my.type, x=x, y=date), data=data.comb, alpha=tree.dupalpha, size=point.size),
					  flipped=T, scaled=T, colour.value.=my.colour.value, colour.break.=my.colour.break))
	
	plot.hist(data.comb, pdf.dup.hist.file)
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--patid", type='character')
op <- add_option(op, "--patid2", type='character')
op <- add_option(op, "--real", type='logical', action='store_true')
op <- add_option(op, "--distmax", type='double')
op <- add_option(op, "--distmin", type='double')
op <- add_option(op, "--distby", type='double')
op <- add_option(op, "--yearstart", type='character')
op <- add_option(op, "--yearend", type='character')
op <- add_option(op, "--yearby", type='double')
op <- add_option(op, "--therapy", type='character')
op <- add_option(op, "--therapy2", type='character')
op <- add_option(op, "--therapy3", type='character')
op <- add_option(op, "--therapy3end", type='character')
op <- add_option(op, "--liktol", type='numeric')
op <- add_option(op, "--plotdups", type='logical', action='store_true')
op <- add_option(op, "--cartoon", type='logical', action='store_true')
op <- add_option(op, "--xtitle", type='character')
op <- add_option(op, "--histbymonth", type='logical', action='store_true')
op <- add_option(op, "--histheight", type='numeric')
op <- add_option(op, "--mincoltime", type='numeric')
op <- add_option(op, "--maxcoltime", type='numeric')
op <- add_option(op, "--histfreqby", type='numeric')
op <- add_option(op, "--dnashapescale", type='numeric')
op <- add_option(op, "--maxvl", type='numeric')
op <- add_option(op, "--vlfile", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--dupshift", type='numeric')
op <- add_option(op, "--nsteps", type='numeric')
op <- add_option(op, "--marklatent", type='logical', action='store_true')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--stats", type='character')
op <- add_option(op, "--regression", type='character')
op <- add_option(op, "--outputfolder", type='character')
op <- add_option(op, "--types", type='character')
op <- add_option(op, "--typevalues", type='character')
op <- add_option(op, "--rainbow", type='logical', action='store_true')
op <- add_option(op, "--black", type='logical', action='store_false', dest="rainbow")
op <- add_option(op, "--colourvalues", type='character')
op <- add_option(op, "--histdupfreqby", type='numeric')
op <- add_option(op, "--data2", type='character')
op <- add_option(op, "--stats2", type='character')
op <- add_option(op, "--regression2", type='character')
op <- add_option(op, "--datacomb", type='character')
op <- add_option(op, "--statscomb", type='character')
op <- add_option(op, "--plotgroups", type='logical', action='store_true')
op <- add_option(op, "--latentedges", type='logical', action='store_true')
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
pat.id <- get.val(args$patid, NA)
pat.id2 <- get.val(args$patid2, NA)
use.real <- get.val(args$real, F) # 0 Training/Censored, 1 RNA/DNA
dist.min <- get.val(args$distmin, NA)
dist.max <- get.val(args$distmax, NA)
dist.by <- get.val(args$distby, NA)
year.by <- get.val(args$yearby, NA)
LIK_TOL <- get.val(args$liktol, 1e-3)
therapy2 <- get.val(args$therapy2, NA)
therapy3 <- get.val(args$therapy3, NA)
therapy3end <- get.val(args$therapy3end, NA)
cartoon <- get.val(args$cartoon, F)
x.title <- get.val(args$xtitle, "Collection Year")
hist.by.month <- get.val(args$histbymonth, F)
hist.height <- get.val(args$histheight, 1.7)
MIN_COL_TIME <- get.val(args$mincoltime, NA)
MAX_COL_TIME <- get.val(args$maxcoltime, NA)
hist.freq.by <- get.val(args$histfreqby, 2)
dna.shape.scale <- get.val(args$dnashapescale, 1)
max.vl <- get.val(args$maxvl, NA)
vl.file <- get.val(args$vlfile, NA)
use.dups <- get.val(args$plotdups, F)
info.file <- get.val(args$info, NA)
dup.shift <- get.val(args$dupshift, 0.01)
nsteps <- get.val(args$nsteps, 1000)
mark.latent <- get.val(args$marklatent, F)
data.file <- get.val(args$data, NA)
stats.file <- get.val(args$stats, NA)
regression.file <- get.val(args$regression, NA)
output.folder <- get.val(args$outputfolder, "plots")
types <- get.val(args$types, "PLASMA,PBMC")
type.values <- get.val(args$typevalues, "16,1,18,5")
use.rainbow <- get.val(args$rainbow, T)
colour.values <- get.val(args$colourvalues, NA)
hist.dup.freq.by <- get.val(args$histdupfreqby, 10)
data.2.file <- get.val(args$data2, NA)
stats.2.file <- get.val(args$stats2, NA)
regression.2.file <- get.val(args$regression2, NA)
comb.data.file <- get.val(args$datacomb, NA)
comb.stats.file <- get.val(args$statscomb, NA)
plot.groups <- get.val(args$plotgroups, F)
latent.edges <- get.val(args$lantentedges, F)
if (use.real) {
	year.start <- get.val(args$yearstart, NA)
	year.end <- get.val(args$yearend, NA)
	THERAPY_START <- as.numeric(as.Date(get.val(args$therapy, NA)))
	
	if (!is.na(therapy2)) {
		therapy2 <- as.Date(therapy2)
	}
	
	if (!is.na(therapy3)) {
		therapy3 <- as.Date(therapy3)
	}
	
	if (!is.na(therapy3end)) {
		therapy3end <- as.Date(therapy3end)
	}
} else {
	year.start <- as.numeric(get.val(args$yearstart, NA))
	year.end <- as.numeric(get.val(args$yearend, NA))
	THERAPY_START <- as.numeric(get.val(args$therapy, NA))
}

if (is.na(data.file)) data.file <- paste0("stats/", pat.id, ".data.csv")
if (is.na(stats.file)) stats.file <- paste0("stats/", pat.id, ".stats.csv")
if (is.na(regression.file)) regression.file <- paste0("stats/", pat.id, ".regression.rds")
if (!is.na(pat.id2)) {
	if (is.na(data.2.file)) data.2.file <- paste0("stats/", pat.id2, ".data.csv")
	if (is.na(stats.2.file)) stats.2.file <- paste0("stats/", pat.id2, ".stats.csv")
	if (is.na(regression.2.file)) regression.2.file <- paste0("stats/", pat.id2, ".regression.rds")
	if (is.na(comb.data.file)) comb.data.file <- paste0("stats/", pat.id, ".comb.data.csv")
	if (is.na(comb.stats.file)) comb.stats.file <- paste0("stats/", pat.id, ".comb.stats.csv")
}
pdf.file <- paste0(output.folder, "/", pat.id, ".pdf")
pdf.disttree.file <- paste0(output.folder, "/", pat.id, ".disttree.pdf")
pdf.colour.disttree.file <- paste0(output.folder, "/", pat.id, ".colour.disttree.pdf")
pdf.colour.tree.file <- paste0(output.folder, "/", pat.id, ".colour.tree.pdf")
pdf.tree.file <- paste0(output.folder, "/", pat.id, ".tree.pdf")
pdf.hist.file <- paste0(output.folder, "/", pat.id, ".hist.pdf")
pdf.histdate.file <- paste0(output.folder, "/",  pat.id, ".histdate.pdf")
pdf.vl.file <- paste0(output.folder, "/", pat.id, ".vl.pdf")
pdf.dup.disttree.file <- paste0(output.folder, "/", pat.id, ".dup.disttree.pdf")
pdf.dup.tree.file <- paste0(output.folder, "/", pat.id, ".dup.tree.pdf")
pdf.dup.hist.file <- paste0(output.folder, "/", pat.id, ".dup.hist.pdf")
pdf.colour.mark.tree.file <- paste0(output.folder, "/", pat.id, ".colour.mark.tree.pdf")
pdf.group.disttree.file <- paste0(output.folder, "/", pat.id, ".group.disttree.pdf")
pdf.group.tree.file <- paste0(output.folder, "/", pat.id, ".group.tree.pdf")
pdf.group.hist.file <- paste0(output.folder, "/", pat.id, ".group.hist.pdf")

tree <- ape::ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "weight", "est.date", "date.diff"), stringsAsFactors=F)
data <- data[match(tree$tip.label, data$tip.label), ]
rownames(data) <- NULL
stats <- read.csv(stats.file, stringsAsFactors=F)
g <- readRDS(regression.file)

if (!is.na(pat.id2)) {
	library(phangorn)
	
	data.2 <- read.csv(data.2.file, col.names=c("tip.label", "type", "censored", "date", "dist", "weight", "est.date", "date.diff"), stringsAsFactors=F)
	data.2 <- data.2[match(tree$tip.label, data.2$tip.label), ]
	rownames(data.2) <- NULL
	stats.2 <- read.csv(stats.2.file, stringsAsFactors=F)
	g.2 <- readRDS(regression.2.file)
	
#	data.2$LM <- 2
#	data$LM <- 1
	data.2$censored[data.2$censored == 1] <- 2
	
	mu <- rep(stats.2[, "Model.Slope"], nrow(tree$edge))
	tr <- function(x) -3 * x ^ 2 / 2 + x / 2 + 1
	anc <- ace(tr(data$censored), tree)
	first.edges <- anc$ace[tree$edge[, 1] - length(tree$tip.label) + 1] > 0
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
	write.table(data, comb.data.file, col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Weight", "Estimated Date", "Date Difference"), row.names=F, sep=",")

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
	write.table(stats, comb.stats.file, row.names=F, col.names=stats.col.names, sep=",")

	colour.break <- c(0, -1, 1, 2)
	colour.value <- c(TRAINING_COLOUR, TRAINING_COLOUR2, CENSORED_COLOUR, CENSORED_COLOUR2)
} else {
	mu <- rep(stats[, "Model.Slope"], nrow(tree$edge))

	colour.break <- c(0, 1)
	colour.value <- c(TRAINING_COLOUR, CENSORED_COLOUR)
}

# shift latent coalescing left
if (latent.edges) {
	tr <- function(x) (as.numeric(x <= 0) - 1/3) * 3
	anc <- ace(tr(data$censored), tree)
	first.edges <- anc$ace[tree$edge[, 1] - length(tree$tip.label) + 1] < 0
	mu[first.edges] <- 1e-8
}

node.dates <- tryCatch(estimate.dates(tree, c(data$date, stats$Estimated.Root.Date, rep(NA, tree$Nnode - 1)), mu, node.mask=length(tree$tip.label) + 1, lik.tol=0, nsteps=nsteps, show.steps=100, opt.tol=1e-16), error=function(e) estimate.dates(tree, data$date, mu, lik.tol=0, nsteps=nsteps, show.steps=100, opt.tol=1e-16))

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
} else {
	if (cartoon) {
		date.ticks <- seq(floor(year.start / year.by) * year.by, year.end, by=year.by)
		x.scale <- scale_x_continuous(name="Collection Year", breaks=date.ticks, limits=c(year.start, year.end))
		y.scale <- scale_y_continuous(name="Collection Year", breaks=date.ticks, limits=c(year.start, year.end))
	} else {
		date.ticks <- seq(floor(year.start / year.by) * year.by, year.end, by=year.by)
		x.scale <- scale_x_continuous(name=x.title, breaks=date.ticks, labels=as.integer(date.ticks / 365.25), limits=c(year.start, year.end))
		y.scale <- scale_y_continuous(name=x.title, breaks=date.ticks, labels=as.integer(date.ticks / 365.25), limits=c(year.start, year.end))
	}
}

type.break <- unlist(lapply(strsplit(types, split=",")[[1]], paste0, c("FALSE", "TRUE")))
type.value <- as.numeric(strsplit(type.values, split=",")[[1]])

x.div <- scale_x_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by))
y.div <- scale_y_continuous(name="Divergence from root", breaks=seq(0.0, dist.max, by=dist.by))

data$my.type <- with(data, paste0(toupper(type), censored > 0))

if (use.rainbow) {
	if (is.na(MIN_COL_TIME))
		MIN_COL_TIME <- min(subset(data, censored <= 0)$date)
	if (is.na(MAX_COL_TIME))
		MAX_COL_TIME <- max(subset(data, censored <= 0)$date)

	data$my.colour <- data$date
	data$my.colour[data$censored > 0] <- "censored"
	my.colour.break <- unique(data$my.colour)
	my.colour.filter <- suppressWarnings(!is.na(as.numeric(my.colour.break)))
	my.colour.value <- rep('black', length(my.colour.break))
	my.colour.scale <- (as.numeric(my.colour.break[my.colour.filter]) - MIN_COL_TIME) / (MAX_COL_TIME - MIN_COL_TIME)
	my.colour.value[my.colour.filter] <- hsv(my.colour.scale * .75, 0.5, 0.5)
} else {
	data$my.colour <- toupper(data$type)
	my.colour.break <- strsplit(types, split=",")[[1]]
	my.colour.value <- if (is.na(colour.values)) colour.blind[1:length(my.colour.break)] else strsplit(colour.values, split=",")[[1]]
}

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type="NODE", censored=NA, date=NA, dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], weight=NA, est.date=NA, date.diff=NA, my.colour=NA, my.type=NA)), node.date=node.dates))

ptree <- phylo4d(tree, all.data=data.all)

if (F) {
print.plot(pdf.file, apply.theme(ggplot(data) +
	geom_point(aes(x=date, y=dist, colour=factor(censored), shape=type), size=6)) +
	theme(legend.position=c(0.02, .98)))
}

if (cartoon) {
	point.size <- 8
	dist.tree.size <- 1
	tree.size <- 1
	scale.offset <- 0.3
	regression.size <- 2
	vl.linesize <- 1
} else {
	point.size <- 6
	dist.tree.size <- .6
	tree.size <- .4
	scale.offset <- 3
	regression.size <- 1
	vl.linesize <- 0.5
}	

#size.value <- point.size * c(1, 1, 1, dna.shape.scale)
#size.value <- size.value[type.mask]

#pdf(pdf.disttree.file)
#apply.theme(ggtree(ptree, colour="#49494980", size=dist.tree.size, ladderize=T) +
#	geom_tippoint(aes(colour=factor(censored), shape=type), size=point.size) +
#	geom_treescale(width=0.02, fontsize=7, offset=scale.offset),
#	flipped=F, scaled=F)
#dev.off()

print.plot(pdf.colour.disttree.file, apply.theme(ggtree(ptree, colour="#49494980", size=dist.tree.size, ladderize=T) +
	geom_tippoint(aes(colour=my.colour, shape=my.type), size=point.size) +
	geom_treescale(width=0.02, fontsize=7, offset=scale.offset),
	flipped=F, scaled=F, colour.value.=my.colour.value, colour.break.=my.colour.break))

#pdf(pdf.tree.file)
#apply.theme(ggtree(ptree, yscale="node.date", colour="#6d4d4180", size=tree.size, ladderize=F) +
#	geom_tippoint(aes(colour=factor(censored), shape=type), size=point.size),
#	flipped=T)
#dev.off()

print.plot(pdf.colour.tree.file, apply.theme(ggtree(ptree, yscale="node.date", colour="#49494980", size=tree.size, ladderize=F) +
	geom_tippoint(aes(colour=my.colour, shape=my.type), size=point.size),
	flipped=T, scaled=T, colour.value.=my.colour.value, colour.break.=my.colour.break))

plot.hist(data, pdf.hist.file)

if (!is.na(vl.file)) {
	data.vl <- read.csv(vl.file, stringsAsFactors=F)
	if (use.real)
		data.vl$Date <- as.Date(data.vl$Date)		
	data.vl$Used <- gsub("(V3| & )", "", data.vl$Used)
	data.vl$my.type <- with(data.vl, paste0(toupper(Type), Censored > 0))
	
	if (is.na(max.vl))
		max.vl <- max(data.vl$VL)
	
	if (use.rainbow) {
		data.vl$my.colour <- as.numeric(data.vl$Date)
		data.vl$my.colour[data.vl$Censored > 0] <- "censored"
	} else {
		data.vl$my.colour <- 'black'
	}
	p.vl <- ggplot(data.vl, aes(x=Date, y=VL))
	if (!is.na(THERAPY_START))
		p.vl <- p.vl + add.therapy(as.numeric(THERAPY_START), Inf, 1)
	if (!is.na(therapy2))
		p.vl <- p.vl + add.therapy(as.numeric(therapy2), as.numeric(THERAPY_START), 2)
	if (!is.na(therapy3))
		p.vl <- p.vl + add.therapy(as.numeric(therapy3), as.numeric(therapy3end), 1)
	p.vl <- p.vl + geom_line(size=vl.linesize)
	p.vl <- p.vl + geom_point2(aes(shape=my.type, colour=my.colour), data=subset(data.vl, Used != ""), size=point.size)
	p.vl <- p.vl + scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*'')))
	
	if (use.real) {
		p.vl <- p.vl + scale_x_date(name="Collection Year", breaks=as.Date(paste0(seq(1984, 2020, by=2), "-01-01")), labels=seq(1984, 2020, by=2))
	} else {
		p.vl <- p.vl + scale_x_continuous(name="Collection Year", breaks=seq(0, 10, by=1))
	}
	
	p.vl <- p.vl + coord_cartesian(ylim=c(10, max.vl))
	p.vl <- p.vl + theme_bw()
	
	p.vl <- p.vl + theme(
		legend.justification=c(1, 1),
		legend.position=0,
		legend.spacing=unit(0, 'cm'),
		legend.margin=margin(.0, .0, .0, .0, 'cm'),
		legend.box.background=element_blank(),
		text=element_text(size=35),
		legend.key.size=unit(1.2, 'cm'),
		axis.text=element_text(size=30, colour='black'),
		axis.text.x=element_text(angle=60, hjust=1),
		legend.text=element_text(size=30),
		panel.grid=element_blank(),
		legend.background=element_blank(),
		legend.key=element_rect(fill="#00000000", colour="#00000000")
	)
	
	if (use.rainbow)
		p.vl <- p.vl + scale_colour_manual(name="", breaks=my.colour.break, limits=my.colour.break, values=my.colour.value)
	else
		p.vl <- p.vl + scale_colour_manual(name="", breaks=c('black'), limits=c('black'), values=c('black'))
	p.vl <- p.vl + scale_shape_manual(name="", breaks=type.break, limits=type.break, values=type.value)
	
	print.plot(pdf.vl.file, p.vl, height=4.2)
}
		
if (mark.latent) {
	data.conf <- as.data.frame(do.call(rbind, lapply(data.all$dist, function(x) unlist(inverse.predict(g, x)))))
	names(data.conf) <- gsub(" ", ".", names(data.conf))
	data.withconf <- as.data.frame(cbind(data.all, data.conf))
	data.withconf$conf.colour <- with(data.withconf, paste0(node.date > Confidence.Limits2, date < Confidence.Limits1))
	conf.colour.break <- c("TRUEFALSE", "FALSETRUE", "FALSEFALSE")
	conf.colour.value <- c('red', 'blue', 'black')
	ptree.withconf <- phylo4d(tree, all.data=data.withconf)

	print.plot(pdf.colour.mark.tree.file,apply.theme(ggtree(ptree.withconf, yscale="node.date", colour="#49494980", size=tree.size, ladderize=F) +
	geom_tippoint(aes(colour=conf.colour, shape=my.type)),
	flipped=T, scaled=T, colour.value.=conf.colour.value, colour.break.=conf.colour.break))
}

if (plot.groups) {
	info.groups <- read.csv(info.file, stringsAsFactors=F)
	
	info.groups <- subset(info.groups, DUPLICATE %in% data.all$tip.label)
	info.split <- split(info.groups, info.groups$TYPE)
	info.groups <- as.data.frame(do.call(rbind, lapply(info.split, function(x) {
		x[sapply(unique(x$DUPLICATE), function(y) which(x$DUPLICATE == y)[1]), ]
	})))
	
	plot.dupes(info.groups, pdf.group.disttree.file, pdf.group.tree.file, pdf.group.hist.file, 0.5, 0.5)
}

if (use.dups) {
	info <- read.csv(info.file, stringsAsFactors=F)
	
	info <- subset(info, DUPLICATE %in% data.all$tip.label)
	
	plot.dupes(info, pdf.dup.disttree.file, pdf.dup.tree.file, pdf.dup.hist.file)
}

