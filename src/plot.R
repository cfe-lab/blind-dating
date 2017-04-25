library(ape)
library(phylobase)
library(ggtree)

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

#source(file.path(source.dir, 'read.info.R'), chdir=T)

get.date <- function(date) {
	as.character(
		if (DATES_AS_DATES)
			as.Date(date, origin=as.Date("1970-01-01"))
		else
			date
	)
}

apply.axes <- function(p, flipped, scaled) {
	if (scaled) {
		if (flipped) {
			p + y.scale +
				scale_x_continuous(name="Divergence from root", breaks=seq(0.0, 0.20, by=.04), limits=c(dist.min, dist.max)) +
				geom_abline(intercept=stats[,"Model.Intercept"], slope=stats[, "Model.Slope"], color="#0060b0b0", linetype=2) +
				geom_hline(yintercept=13361, colour="#60600080", linetype=2)
		} else {
			p + x.scale +
				scale_y_continuous(name="Divergence from root", breaks=seq(0.0, 0.20, by=.04), limits=c(dist.min, dist.max)) +
				geom_abline(intercept=-stats[,"Model.Intercept"]/stats[,"Model.Slope"], slope=1/stats[, "Model.Slope"], color="#0060b0b0", linetype=2) +
				geom_vline(xintercept=13361, colour="#60600080", linetype=2)
		}
	} else {
		if (flipped) {
			p + scale_x_continuous(name="Divergence from root", breaks=seq(0.04, 0.20, by=.04), limits=c(dist.min, dist.max)) +
				theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
		}
		else {
			p + scale_y_continuous(name="Divergence from root", breaks=seq(0.04, 0.20, by=.04), limits=c(dist.min, dist.max)) +
				theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
		}
	}
}

apply.theme <- function(p, flipped=F, scaled=T) {
	apply.axes(p + theme_bw() +
		theme(
			legend.position=c(.15,.8),
			panel.grid.major = element_blank(),
		    panel.grid.minor = element_blank()
		) +
		scale_colour_manual(name="", breaks=type.break, labels=type.label, values=type.value, limits=type.break) +
		scale_shape_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(16, 18), limits=c(0, 1)) +
		scale_size_manual(name="", breaks=c(0, 1), labels=c("Training", "Censored"), values=c(2.5, 3), limits=c(0, 1)) +
		guides(colour=type.guide),
		flipped, scaled)
}

dist.min <- -.005
dist.max <- .27
LIK_TOL <- 1e-1

args <- commandArgs(trailingOnly = T)

tree.file <- args[1]
pat.id <- args[2]
type <- as.numeric(args[3]) # 0 Training/Censored, 1 RNA/DNA

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")
pdf.file <- paste0("plots/", pat.id, ".pdf")
pdf.disttree.file <- paste0("plots/", pat.id, ".disttree.pdf")
pdf.tree.file <- paste0("plots/", pat.id, ".tree.pdf")
pdf.hist.file <- paste0("plots/", pat.id, ".hist.pdf")

tree <- ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
stats <- read.csv(stats.file, , stringsAsFactors=F)

node.dates <- estimate.dates(tree, data$date, 1/stats[,"Model.Slope"], lik.tol=LIK_TOL, show.steps=1000, nsteps=0)

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type=rep("NODE", tree$Nnode), censored=rep(NA, tree$Nnode), date=rep(NA, tree$Nnode), dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], est.date=rep(NA, tree$Nnode), date.diff=rep(NA, tree$Nnode))), node.date=node.dates))

ptree <- phylo4d(tree, all.data=data.all)

if (type == 1) {
	x.scale <- scale_x_continuous(name="Year", breaks=c(9131, 10957, 12784, 14610, 16436), labels=c("1995", "2000", "2005", "2010", "2015"), limits=c(8766, 17532))
	y.scale <- scale_y_continuous(name="Year", breaks=c(9131, 10957, 12784, 14610, 16436), labels=c("1995", "2000", "2005", "2010", "2015"), limits=c(8766, 17532))
	type.guide <- guide_legend(override.aes=list(shape=15, size=5))
	type.break <- c("PLASMA", "PBMC", "PBMC (cultured)", "WHOLE BLOOD")
	type.label <- c("Plasma RNA", "PBMC DNA", "Cultured PBMC DNA", "Whole Blood DNA")
	type.value <- c('black', 'red', 'cyan', 'darkgreen')
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