library(ape)
library(ggtree)
library(phylobase)
source("~/git/node.dating/src/node.dating.R")

args <- commandArgs(trailingOnly = T)

tree.file <- args[1]
pat.id <- args[2]
time.tree.file <- args[3]

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")
time.data.file <- paste0("stats/", pat.id, ".timedata.csv")

tree <- ape::ladderize(read.tree(tree.file))
data <- read.csv(data.file, col.names=c("tip.label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)
stats <- read.csv(stats.file, stringsAsFactors=F)

mu <- stats$Model.Slope

stree <- drop.tip(tree, which(data$censored == 1))
while (length(grep("^N", stree$tip.label)) > 0) stree <- drop.tip(stree, grep("^N", stree$tip.label))

m <- match(c(stree$tip.label, stree$node.label), c(tree$tip.label, tree$node.label))

sdates <- estimate.dates(stree, data$date[m[1:length(stree$tip.label)]], mu, nsteps=0, show.steps=1000, lik.tol=1e-5)

new.dates <- rep(NA, length(tree$tip.label) + tree$Nnode)
new.dates[m] <- sdates
plotting.dates <- new.dates
plotting.dates[1:length(tree$tip.label)] <- data$date

while (any(is.na(new.dates))) {
	new.dates[is.na(new.dates)] <- unlist(lapply(which(is.na(new.dates)), function (i) {
		e <- tree$edge[, 2] == i
		p <- tree$edge[e, 1]
		if (is.na(new.dates[p]))
			NA
		else
			new.dates[p] + tree$edge.length[e] / mu
	}))
}

plotting.dates <- estimate.dates(tree, plotting.dates, mu, m, nsteps=0, show.steps=1000, lik.tol=1e-5)

data.all <- as.data.frame(cbind(rbind(data, data.frame(tip.label=paste0("N.", 1:tree$Nnode), type="NODE", censored=NA, date=NA, dist=node.depth.edgelength(tree)[1:tree$Nnode + nrow(data)], est.date=NA, date.diff=NA, ci.high=NA, ci.low=NA)), new.date=new.dates, plotting.date=plotting.dates))

ptree <- phylo4d(tree, all.data=data.all)





time.tree <- tree
time.tree$edge.length <- apply(tree$edge, 1, function(x) new.dates[2] - new.dates[1])
write.tree(time.tree, time.tree.file)

data <- as.data.frame(cbind(data, new.date=new.dates[1:length(tree$tip.label)], new.date.diff=new.dates[1:length(tree$tip.label)]-data$date))
write.table(data, time.data.file, col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Estimated Date", "Date Difference", "New Dating", "New Dating Differnce"), row.names=F, sep=",")