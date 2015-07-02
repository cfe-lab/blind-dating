library(ape)
library(parallel)

#source('ape.patches.R')
source('include/rtt.R')
source('include/test.R')
source('include/raxml.R')
source('queue.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) (gsub("(.+)_([0-9\\.]+)$", "\\1", x, perl=T))

n.simulated <- 50

trees.read <- function(base.path) {
	ml.tree.read <- function(i, path){
		read.tree(sprintf("%s/HIV_ml_%d_out.nwk", path, i))
	}
	mclapply(1:n.simulated, ml.tree.read, base.path)
}

hiv.rna.read <- function(){
	trs <- dir('./7_ml_tree')[grep('_rna.nwk', dir('./7_ml_tree'))]
	ml.tree.read <- function(tr) {
		read.tree(paste('./7_ml_tree', tr, sep='/'))
	}
	lapply(trs, ml.tree.read)
}

#trees <- mclapply(1:n.simulated, ml.tree, mc.cores=1)
#trees <- hiv.rna.read()
trees <- trees.read("trees/")
n.simulated <- length(trees)
run_name <- "HIV RNA"
species <- "Simulated"
n.runs <- 1

#for(remove in c(1, 5, 10, 25)){
pdf(sprintf('hist.pdf', run_name, remove), width=11.5, height=8.5)
par(mfrow=c(1, 2), pty="s")

## Show an example of one of these simulated trees
tr <- trees[[1]]
tip.dates <- extract_dates(tr$tip.label)


sub.lengths <- tree.ed(tr)

to.remove <- sort(choose.tips.to.remove(tr, tip.dates, 25))
real <- tip.dates[to.remove]
tip.dates[to.remove] <- NA

missing.indices <- which(is.na(tip.dates))
valid.indices <- which(!is.na(tip.dates))
subs <- sub.lengths[missing.indices]

tr <- rtt(tr, tip.dates)
model <- build.tree.lm(tr, tip.dates, valid.indices)

plot(
	tip.dates[valid.indices], sub.lengths[valid.indices], 
	xlab="Time", 
	ylab="Expected Number of Subs.")
mtext("A) Simulated Data", side=3, adj=0, line=1.1, cex=1.5, font=2); 
points(real, subs, col="red")
abline(model)
##

plot(c(1,1), xlim=c(-0.1,0.1), ylim=c(0, 25), xlab="Normalized Error", ylab="Density")
add <- F
errs <- c()
dens <- queue(TRUE)

for(remove in c(25)) {
	for(i in 1:n.simulated){
		if(length(trees[[i]]$tip.label) <= remove) { next }
		tip.dates <- extract_dates(trees[[i]]$tip.label)

		min_date <- min(tip.dates)
		max_date <- max(tip.dates)
		tip.dates <- (tip.dates - min_date)/(max_date - min_date)
		rerr <- c()


		for(j in 1:n.runs) {

			err  <- test.rtt.err(trees[[i]], tip.dates, remove)
			rerr <- c(rerr,err)
			errs <- c(errs, err)
		}

		d <- density(rerr, bw=0.015) #// adjust = 10)
		enqueue(dens, d)
		polygon(d, 
			col=rgb(0.5,0.0,0.9,1/(n.simulated)), xlim=c(-0.2,0.2),  border=rgb(1,1,1,0), add=T)
			# add <- T
	}
}

abline(h = 0, col = "black", lwd = 1)
abline(v = median(errs), col = "black", lwd = 2) 

dev.off()


