library(ape)
library(parallel)

setwd('~/git/blind-dating/data/simulated/')

source('../common/rtt.R')
source('../common/test.R')
source('../common/queue.R')

# Auxilary functions to read data out of tips
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

# Read al the trees
trees <- trees.read("trees/")
species <- "Simulated"
n.runs <- 1

#
#pdf(sprintf('hist.pdf'), width=11.5, height=8.5)
par(mfrow=c(1, 2), pty="s")
par(mar=c(5.1, 5.1, 2.1, 2.1))

# Show an example of one of these simulated trees
tr <- trees[[1]]

# 
tip.dates <- extract_dates(tr$tip.label)
to.remove <- sort(choose.tips.to.remove(tr, tip.dates, 50))
real <- tip.dates[to.remove]
tip.dates[to.remove] <- NA

missing.indices <- which(is.na(tip.dates))
valid.indices <- which(!is.na(tip.dates))

tr <- rtt(tr, tip.dates)
sub.lengths <- tree.ed(tr)
subs <- sub.lengths[missing.indices]

model <- build.tree.lm(tr, tip.dates, valid.indices)

plot(
	tip.dates[valid.indices], sub.lengths[valid.indices], 
	xlab="Simulation Time", 
	ylab="Expected Number of Subs.", pch=20, cex.lab=1.4, cex.axis=1.2, tck=.01,  axes=F)
#mtext("A) Simulated Data", side=3, adj=0, line=1.1, cex.lab=1.5, font=2); 
points(real, subs, col="red", pch=5,  cex=1.2)
abline(model)
legend(590, 0.160, c("Calibration dates", "Censored dates"), col = c("black", "red"), lty = c(-1, -1), pch = c(20, 5), bg = par("bg"), cex=1.2)#, merge=TRUE)

axis(side=1, at=seq(0, 3000, by=100), tck=.01, cex.axis=1.2)
axis(side=2, at=seq(0, 100, by=0.1), tck=.01, cex.axis=1.2)
box()

par(mar=c(5.1, 5.1, 2.1, 2.1))
plot(c(-1001,-1100), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F, cex.lab=1.4, cex.axis=1.2)
axis(side=1, at=seq(-3, 3, by=0.5), tck=.01, cex.axis=1.2)
axis(side=2, at=seq(0, 3, by=0.5), tck=.01, cex.axis=1.2)
box()

add <- F
errs <- c()
dens <- queue(TRUE)

for(i in 1:n.simulated){
	tip.dates <- extract_dates(trees[[i]]$tip.label)

	min_date <- min(tip.dates)
	max_date <- max(tip.dates)
	tip.dates <- (tip.dates - min_date)/(max_date - min_date)
	rerr <- c()

	for(j in 1:n.runs) {

		err  <- test.rtt.err(trees[[i]], tip.dates, 50)
		rerr <- c(rerr,err)
		errs <- c(errs, err)
	}

	d <- density(rerr, bw=0.15) #// adjust = 10)
	enqueue(dens, d)

	polygon(d, 
		col=rgb(0.5,0.0,0.9,1/(n.simulated)), xlim=c(-0.2,0.2),  border=rgb(1,1,1,0))#, add=T)
		# add <- T
}

abline(h = 0, col = "black", lwd = 1)
abline(v = mean(errs), col = "black", lwd = 2) 
abline(v = median(errs), col = "red", lwd = 2) 

print(sum(errs*errs)/length(errs))
print(median(errs))
print(mean(errs))

#dev.off()


