library(ape)
library(parallel)

#source('ape.patches.R')
setwd('~/git/blind-dating/data/ancre/')
source('../common/rtt.R')
source('../common/test.R')
source('../common/raxml.R')
source('../common/fasttree.R')
source('../common/queue.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))

hiv.rna.read <- function(){
	trs <- dir('./tree')[grep('.tre', dir('./tree'))]
	ml.tree.read <- function(tr) {
		drop.tip(root(read.tree(paste('./tree', tr, sep='/')), "REFERENCE"), "REFERENCE")
	}
	lapply(trs, ml.tree.read)
}

extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) (gsub("(.+)_([0-9\\.]+)$", "\\1", x, perl=T))
extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

trees <- hiv.rna.read()
n.simulated <- length(trees)
n.runs <- 10

#pdf(sprintf('hist.pdf'), width=11.5, height=8.5)
#par(mfrow=c(1, 2), pty="s")
par(mar=c(5,4,2,2)+0.1)

## Show an example of one of these simulated trees
tree <- trees[[1]]

tip.dates <- extract_dates(tree$tip.label)
print(tip.dates)

# Extract the data from the tip labels
tip.dates <- extract_dates(tree$tip.label)
n <- length(tree$tip.label)

spercent <- floor(n*0.25)

tip.plasma <- rep(T, n)
tip.plasma[sample(1:n, spercent)] <- F
tip.pbmc <- !tip.plasma

plasma.dates <- tip.dates
plasma.dates[tip.pbmc] <- NA

tree <- rtt(tree, plasma.dates)
plasma.dates <- tip.dates[tip.plasma]
pbmc.s.dates <- tip.dates[tip.pbmc] 

distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
plasma.dists <- distances[tip.plasma]
pbmc.dists <- distances[tip.pbmc]

model <- glm(plasma.dists ~ plasma.dates)

a<-model$coefficients[[1]]
b<-model$coefficients[[2]]

plot(
	jitter(plasma.dates), plasma.dists, 
	xlab="Time (days)", 
	ylab="Expected Number of Subs.", pch=20,  cex=1.2, tck=.01,  axes=F)
#mtext("C) RNA Only", side=3, adj=0, line=1.1, cex=1.5, font=2)
points(jitter(pbmc.s.dates), pbmc.dists, col="red", pch=5,  cex=1.2)
abline(model)
##

legend(0, 0.28, c("Calibration dates", "Censored dates"), col = c("black", "red"), lty = c(-1, -1), pch = c(20, 5), merge = TRUE, bg = par("bg"))#, cex=1.2)

axis(side=1, at=seq(0, 3000, by=100), tck=.01)
axis(side=2, at=seq(0, 100, by=0.05), tck=.01)
box()



par(mar=c(5.1, 6.1, 4.1, 2.1))
plot(c(-1001,-1100), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F)
axis(side=1, at=seq(-3, 3, by=0.5), tck=.01)
axis(side=2, at=seq(0, 3, by=0.5), tck=.01)
box()

add <- F
errs <- c()

for(i in 1:n.simulated){
	for(k in 1:n.runs*3) {
		tree <- trees[[i]]

		tip.dates <- extract_dates(tree$tip.label)
		tip.types <- types(tree$tip.label)

		n <- length(tree$tip.label)
		spercent <- floor(n*0.25)
		tip.plasma <- rep(T, n)
		tip.plasma[sample(1:n, spercent)] <- F
		tip.pbmc <- !tip.plasma

		if(length(which(tip.pbmc==F)) <= 2){ next }

		min_date <- min(tip.dates)
		max_date <- max(tip.dates)

		tip.dates <- (tip.dates - min_date)/(max_date - min_date)

		plasma.dates <- tip.dates
		plasma.dates[tip.pbmc] <- NA
			
		tree <- rtt(tree, plasma.dates)

		plasma.dates <- tip.dates[tip.plasma]
		pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates

		distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
		plasma.dists <- distances[tip.plasma]
		pbmc.dists <- distances[tip.pbmc]

		model <- glm(plasma.dists ~ plasma.dates)

		a<-model$coefficients[[1]]
		b<-model$coefficients[[2]]

		err <- pbmc.s.dates - (pbmc.dists/b - a/b)
		d <- density(err, bw=0.15)
		polygon(d, col=rgb(0.5, 0.0, 0.9, 1/(n.simulated*n.runs)), xlim=c(-0.1,0.2),  border=rgb(1,1,1,0))
		errs <- c(err, errs)
	}
}

abline(h = 0, col = "black", lwd = 1)
abline(v = mean(errs), col = "black", lwd = 2) 
abline(v = median(errs), col = "red", lwd = 2) 

print(sum(errs*errs)/length(errs))
print(sprintf("Mean %s, median %s", mean(errs), median(errs)))

warnings()
dev.off()


