library(ape)
library(parallel)

#source('ape.patches.R')
source('include/rtt.R')
# source('include/test.R')
source('include/raxml.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))

hiv.rna.read <- function(){
	trs <- dir('./trees.good')[grep('.tre', dir('./trees.good'))]
	ml.tree.read <- function(tr) {
		drop.tip(read.tree(paste('./trees.good', tr, sep='/')), "REFERENCE")
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

pdf(sprintf('hist.pdf'), width=11.5, height=8.5)
par(mfrow=c(1, 2), pty="s")

# ## Show an example of one of these simulated trees
tree <- trees[[1]]

# tip.dates <- extract_dates(tree$tip.label)
# print(tip.dates)

# # Extract the data from the tip labels
# tip.dates <- extract_dates(tree$tip.label)
n <- length(tree$tip.label)

# Extract the data from the tip labels
tip.dates <- extract_dates(tree$tip.label)
tip.types <- types(tree$tip.label)

# Mark the PBMC and PLASMA cells
tip.pbmc <- tip.types == "PBMC"
tip.plasma <- tip.types == "PLASMA"

plasma.dates <- tip.dates
plasma.dates[tip.pbmc] <- NA

# tree <- rtt(tree, plasma.dates)
plasma.dates <- tip.dates[tip.plasma]
pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates

distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
plasma.dists <- distances[tip.plasma]
pbmc.dists <- distances[tip.pbmc]

model <- glm(plasma.dists ~ plasma.dates)

a<-model$coefficients[[1]]
b<-model$coefficients[[2]]

# # (pbmc.dists/b - a/b)

plot(
	plasma.dates, plasma.dists, 
	xlab="Time (days)", 
	ylab="Expected Number of Subs.")
mtext("D) DNA and RNA", side=3, adj=0, line=1.1, cex=1.5, font=2); 
points(pbmc.s.dates, pbmc.dists, col="red")
abline(model)
# ##

plot(c(100,100), xlim=c(-1,1), ylim=c(0, 2.5), xlab="Normalized Error", ylab="Density")
# add <- F
errs <- c()

for(i in 1:n.simulated) {
	for(k in 1:1) {
		tree <- trees[[i]]

		tip.dates <- extract_dates(tree$tip.label)
		tip.types <- types(tree$tip.label)


		tip.pbmc <- tip.types == "PBMC"
		tip.plasma <- tip.types == "PLASMA"

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
		d <- density(err, bw=0.2)
		polygon(d, col=rgb(0.5, 0.0, 0.9, 1/(n.simulated)), xlim=c(-0.1,0.2),  border=rgb(1,1,1,0))
		errs <- c(err, errs)
	}
}

abline(h = 0, col = "black", lwd = 1)
abline(v = mean(errs), col = "black", lwd = 2) 
abline(v = median(errs), col = "red", lwd = 2) 
warnings()
dev.off()


