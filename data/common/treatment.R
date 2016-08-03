library(plotrix)
library(lme4)

make.hist <- function(pdf.file, df) {
	pdf(pdf.file, width=7, height=7)
	par(mfrow=c(1, 1), pty="s")
	dev.hist = dev.cur()
	
	par(mar=c(5, 5, 2, 2), ps=12)
	plot(c(-001,-100), xlim=c(min(df$Error), max(df$Error)), ylim=c(0, 10), xlab="Difference", ylab="Density", axes=F, cex.lab=1)
	axis(side=1, at=seq(floor(min(df$Error)), ceiling(max(df$Error)), by=100), tck=.01, cex.axis=1)
	axis(side=2, at=seq(0, 10, by=0.5), tck=.01, cex.axis=1)
	box()
	
	polygon(df$Error, col=rgb(0.5, 0.0, 0.9, 1), xlim=c(-0.1,0.2), border=rgb(1,1,1,0))
			
	dev.off()
}

args <- commandArgs(trailingOnly=T)

pat.id <- args[1]
tree.path <- args[2]
treatment.start <- as.numeric(args[3])
treatment.end <- as.numeric(args[4])

# plot
#source(file.path(source.dir, 'rtt.R'), chdir=T)
#source(file.path(source.dir, 'test.R'), chdir=T)
#source('../common/raxml.R', chdir=T)
#source('../common/queue.R', chdir=T)
library(ape)

args <- commandArgs(trailing = T)


tree.path <- file.path(tree.path, paste0("patient_", pat.id, ".tre"))
data.type <- 1
use.rtt <- 2
out.pdf <- paste0("patient_", pat.id, ".pdf")
leg.lab <- 0

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)?$", "\\2", x, perl=T)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) (gsub("(.+)_([0-9\\.]+)$", "\\1", x, perl=T))
extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
extract_real_pbmc_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)_PBMC_([0-9\\.]+)$", "\\2", x, perl=T))
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

tree <- read.tree(tree.path)
tree <- if (sum("REFERENCE"==tree$tip.label) > 0) {
			if (use.rtt)			
				drop.tip(tree, "REFERENCE")
			else
				drop.tip(root(tree, "REFERENCE"), "REFERENCE")
		} else {
			tree
		}

#pdf(out.pdf, width=8, height=8)

##
		tip.dates <- extract_dates_pp(tree$tip.label)
		tip.types <- types(tree$tip.label)
		
		if (data.type == 0) {
			n.tips = length(tree$tip.label)
			tips <- sort(sample(1:n.tips, n.tips*.5))
			tip.types = rep("PLASMA", n.tips)	
			tip.types[tips] = "PBMC"
		} else {	
			tip.types <- types(tree$tip.label)
		}
		
		# Mark the PBMC and PLASMA cells
		tip.pbmc <- tip.types == "PBMC"
		tip.plasma <- tip.types == "PLASMA"

		# Normalize dates
		min_date <- min(tip.dates)
		max_date <- max(tip.dates)
		
		scale.val <- max_date - min_date
		
		# Setup plasma dates and re-root tree
		plasma.dates <- tip.dates
		
		if (use.rtt == 1)
			plasma.dates[tip.pbmc] <- NA
		
		if (use.rtt == 2) {
			tree <- rtt(tree, plasma.dates)
			
			# reset values after rerooting
			tip.dates <- extract_dates_pp(tree$tip.label)
			tip.types <- types(tree$tip.label)
		
			if (data.type == 0) {
				n.tips = length(tree$tip.label)
				tips <- sort(sample(1:n.tips, n.tips*.5))
				tip.types = rep("PLASMA", n.tips)	
				tip.types[tips] = "PBMC"
			} else {	
				tip.types <- types(tree$tip.label)
			}
			
			# Mark the PBMC and PLASMA cells
			tip.pbmc <- tip.types == "PBMC"
			tip.plasma <- tip.types == "PLASMA"	
		}	

		plasma.dates <- tip.dates[tip.plasma]
		pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates
				
		distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
		plasma.dists <- distances[tip.plasma]
		pbmc.dists <- distances[tip.pbmc]
		
		treatment.stage <- as.integer(tip.dates > treatment.start) + as.integer(tip.dates > treatment.end)
			
		model <- lmer(distances ~ tip.dates + (tip.dates | treatment.stage))
				
#		a<-model$coefficients[[1]]
#		b<-model$coefficients[[2]]
					
		# make plot showing latency date estimation
		par(cex=1, mar=c(5,5,2,2))
		if (data.type==2) {
			tip.real.dates <- extract_real_pbmc_dates(tree$tip.label)
						
			pbmc.r.dates <- tip.real.dates[tip.pbmc] #real dates
		
			plot.xlim <- c(min(c(tip.dates, pbmc.r.dates)), max(c(tip.dates, pbmc.r.dates)))
		} else {
			plot.xlim <- c(min(tip.dates), max(tip.dates))
		}
		plot.ylim <- c(min(c(plasma.dists, pbmc.dists)), max(c(plasma.dists, pbmc.dists)))
		
		y.intersp=1
		cex=1
				
		plot(plasma.dates, plasma.dists, xlab="Time (days)", ylab=sprintf("Divergence from root"), xlim=plot.xlim, ylim=plot.ylim,  pch=20, cex.lab=1.2, cex.axis=1, col='black', cex=cex)
		points(pbmc.s.dates, pbmc.dists, col='red', pch=5, lty=2, cex=cex)
#		abline(model, lty=2)
		
		latent <- function(x, r=75) {
			draw.circle(pbmc.s.dates[x], pbmc.dists[x], radius=r, col="#00000000", border="#003ecc")
			lines(c((pbmc.dists[x] - a) / b, pbmc.s.dates[x] - r), y=rep(pbmc.dists[x], 2), col="#003ecc", lty=2)
		}
		
		legend.labels <- if (data.type == 2) {
				l.x <- plot.xlim[2] - (plot.xlim[2] - plot.xlim[1])*.395
				c("Sample dates (uncensored)", "Integration dates (censored)", "Sample dates (censored)")
			} else {
				if (leg.lab) {
					l.x <- plot.xlim[2] - (plot.xlim[2] - plot.xlim[1])*.3825
					c("Sample dates (uncensored)", "Sample dates (censored)")
				} else {
					l.x <- plot.xlim[2] - (plot.xlim[2] - plot.xlim[1])*.295
					c("Sample dates (RNA)", "Sample dates (DNA)")
				}
			}
		
		if (data.type == 2) {		
			points(pbmc.r.dates, pbmc.dists, col="red", pch=18,  cex=cex)
	
			for (i in 1:length(pbmc.dists)) {
				lines(x=c(pbmc.s.dates[i], pbmc.r.dates[i]), y=rep(pbmc.dists[i], 2), col='red')
			}
					
			l.y <- (plot.ylim[2] - plot.ylim[1])*.12 + plot.ylim[1]
					
			legend(l.x, l.y, legend.labels, col = c("black", "red", "red"), pch = c(20, 18, 5), y.intersp=y.intersp, cex=cex)
		} else {
			l.y <- (plot.ylim[2] - plot.ylim[1])*.08 + plot.ylim[1]
		
			legend(l.x, l.y, legend.labels, col = c("black", "red"), pch = c(20, 5), y.intersp=y.intersp, cex=cex)
		}
				
abline(v=treatment.start, col=rgb(0, .8, 0, 1), lty=2, lwd=1)
abline(v=treatment.end, col=rgb(0, .8, 0, 1), lty=2, lwd=1)
			
		
#dev.off()


# read data
df <- read.csv(paste0("data.rttall.csv"), header=T)
good <- read.csv("good_fits.txt", header=F)
names(good) <- c("id", "name")
df.pat <- df[df$Patient == good[good$name == as.numeric(pat.id), 1], ]

# binomial test
df.pat$Is.Latent <- df.pat$Error > 0
bin.test <- binom.test(sum(df.pat$Is.Latent), length(df.pat$Is.Latent), alternative="greater")
bin.test

# histogram
make.hist(paste0("hist_", pat.id, ".pdf"), df.pat)

