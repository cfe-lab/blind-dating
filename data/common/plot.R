library(ape)
library(parallel)

source('../common/rtt.R')
source('../common/test.R')
source('../common/raxml.R')
source('../common/queue.R')

args <- commandArgs(trailing = T)


tree.dir <- args[1]
data.type <- as.integer(args[2])	# 0 = no dna, 1 = mixed (latency unknown), 2 = mixed (latency known)
use.rtt <- as.integer(args[3])		# 0 = no, 1 = yes (only plasma), 2 = yes (all)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)?$", "\\2", x, perl=T)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) (gsub("(.+)_([0-9\\.]+)$", "\\1", x, perl=T))
extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
extract_real_pbmc_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)_PBMC_([0-9\\.]+)$", "\\2", x, perl=T))
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

trees.read <- function(base.path){
	trs <- dir(base.path)[grep('.tre', dir(base.path))]
	
	print(trs)
	
	ml.tree.read <- function(tr) {	
		tree = read.tree(paste(base.path, tr, sep='/'))
	
		if (sum("REFERENCE"==tree$tip.label) > 0) {
			if (use.rtt)			
				drop.tip(tree, "REFERENCE")
			else
				drop.tip(root(tree, "REFERENCE"), "REFERENCE")
		}
		else
			tree
			
	}
	
	lapply(trs, ml.tree.read)
}

plot.hist <- function(err, colour, pos) {
		par(mfg=c(1, pos))
		
		d <- density(err, bw=0.15)
		polygon(d, col=colour, xlim=c(-0.1,0.2),  border=rgb(1,1,1,0))
}

trees <- trees.read(tree.dir)
run_name <- "HIV RNA"
species <- "Simulated"
n.runs <- 1

if (data.type == 2) {
	pdf('plot.pdf', width=17.5, height=8.5)
	dev.plot = dev.cur()

	pdf('hist.pdf', width=11.5, height=8.5)
	par(mfrow=c(1, 2), pty="s")
	dev.hist = dev.cur()
} else {
	pdf('plot.pdf', width=11.5, height=8.5)
	dev.plot = dev.cur()
	
	pdf('hist.pdf', width=6.5, height=8.5)
	par(mfrow=c(1, 1), pty="s")
	dev.hist = dev.cur()
}

par(mar=c(5.1, 6.1, 4.1, 2.1))
plot(c(-001,-100), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F)
axis(side=1, at=seq(-3, 3, by=0.5), tck=.01)
axis(side=2, at=seq(0, 3, by=0.5), tck=.01)
box()

if (data.type == 2) {
	par(mar=c(5.1, 6.1, 4.1, 2.1))
	plot(c(-201,-300), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F)
	axis(side=1, at=seq(-3, 3, by=0.5), tck=.01)
	axis(side=2, at=seq(0, 3, by=0.5), tck=.01)
	box()
}


add <- F
errs <- c()
serrs <- c()
rerrs <- c()
srerrs <- c()
dens <- queue(TRUE)

# for data
patients = c()
seq.ids = c()
times <- c()
dists <- c()
r.times <- c()

# for stats
rmses <- c()
medians <- c()
scales <- c()
#slows <- c()

count <- 1

#tpertime <- c()
#glmslope <- c()
#plasmarmse <- c()

for(tree in trees){
	print(count)
	
#	for(k in 1:n.runs) {
		#####
		#
		#####
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
		
		scale_val <- max_date - min_date

#		tip.dates <- (tip.dates - min_date)/(max_date - min_date)
		
		# Setup plasma dates and re-root tree
		plasma.dates <- tip.dates
		
		if (use.rtt == 1)
			plasma.dates[tip.pbmc] <- NA
		
		if (use.rtt) {
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
	
		model <- glm(plasma.dists ~ plasma.dates)
				
		a<-model$coefficients[[1]]
		b<-model$coefficients[[2]]
		
		err <- (pbmc.s.dates - (pbmc.dists/b - a/b)) / scale_val
		 		
		dev.set(dev.hist)
		
		plot.hist(err, rgb(0.5, 0.0, 0.9, 1/(length(trees))), 1)
		
		errs <- c(errs, err)
		
		# fill data data.frame information
		patients <- c(patients, rep(count, sum(tip.pbmc)))
		seq.ids <- c(seq.ids, gsub("^(.+)_([0-9\\.]+)_PBMC_([0-9\\.]+)$", "\\1", tree$tip.label[tip.pbmc], perl=T))
		times <- c(times, pbmc.s.dates)
		dists <- c(dists, pbmc.dists)
		serrs <- c(serrs, pbmc.s.dates - (pbmc.dists/b - a/b))

		if (data.type == 2) {
			tip.real.dates <- extract_real_pbmc_dates(tree$tip.label)
						
			pbmc.r.dates <- tip.real.dates[tip.pbmc] #real dates
			rerr <- (pbmc.r.dates - (pbmc.dists/b - a/b)) / scale_val
			
			plot.hist(rerr, rgb(0.5, 0.0, 0.9, 1/(length(trees))), 2)
			
			rerrs <- c(rerr, rerrs)
			
			r.times <- c(r.times, pbmc.r.dates)
			srerrs <- c(srerrs, pbmc.r.dates - (pbmc.dists/b - a/b))
		}		
		dev.set(dev.plot)
		
		
		if (data.type == 2) {
			par(mfrow=c(1, 3), pty="s")
		} else {
			par(mfrow=c(1, 2), pty="s")
		}
		
		# make plot showing latency date estimation
		par(cex=1, mar=c(5,5,2,2)+0.1)
		if (data.type==2) {			
			plot.xlim <- c(min(c(tip.dates, pbmc.r.dates)), max(c(tip.dates, pbmc.r.dates)))
		}
		else {			
			plot.xlim <- c(min(tip.dates), max(tip.dates))
		}
		plot.ylim <- c(min(c(plasma.dists, pbmc.dists)), max(c(plasma.dists, pbmc.dists)))
		
		plot(plasma.dates, plasma.dists, xlab="Time", ylab=sprintf("Divergence from root"), xlim=plot.xlim, ylim=plot.ylim,  pch=20, cex.lab=1.3, cex.axis=1.1, col='grey')
		points(pbmc.s.dates, pbmc.dists, col='red', pch=5, lty=2)
		abline(model, lty=2)
		
#		tryfx <- function() {
#			model2 <- glm(pbmc.dists ~ pbmc.s.dates)
#			abline(model2, lty=2, col='red')
#		}
		
#		try(tryfx())

#		tpertime <- c(tpertime, length(tree$tip.label) / max(plasma.dates))
#		glmslope <- c(glmslope, b)
		
#		perr <- plasma.dates - (plasma.dists/b - a/b)
		
#		plasmarmse <- c(plasmarmse, sqrt(sum(perr*perr))/length(perr))

		if (data.type == 2) {
		
			points(pbmc.r.dates, pbmc.dists, col="red", pch=18,  cex=1)
	
			for (i in 1:length(pbmc.dists)) {
				lines(x=c(pbmc.s.dates[i], pbmc.r.dates[i]), y=rep(pbmc.dists[i], 2), col='red')
			}
	
#			legend(1200, 0.36, c("Sample dates (RNA)", "Latency dates (DNA)", "Sample dates (DNA)"), col = c("grey", "red", rgb(1,0,0,0.5)), pch = c(20, 18, 5), cex=1)
		}
		else {
#			legend(1200, 0.36, c("Sample dates (RNA)", "Sample dates (DNA)"), col = c("grey", rgb(1,0,0,0.5)), pch = c(20, 5), cex=1)
		}
		
		par(mar=c(5.1, 6.1, 4.1, 2.1))
		plot(c(-1001,-1100), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F)
		axis(side=1, at=seq(-3, 3, by=0.5), tck=.01)
		axis(side=2, at=seq(0, 3, by=0.5), tck=.01)
		box()

		if (data.type == 2) {
			par(mar=c(5.1, 6.1, 4.1, 2.1))
			plot(c(-1101,-1200), xlim=c(-1.0,1.0), ylim=c(0, 3), xlab="Normalized Error", ylab="Density", axes=F)
			axis(side=1, at=seq(-3, 3, by=0.5), tck=.01)
			axis(side=2, at=seq(0, 3, by=0.5), tck=.01)
			box()
		}
		
		rmse_val <- sqrt(sum(err*err)/length(err))
		median_val <- median(err)
#		slow_val <- sum(err > 0)/length(err)
		
		rmses <- c(rmses, rmse_val)
		medians <- c(medians, median_val)
		scales <- c(scales, scale_val)
#		slows <- c(slows, slow_val)
		
		plot.hist(err, rgb(0.5, 0.0, 0.9, 1), 2)
		
		abline(h = 0, col = "black", lwd = 1)
		abline(v = median_val, col = "red", lwd = 2)
#		abline(v = median_val - rmse_val, col = "black", lty = 2) 
#		abline(v = median_val + rmse_val, col = "black", lty = 2) 

		
		if (data.type == 2) {
			mse_val <- sqrt(sum(rerr*rerr)/length(rerr))
			median_val <- median(rerr)
#			slow_val <- sum(rerr > 0)/length(rerr)
		
			rmses <- c(rmses, rmse_val)
			medians <- c(medians, median_val)
			scales <- c(scales, scale_val)
#			slows <- c(slows, slow_val)
		
			plot.hist(rerr, rgb(0.5, 0.0, 0.9, 1), 3)
			
			abline(h = 0, col = "black", lwd = 1)
			abline(v = median_val, col = "red", lwd = 2)
#			abline(v = median_val - rmse_val, col = "black", lty = 2) 
#			abline(v = median_val + rmse_val, col = "black", lty = 2) 
		}
#	}

	count <- count + 1
}

# write data data.frame
if (data.type == 2) {
	df2 <- data.frame(Patient=patients, Seq.ID=seq.ids, Time=times, Distance=dists, Error=serrs, Normal.Error=errs, Real.Time=r.times, Real.Error=srerrs, Real.Normal.Error=rerrs)
} else {
	df2 <- data.frame(Patient=patients, Seq.ID=seq.ids, Time=times, Distance=dists, Error=serrs, Normal.Error=errs)
}

write.csv(df2, file="data.csv", row.names=FALSE)


dev.set(dev.hist)

rmse_val <- sqrt(sum(errs*errs)/length(errs))
median_val <- median(errs)
#slow_val <- sum(errs > 0) / length(errs)

rmses <- c(rmses, rmse_val)
medians <- c(medians, median_val)
scales <- c(scales, 1)
#slows <- c(slows, slow_val)
				
print(sprintf("Test: median %f, RMSE %f",
	median_val, rmse_val))

par(mfg=c(1, 1))
abline(h = 0, col = "black", lwd = 1)
abline(v = median_val, col = "red", lwd = 2) 
#abline(v = median_val - rmse_val, col = "black", lty = 2) 
#abline(v = median_val + rmse_val, col = "black", lty = 2) 


if (data.type == 2) {
	rmse_val <- sqrt(sum(rerrs*rerrs)/length(rerrs))
	median_val <- median(rerrs)
#	slow_val <- sum(rerrs > 0) / length(rerrs)
	
	rmses <- c(rmses, rmse_val)
	medians <- c(medians, median_val)
	scales <- c(scales, 1)
#	slows <- c(slows, slow_val)
	
	print(sprintf("Overall (real): median %f, rmse %f",
		median_val, rmse_val))

	par(mfg=c(1, 2))
	abline(h = 0, col = "black", lwd = 1)
	abline(v = median_val, col = "red", lwd = 2)
#	abline(v = median_val - rmse_val, col = "black", lty = 2) 
#	abline(v = median_val + rmse_val, col = "black", lty = 2) 
}


#df <- data.frame(Median=medians, RMSE=rmses, Scale=scales, Percent.Slow=slows)
df <- data.frame(Median=medians, RMSE=rmses, Scale=scales)

write.csv(df, file="stats.csv", row.names=FALSE)


#df2 <- data.frame(TipsPerTime=tpertime, GLMSlope=glmslope, PlasmaRMSE=plasmarmse)

#write.csv(df2, file="fit.csv", row.names=FALSE, col.names=TRUE, sep='\t')
	
warnings()

dev.off()
dev.off()
