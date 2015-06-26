library(Hmisc)

test.rtt <- function(t, tip.dates, remove = 1, random = T) {
	to.remove <- if(!random) 
	remove
	else
	sample(length(tip.dates), remove, replace = F)

	mdates <- tip.dates
	mdates[to.remove] <- (NA)

	if(length(unique(mdates[!is.na(mdates)])) <= 1) {
		return(NULL)
	}

	t <- rtt(t, mdates)

	
	
	sum((t$tip.dates - tip.dates)^2)/remove
}

test.n.plot.rtt <- function(t, tip.dates, remove = 10, random = T, name="") {
	to.remove <- if(!random) 
	remove
	else
	sample(length(tip.dates), remove, replace = F)

	mdates <- tip.dates
	mdates[to.remove] <- (NA)
	if(length(unique(mdates[!is.na(mdates)])) <= 1) {
		return(NULL)
	}
	print(mdates)
	t <- rtt(t, mdates)
	msre <- sum((t$tip.dates - tip.dates)^2)/remove
	plot.branchlen.date(t, mdates, tip.dates, msre, name)
}

plot.branchlen.date <- function(tr, tip.dates, tip.dates.recon=NULL, err=0, name="") {
	missing.indices <- which(is.na(tip.dates))
	valid.indices <- which(!is.na(tip.dates))

	tip.lengths <- node.depth.edgelength(tr)

	distances <- tip.lengths[valid.indices]
	times <- tip.dates[valid.indices]
	model <- lm(times ~ distances)

	par(mfrow=c(1, 2), pin=c(5,4), mar=c(10, 8, 8, 2) + 0.1)
	plot(distances, times, main="Time vs Distance", sub=sprintf("(Mean square error: %s)", signif(err, digits=6)),
		xlab="Expected Substitutions", ylab="Time")
	abline(model)

	if(!is.null( tip.dates.recon)){
		distancesm <- tip.lengths[missing.indices]
		timesm <- tip.dates.recon[missing.indices]

		ptimes <- predict(model, data.frame(distances=tip.lengths[missing.indices]))
		delta <- timesm-ptimes

		errbar(distancesm, timesm, timesm, ptimes,cap=0,col="red", add=T)
		text(distancesm, timesm+(delta/2), labels=signif(abs(delta),digits=4), cex=0.5)
		points(distancesm, timesm, col="red")
	}
	par(pin=c(4,6), mar=c(5, 4, 4, 2) + 0.1)
	tr$tip.label[missing.indices] <- rep("???", length(missing.indices))
	plot(tr,cex=0.75,main=name)
}