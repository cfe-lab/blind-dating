library(Hmisc)
source("../common/queue.R")

tree.ed <- function(tree) {
	l <- length(tree$tip.label)
	node.depth.edgelength(tree)[1:l]
}

# Builds a linear model from a given tree
# discarding elements in discard (a list of 
# 1:tip.length of truth values)
build.tree.lm <- function(tree, tip.dates, keep=NULL) {
	
	# Input is time
	time <- tip.dates
	response <- tree.ed(tree)

	if(!is.null(keep)) {
		time <- tip.dates[keep]
		response <- response[keep]
	}

	glm(response ~ time)
}

# Predicts time from distance
predict.time <- function(model, dists) {
	a<-model$coefficients[[1]]
	b<-model$coefficients[[2]]

	(dists/b - a/b)
}

choose.tips.to.remove <- function(t, tip.dates, remove = 1, random = T) {

	to.remove <- if(!random) 
	remove
	else
	sample(length(tip.dates), remove, replace = F)
}

##
get.affected.tips <- function(t, edge) {
	to.remove <- queue(T)
	tips <- queue(T)

	is.tip <- function(t, i) {
		if( i <= t$Nnode+1){
			return(T);
		}
		return(F);
	}

	enqueue(to.remove, t$edge[edge, 2])
	while(length(to.remove$q) > 0) {
		node <- dequeue(to.remove)

		if(is.tip(t, node)){
			enqueue(tips, node)
		} else {
			for(next_edg in t$edge[t$edge[,1] == node,2]){
				enqueue(to.remove, next_edg)
			}
		}
	}
	unlist(tips$q)
}

latentize.tree <- function(tr, scale, rename_funct=NULL) {
	tips <- c()

	for(s in scale){
		# Choose a branch at random
		br <- sample(1:(tr$Nnode*2), 1)
		# print(br)

		# Scale that branch
		tr$edge.length[br] = tr$edge.length[br] * s
		tips <- unique(c(tips, get.affected.tips(tr, br)))
	}

	if(!is.null(rename_funct)) {
		for(tip in tips) {
			tr$tip.label[tip] <- rename_funct(tr$tip.label[tip])
		}
	}
	tr
}

#### 
# OLD STUFF
####

test.rtt.err <- function(t, tip.dates, remove = 1, random = T) {
	to.remove <- if(!random) 
	remove
	else
	sample(length(tip.dates), remove, replace = F)

	mdates <- tip.dates
	mdates[to.remove] <- (NA)

	if(length(unique(mdates[!is.na(mdates)])) <= 1) {
		return(NULL)
	}

	t <- rtt(t, mdates, ncpu=8)

	(t$tip.dates - tip.dates)
}

test.rtt <- function(t, tip.dates, remove = 1, random = T) {
	sum(test.rtt.err(t, tip.dates, remove, random)^2)/remove
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
	model <- lm(distances ~ times) # times ~ distances)

	#  
	a<-model$coefficients[[1]]
  	b<-model$coefficients[[2]]

	par(mfrow=c(1, 2), pin=c(5,4), mar=c(10, 8, 8, 2) + 0.1)
	plot(distances, times, main="Time vs Distance", sub=sprintf("(Mean square error: %s)", signif(err, digits=6)),
		xlab="Expected Substitutions", ylab="Time")

	abline(model)

	if(!is.null( tip.dates.recon)){
		distancesm <- tip.lengths[missing.indices]
		timesm <- tip.dates.recon[missing.indices]

		ptimes <- (tip.lengths[missing.indices]/b - a/b)
		delta <- timesm-ptimes

		errbar(distancesm, timesm, timesm, ptimes,cap=0,col="red", add=T)
		text(distancesm, timesm+(delta/2), labels=signif(abs(delta),digits=4), cex=0.5)
		points(distancesm, timesm, col="red")
	}
	par(pin=c(4,6), mar=c(5, 4, 4, 2) + 0.1)
	tr$tip.label[missing.indices] <- rep("???", length(missing.indices))
	plot(tr,cex=0.75,main=name)
}