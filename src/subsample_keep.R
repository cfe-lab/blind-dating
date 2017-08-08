library(seqinr)

get.time.points <- function(i) {
	s <- sampler[, i]
	
	which(dates %in% s)
}

args <- commandArgs(trailingOnly = T)

INFO_FILE <- args[1]
SUFFIX <- args[2]
METHOD <- as.integer(args[3]) # 0: N RNA time points 1: all but 1 time point, 2: 2 RNA time points, -2: all but 2 time points
REPS <- as.integer(args[4]) # number of replicates generated for method 0
TIME_POINTS <- as.integer(args[5]) # number of time points, N, for method 0
SEED <- as.integer(args[6])

set.seed(SEED)

cat("Reading...\n")
info <- read.csv(INFO_FILE, stringsAsFactors=F)
types <- info$CENSORED
dates <- info$COLDATE
time.points <- unique(dates)
n.points <- length(time.points)

cat("Choosing...\n")
method <- if (METHOD == 0) {
	sampler <- combn(time.points, TIME_POINTS)
	s <- sampler[, sample(ncol(sampler), REPS)] 
	reps <- REPS

	get.time.points
} else if (METHOD == 1) {
	sampler <- combn(time.points, n.points - 1)
	reps <- ncol(sampler)
	
	get.time.points
} else if (METHOD == 2) {
	sampler <- combn(time.points, 2)
	reps <- ncol(sampler)
	
	get.time.points
} else {
	sampler <- combn(time.points, n.points - 2)
	reps <- ncol(sampler)
	
	get.time.points
}

cat("Filtering...\n")
filters <- lapply(1:reps, method)
cat("Writing...\n")
suppress <- lapply(1:reps, function(i) {
	info.i <- info
	info.i$CENSORED[-filters[[i]]] <- 1
	info.i$CENSORED[filters[[i]]] <- 0
	write.csv(info.i, paste0(SUFFIX, i, ".csv"), row.names=F)
})