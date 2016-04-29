#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
library(NELSI)
library(Rmpfr)
#library(phangorn)

#source('../common/rtt.R')
#source('../common/raxml.R')
#source('../common/queue.R')

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

source(file.path(source.dir, 'test.R'), chdir=T)

args <- commandArgs(trailingOnly = T)

unlatency.rate <- as.double(args[1])
latency.rate <- as.double(args[2])
r.seed <- as.integer(args[3])
indelible.seed <- args[4]

if (r.seed != 0) {
	set.seed(r.seed)
}

file.remove("latency.csv")

# Number of guide trees to create
# TO DO: quantify variance
n.trees <- 50
n.partitions <- n.trees
n.replicates <- 1
n.tips <- 100

# From BEAST
clock.rate <- 0.0001964
noise.rate <- 0.00001417

# birth-death model parameters, from BEAST
sampprob <-c(0.005237)
lambda <- c(0.05116)
mu <- c(0.05006)
times<-c(0)

sim.params <- list(rate = clock.rate, noise = noise.rate)
sim.clockmodel <- simulate.clock

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

date.branches <- function(s.tree) {
	tree <- s.tree$phylogram
	subsr <- s.tree$tree.data.matrix[,6]
	times <- s.tree$tree.data.matrix[,7]

	# Calculate the cumulative times
	tree$edge.length <- times
	times <- node.depth.edgelength(tree)

	tree$edge.length <- subsr

	dates <- unlist(Map(toString,times))[1:n.tips]
	tree$tip.label <- paste(tree$tip.label, dates, sep='_')
	tree
}

# move choosing tips outside of making tips latent to sync RNG of latent and unlatent methods
choose.tips <- function(
	s.tree,
	percent = 0.5
) {
  tr <- s.tree$phylogram

  # #  
  n.tips <- length(tr$tip.label)
  n.mod <- as.integer(n.tips*percent)  # number of tips to modify
  delta <- rep(0, n.tips)
  types <- rep("PLASMA", n.tips)
 
  # #
  tips <- sort(sample(1:n.tips, n.mod))  # indices of tips to modify
  
  # #
  s.tree.new <- s.tree
  s.tree.new$latent.tips <- tips
  
  s.tree.new
}

# TO DO: probably do a range of values here
make.latent <- function(
	s.tree
#	percent=0.5
#	latency.rate=0.0001
) {
#  rate = clock.rate
#  noise = noise.rate
  tr <- s.tree$phylogram

  # #  
  n.tips <- length(tr$tip.label)
#  n.mod <- as.integer(n.tips*percent)  # number of tips to modify
#  delta <- rep(0, n.tips)
  types <- rep("PLASMA", n.tips)
 
  # #
  tips <- s.tree$latent.tips
#  tips <- sort(sample(1:n.tips, n.mod))  # indices of tips to modify
#  scale <- rbeta(n.mod, 1, 100)  # left-skew (median about 0.006)

  # #
  edge.length <- s.tree$tree.data.matrix[,7]
#  edges <- which(tr$edge[,2] %in% tips)  # post-traversal indices of edges to modify
  edges <- unlist(lapply(tips, function(x) {which(tr$edge[,2] == x)}))  # post-traversal indices of edges to modify (preserves order of tips)
#  edge.length <- tr$edge.length  # branch lengths in unit time
  
  # latency distribution
  rlatent <- function(sc) {
    u <- mpfr(rexp(length(sc), unlatency.rate), 20)
    print(u)
    t.0 <- (sc - u) * ((sc - u) > 0)
#   print(t.0)
    r <- exp(-latency.rate * t.0)-exp(-latency.rate * mpfr(sc, 20))
#    print(r)
    x <- mpfr(runif(length(sc)), 20)
#    print(x)
    
    as.double(-log(exp(-latency.rate * t.0) - r * x) / latency.rate)
  }
    
  types[tips] <- "PBMC"
 
  edge.mod <- edge.length
#  edge.mod[edges] <- edge.mod[edges]*scale
#  print(edge.mod[edges])
  
  edge.mod[edges] <- if (latency.rate > 1E-100 && unlatency.rate > 1E-100)
      rlatent(edge.mod[edges])
    else 
      if (unlatency > 1E-100) {
        edge.mod[edges] <- edge.mod[edges] - 100
        bad <- edge.mod[edges] <= 0
        edge.mod[edges][bad] <- .1 * (edge.mod[edges][bad] + 100)
        edge.mod[edges]
      }
      else
        edge.mod[edges]
#  print(edge.mod[edges])
  
        
#  delta <- edge.length - edge.mod
  write.table(data.frame(id=tr$tip.label[tips], length=edge.length[edges], latency=edge.length[edges]-edge.mod[edges]), "latency.csv", append=T, row.names=F, sep=",")
  
  # # 
  latent <- edge.mod
  actual <- edge.length 
  
  actual.nozero <- actual
  actual.nozero[actual.nozero == 0] <- 1

  # # .
  # Assume evolution rate is constant along the last edge
  latent.evo <- s.tree$tree.data.matrix[,6]*latent/actual.nozero
  actual.evo <- s.tree$tree.data.matrix[,6]
#  err <- rnorm(length(tr$edge.length), mean = 0, sd = noise)
#  latent.evo <- abs(latent * rate + err)  # convert time to exp. sub'ns
#  actual.evo <- abs(actual * rate + err)
  
  # # #
  tr$edge.length <- latent
  times.l <- node.depth.edgelength(tr)[1:n.tips]
  tr$edge.length <- latent.evo
  evo.l <- node.depth.edgelength(tr)[1:n.tips]
  
  tr$edge.length <- actual
  times.a <- node.depth.edgelength(tr)[1:n.tips]
  tr$edge.length <- actual.evo
  evo.a <- node.depth.edgelength(tr)[1:n.tips]
  
  # 
  tr$tip.label <- paste(tr$tip.label, times.l, types, times.a, sep="_")
  tr$edge.length <- latent.evo
  
  tr
}

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob, rho=0, timestop=0)
trees <- lapply(trees, function(x) {unroot(x[[1]])})
sim.trees <- lapply(trees, sim.clockmodel, params=sim.params)

#if (latent) {
sim.trees <- lapply(sim.trees, choose.tips)
trees <- lapply(sim.trees, make.latent)
#} else {
#	trees <- lapply(sim.trees, date.branches)
#}

print(length(trees))

indel_control <- sprintf(
"
[TYPE] NUCLEOTIDE 2

[SETTINGS]
  [output]                   FASTA 
  [randomseed]               %s

[MODEL]    HKY_HIV
  [submodel] HKY 8.5               
  [statefreq] 0.42 0.15 0.15 0.28
"
, indelible.seed)


for(i in 1:n.trees) {
	tree_dat <- write.tree(trees[[i]])
	print(tree_dat)
	indel_control <- paste0(indel_control, sprintf("[TREE] tree_%d %s \n", i, tree_dat))
}

#index <- sample(1:n.trees, n.partitions, replace = (n.partitions > n.trees))
index <- 1:n.trees
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%d [tree_%d HKY_HIV 700] \n", i, index[i]))
}

indel_control <- paste0(indel_control, "[EVOLVE] \n")
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("    pHKY_%d %d HIV_%d_out \n", i, n.replicates, i))
}

write(indel_control, 'control.unfixed.txt')
