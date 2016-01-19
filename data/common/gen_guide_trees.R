#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
library(NELSI)

source('../common/rtt.R')
source('../common/test.R', chdir = T)
source('../common/raxml.R')
source('../common/queue.R')

args <- commandArgs(trailingOnly = T)

latent <- as.integer(args[1])
r.seed <- as.integer(args[2])
indelible.seed <- args[3]

if (r.seed != 0) {
	set.seed(r.seed)
}

# Number of guide trees to create
# TO DO: quantify variance
n.trees <- 50
n.partitions <- 100
n.replicates <- 1
n.tips <- 100

# env evolution rate from  Alizon and Fraser, "Within-host and between-host evolutionary rates across the HIV-1 genomes" (2013) and
# Perelson et al., "HIV-1 Dynamics in Vivo: Virion Clearance Rate, Infected Cell Life-Span, and Viral Generation Time" (1996)
clock.rate <- 0.0003 # 0.00010296
noise.rate <- 0.0001 # 0.00009050

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

# TO DO: probably do a range of values here
make.latent <- function(
	tr, 
	percent=0.5,
	latency.rate=0.001
) {
  rate = clock.rate
  noise = noise.rate

  # #  
  n.tips <- length(tr$tip.label)
  n.mod <- as.integer(n.tips*percent)  # number of tips to modify
  delta <- rep(0, n.tips)
  types <- rep("PLASMA", n.tips)

 
  # #
  tips <- sort(sample(1:n.tips, n.mod))  # indices of tips to modify
  scale <- rbeta(n.mod, 1, 100)  # left-skew (median about 0.006)

  # #
  edges <- which(tr$edge[,2] %in% tips)  # post-traversal indices of edges to modify
  edge.length <- tr$edge.length  # branch lengths in unit time

  
  # I think it is more realistic to use a conditional exponential
  # distribution:
  # let Y be the waiting time to lineage sampling
  # let t be the waiting time to the lineage going latent
  # then the probability distribution of (t) is
  #   m exp(-mt) / (1-exp(-L y))
  # where L is the sampling rate
  # and m is the latency rate
  #v <- rexpo()
 # v <- lapply(edge.length[edges], function(x) {rexp(1, latency.rate/(1-exp(-clock.rate*x))}
  
  types[tips] <- "PBMC"
 
  edge.mod <- edge.length 
  edge.mod[edges] <- edge.mod[edges]*scale
  delta <- edge.length - edge.mod
  
  # # 
  latent <- edge.length - delta  # = edge.mod
  actual <- edge.length 

  # # .
  err <- rnorm(length(tr$edge.length), mean = 0, sd = noise)
  latent.evo <- abs(latent * rate + err)  # convert time to exp. sub'ns
  actual.evo <- abs(actual * rate + err)
  
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


# TO DO: need to fit
# Values from Perelson et al., "HIV-1 Dynamics in Vivo: Virion Clearance Rate, Infected Cell Life-Span, and Viral Generation Time" (1996)
sampprob <-c(.99)
lambda <- c(.007) # c(0.761843113)
mu <- c(.007) # c(0.0333333)
times<-c(0)

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob,rho=0,timestop=0)
trees <- lapply(trees, function(x) {unroot(x[[1]])})
sim.trees <- lapply(trees, sim.clockmodel, params=sim.params)

if (latent)
	trees <- lapply(trees, make.latent)
else
	trees <- lapply(sim.trees, date.branches)

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

index <- sample(1:n.trees, n.partitions, replace = (n.partitions > n.trees))
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%d [tree_%d HKY_HIV 700] \n", i, index[i]))
}

indel_control <- paste0(indel_control, "[EVOLVE] \n")
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("    pHKY_%d %d HIV_%d_out \n", i, n.replicates, i))
}

write(indel_control, 'control.unfixed.txt')
