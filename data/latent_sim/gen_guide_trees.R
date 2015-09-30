#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
library(NELSI)

setwd('/Users/art/git/blind-dating/data/latent_sim/')
source('../common/rtt.R')
source('../common/test.R')
source('../common/raxml.R')
source('../common/queue.R')

# Number of guide trees to create
n.trees <- 10
n.partitions <- 100
n.replicates <- 1
n.tips <- 100

clock.rate <- 0.0003

# sim.params <- list(rate = 0.0003, noise = 0.001)
sim.params <- list(rate = clock.rate, noise = 0.0001)
sim.clockmodel <- simulate.clock


trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

make.latent <- function(tr, percent=0.5, rate=0.0003, noise=0.0001) {
  # #
  n.tips <- length(tr$tip.label)
  n.mod <- as.integer(n.tips*percent)  # number of tips to modify
  delta <- rep(0, n.tips)
  types <- rep("PLASMA", n.tips)
  
  # #
  tips <- sort(sample(1:n.tips, n.mod))  # indices of tips to modify
  scale <- rbeta(n.mod, 1, 100)
  types[tips] <- "PBMC"
  
  # #
  edges <- which(tr$edge[,2] %in% tips)
  edge.length <- tr$edge.length 
  
  edge.mod <- edge.length 
  edge.mod[edges] <- edge.mod[edges]*scale
  delta <- edge.length - edge.mod
  
  
  # # 
  latent <- edge.length - delta
  actual <- edge.length 

  # # 
  err <- rnorm(length(tr$edge.length), mean = 0, sd = noise)
  latent.evo <- abs(latent * rate + err)
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


sampprob <-c(.99)
lambda <- c(.007)
mu <- c(.007)
times<-c(0)

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob,rho=0,timestop=0)
trees <- lapply(trees, function(x)  {make.latent(unroot(x[[1]]))})



indel_control <- 
"
[TYPE] NUCLEOTIDE 2

[SETTINGS]
  [output]                   FASTA 
  [randomseed]               1568745

[MODEL]    HKY_HIV
  [submodel] HKY 8.5               
  [statefreq] 0.42 0.15 0.15 0.28
"

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

write(indel_control, 'control.txt')
