#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
library(NELSI)

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

args <- commandArgs(trailingOnly = T)

r.seed <- as.integer(args[1])
indelible.seed <- args[2]

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

# From GLM
clock.rate <- 3.94e-5
noise.rate <- 2.00e-6

# birth-death model parameters, from BEAST
lambda <- 8.87e-3
mu <- 7.75e-3
sampprob <- 0.115
times<-c(0)

sim.params <- list(rate = clock.rate, noise = noise.rate)
sim.clockmodel <- simulate.clock

cat("Building trees...\n")

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob, rho=0, timestop=0)
trees <- lapply(trees, function(x) {unroot(x[[1]])})
sim.trees <- lapply(trees, sim.clockmodel, params=sim.params)

suppress <- lapply(1:n.trees, function(i) {data <- data.frame(PATIENT=paste0("SIM_", i), SEQID=trees[[i]]$tip.label, FULLSEQID=trees[[i]]$tip.label, COLDATE=node.depth.edgelength(tree)[1:50], CENSORED=0, KEPT=1, DUPLICATE="", NOTE=""); write.csv(data, paste0("info/SIM_", i, ".csv"), row.names=F)})

tree <- sim.tree$phylogram
tree$edge.length <- sim.tree$tree.data.matrix[, 6]


cat("Preparing INDELible file...\n")

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


for (i in 1:n.trees) {
	tree_dat <- write.tree(trees[[i]])
	indel_control <- paste0(indel_control, sprintf("[TREE] tree_%d %s \n", i, tree_dat))
}

#index <- sample(1:n.trees, n.partitions, replace = (n.partitions > n.trees))
index <- 1:n.trees
for (i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%d [tree_%d HKY_HIV 700] \n", i, index[i]))
}

indel_control <- paste0(indel_control, "[EVOLVE] \n")
for (i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("    pHKY_%d %d HIV_%d_out \n", i, n.replicates, i))
}

write(indel_control, 'control.unfixed.txt')
