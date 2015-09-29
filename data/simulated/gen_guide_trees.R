#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
library(NELSI)

# Number of guide trees to create
n.trees <- 10
n.partitions <- 50
n.replicates <- 1
n.tips <- 100

sim.params <- list(rate = 0.0003, noise = 0.0001)
sim.clockmodel <- simulate.clock

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

sampprob <-c(.99)
lambda <- c(.007)
mu <- c(.007)
times<-c(0)

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob,rho=0,timestop=0)
trees <- lapply(trees, function(x) {unroot(x[[1]])})
sim.trees <- lapply(trees, sim.clockmodel, params=sim.params)
trees <- lapply(sim.trees, date.branches)

indel_control <- 
"
[TYPE] NUCLEOTIDE 2

[SETTINGS]
  [output]                   FASTA 
  [randomseed]               1568746

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
