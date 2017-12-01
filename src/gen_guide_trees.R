#!/usr/bin/Rscript

# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)
#library(NELSI)


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

suffix <- args[1]
r.seed <- as.integer(args[2])
indelible.seed <- args[3]
reference <- args[4]

if (r.seed != 0) {
	set.seed(r.seed)
}

indel_control <- ""

if (suffix == "" || suffix == "head") {
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
}

if (suffix != "head") {
# Number of guide trees to create
# TO DO: quantify variance
n.tips <- 100

# From GLM
clock.rate <- 1.96e-4
noise.rate <- 1.42e-5

# birth-death model parameters, from BEAST
lambda <- 5.12e-2
mu <- 5.01e-2
sampprob <- 5.24e-3
times<-c(0)

sim.params <- list(rate = clock.rate, noise = noise.rate)
sim.clockmodel <- simulate.clock

cat("Building tree...\n")
 
rooted.tree <- sim.bdsky.stt(n.tips, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob, rho=0, timestop=0)[[1]]
tree <- unroot(rooted.tree)
stree <- sim.clockmodel(tree, params=sim.params)

data <- data.frame(PATIENT=paste0("SIM_", suffix), SEQID=tree$tip.label, FULLSEQID=tree$tip.label, COLDATE=node.depth.edgelength(tree)[1:n.tips], TYPE="PLASMA", CENSORED=0, KEPT=1, DUPLICATE="", NOTE="")
write.csv(data, paste0("info/SIM_", suffix, ".csv"), row.names=F)

write.tree(tree, sprintf("trees.ori/SIM_%s.time.nwk", suffix))

tree <- stree$phylogram
tree$edge.length <- stree$tree.data.matrix[, 6]
#tree$edge.length <- tree$edge.length * clock.rate

cat("Writing tree to file...\n")

tree_dat <- write.tree(tree)
write.tree(tree, sprintf("trees.ori/SIM_%s.nwk", suffix))
indel_control <- paste0(indel_control, sprintf("[TREE] tree_%s %s \n", suffix, tree_dat))

indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%s [tree_%s HKY_HIV %s] \n", suffix, suffix, reference))


indel_control <- paste0(indel_control, sprintf("    pHKY_%s %d simulated/SIM_%s \n", suffix, 1, suffix))
}

write(indel_control, sprintf('control.unfixed.%s.txt', suffix))