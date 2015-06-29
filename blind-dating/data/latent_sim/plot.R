library(ape)
library(parallel)

#source('ape.patches.R')
source('include/rtt.R')
source('include/test.R')
source('include/raxml.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) gsub("(u[0-9]+_t[0-9]+)_([0-9\\.]+)$", "\\1", x, perl=T)

n.simulated <- 50

ml.tree <- function(i){
	dna.file <- sprintf("__tmp%d.dna",i)
	unlink(dna.file)

	simenv.dna <- read.FASTA(sprintf("simulated/HIV_%d_out_TRUE.fas", i))
	#simenv.dna <- read.dna(sprintf("simulated/HIV_SIM_%d.phy", i))
	uid <- sample(1:length(names(simenv.dna)), length(names(simenv.dna)), replace=F)
	# read.FASTA is buggy, and creates the DNABin structure incorrectly
	# BUT, it's so broken that we can mess around with the names, then...
	for(j in 1:length(names(simenv.dna))) { 
		# give each bin a unique name 
		names(simenv.dna)[j] <- sprintf("u%s_%s", uid[j], names(simenv.dna)[j])
	}

	# ... output a file (which works for some reason)
	write.dna(simenv.dna, dna.file)
	# and re-read it as a properly formatted DNA Bin
	simenv.dna <- read.dna(dna.file)

	tree <- raxml(simenv.dna, N=100, parsimony.seed=10000, bootstrap.seed=1000)
	write.tree(tree, sprintf("trees/HIV_ml_%d_out.nwk", i))

	unlink(dna.file)
	tree
}

tr <-  read.tree("~/Desktop/paper/blind-dating/data/latent_sim/trees/HIV_ml_13_out.nwk")
plot(tr)

prep.tree <-function(tr){
	latentize.tree(tr, 0.5, function(name){
			tag <- extract_tag(name)
			date <- extract_dates(name)
			ret <- sprintf("%s_PBMC_%s", tag, date)
			ret
	})
}

trees.read <- function(base.path) {
	ml.tree.read <- function(i, path){
		read.tree(sprintf("%s/HIV_ml_%d_out.nwk", path, i))
	}
	mclapply(1:n.simulated, ml.tree.read, base.path)
}

hiv.rna.read <- function(){
	trs <- dir('./7_ml_tree')[grep('_rna.nwk', dir('./7_ml_tree'))]
	ml.tree.read <- function(tr) {
		read.tree(paste('./7_ml_tree', tr, sep='/'))
	}
	lapply(trs, ml.tree.read)
}

#trees <- mclapply(1:n.simulated, ml.tree, mc.cores=1)
#trees <- hiv.rna.read()
trees <- trees.read("trees/")
n.simulated <- length(trees)
run_name <- "HIV RNA"
species <- "Simulated"
n.runs <- 1

#for(remove in c(1, 5, 10, 25)){
pdf(sprintf('hist.pdf', run_name, remove), width=11.5, height=8.5)
add <- F

for(remove in c(25)) {

	for(i in 1:n.simulated){

		if(length(trees[[i]]$tip.label) <= remove) { next }
		tip.dates <- extract_dates(trees[[i]]$tip.label)

		for(j in 1:n.runs) {
			print(remove)
			err  <- test.rtt.err(trees[[i]], tip.dates, remove)
			hist(err, 
				add=add, 
				main="Histogram of Difference ",  
				col=rgb(0.4,0.0,0.2,0.05), 
				xlim=c(-0.2,0.2), 
				freq = T,
				breaks=10)
			# dens <- density(err, adjust=10)
			# lines(dens)

			add <- T
		# 	# if(!is.null(err)) {
		# 	# 	mse[(j-1)*n.simulated + i] <- err
		# 	# }
		}
	}

# 	par(mfrow = c(1, 1))
# 	mean_mse <- mean(mse)
# 	hist(mse, breaks=100,  main="Histogram of MSE", sub=sprintf("(Average MSE [in red]: %s)", signif(mean_mse, digits=6)))
# 	abline(v = mean_mse, col = "red", lwd = 2)
 }
	dev.off()


