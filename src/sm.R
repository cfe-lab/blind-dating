library(optparse)
library(ape)
library(slatkin.maddison)
library(parallel)

op <- OptionParser()
op <- add_option(op, "--trees", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--log", type='character')
op <- add_option(op, "--bootstrap", type='numeric', default=1000)
op <- add_option(op, "--threads", type='numeric', default=2)
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
log.file <- args$log
bootstrap <- args$bootstrap - 1
threads <- args$threads

trees <- read.nexus(tree.file)

if (class(trees) == "phylo") {
	trees <- list(trees)
	names(trees) <- "tree"
}
info <- read.csv(info.file)

data <- do.call(
	rbind,
	mclapply(
		1:length(trees),
		function(i) {
			info <- info[match(trees[[i]]$tip.label, info$FULLSEQID), ]
			test <- sm_test(as.factor(info$TYPE), trees[[i]], bootstrap)
			x <- cbind(data.frame(patient=names(trees)[i]), t(test))
			write.csv(x, paste(log.file, names(trees)[i], "csv", sep='.'))
			x
		},
		mc.cores=threads
	)
)

write.csv(data, log.file, row.names=F)