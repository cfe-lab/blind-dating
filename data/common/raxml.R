library(ape)

raxml <- function(seq, parsimony.seed=NULL, bootstrap.seed=NULL, executable='~/Binaries/raxmlHPC', threads=12, N=10, name='raxml', clear=TRUE, model="GTRGAMMA") {
	if(class(seq) != "DNAbin") stop('seq should be of class \'DNAbin\'')
	if(is.null(parsimony.seed)) {
		parsimony.seed <- as.integer(sample(2**31,1))
		warning(sprintf('parsimony.seed should be fixed for debugging!\nSetting to %d!', parsimony.seed))
	}

	if(is.null(bootstrap.seed)) {
		bootstrap.seed <- as.integer(sample(2**31,1))
		warning(sprintf('bootstrap.seed should be fixed for debugging!\nSetting to %d!', bootstrap.seed))
	}

	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}

	wd <- path.expand(getwd())
	executable <- path.expand(executable)
	r.name <- paste0(name,"_R")
	b.name <- paste0(name,"_B")

	dnafile <- tempfile()
	write.dna(seq, dnafile)

	cmd <- sprintf('%s -T %d -b %d -m %s -p %d -N %d -s %s -n %s -w %s -O', 
					executable, threads, bootstrap.seed, model, parsimony.seed, N, dnafile, r.name, wd)

	cmd2 <- sprintf('%s -T %d -f d -m %s  -s %s -N %d -n %s -w %s -p %d -O', 
					executable, threads, model,  dnafile, N, b.name, wd, parsimony.seed)
	
	cmd3 <- sprintf('%s -T %d -f b -n %s -m %s -t %s/RAxML_bestTree.%s -z %s/RAxML_bootstrap.%s -s %s -w %s -O', 
					executable, threads, name, model, wd, b.name, wd, r.name, dnafile, wd)

	system(cmd)
	system(cmd2)
	system(cmd3)

	tr <- read.tree(sprintf('RAxML_bipartitions.%s', name))
	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}
	unlink(dnafile)
	tr
}