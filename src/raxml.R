library(treeio)
library(ape)

raxml <- function(dnafile, parsimony.seed=NULL, bootstrap.seed=NULL, executable='raxmlHPC', threads=10, N=10, name='raxml', clear=TRUE, model="GTRGAMMA", tmp=TRUE, tmp.dir=if (tmp) tempdir() else NULL) {
	if(is.null(parsimony.seed)) {
		parsimony.seed <- as.integer(sample(2**31,1))
		warning(sprintf('parsimony.seed should be fixed for debugging!\nSetting to %d!', parsimony.seed))
	}

	if(is.null(bootstrap.seed)) {
		bootstrap.seed <- as.integer(sample(2**31,1))
		warning(sprintf('bootstrap.seed should be fixed for debugging!\nSetting to %d!', bootstrap.seed))
	}

	dnafile <- normalizePath(dnafile)

	if (tmp) {
		old.wd = normalizePath(getwd())
		new.wd = tmp.dir
		setwd(new.wd)
	}
	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}

	wd <- normalizePath(getwd())
	r.name <- paste0(name,"_R")
	b.name <- paste0(name,"_B")
	threads <- if (threads > 1) sprintf("-T %d", threads) else "" 

	cmd <- sprintf('%s %s -b %d -m %s -p %d -N %d -s %s -n %s -w %s -O', 
					executable, threads, bootstrap.seed, model, parsimony.seed, N, dnafile, r.name, wd)

	cmd2 <- sprintf('%s %s -f d -m %s  -s %s -N %d -n %s -w %s -p %d -O', 
					executable, threads, model,  dnafile, N, b.name, wd, parsimony.seed)
		
	cmd3 <- sprintf('%s %s -f b -n %s -m %s -t %s/RAxML_bestTree.%s -z %s/RAxML_bootstrap.%s -s %s -w %s -O', 
					executable, threads, name, model, wd, b.name, wd, r.name, dnafile, wd)

	cat("First\n")
	system(cmd)
	cat("Second\n")
	system(cmd2)
	cat("Third\n")
	system(cmd3)

	tr <- read.tree(sprintf('RAxML_bipartitions.%s', name))
	
	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}
	if(tmp) {
		setwd(old.wd)
		if (clear)
			unlink(new.wd)
	}
	
	tr
}
