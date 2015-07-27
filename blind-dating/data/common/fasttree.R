library(ape)

fasttree <- function(seq, model='GTR', executable='~/Binaries/FastTreeMP') {
	if(class(seq) != "DNAbin") stop('seq should be of class \'DNAbin\'')

	if(model == 'GTR') {
		model <- '-gtr'
	} else if(model == 'WAG') {
		model <- '-wag'
	} else {
		model <- ''
	}

	# Setup paths 
	executable <- path.expand(executable)
	dnafile <- tempfile()
	treefile <- tempfile()

	# Write tree
	write.dna(seq, dnafile)


	cmd <- sprintf('%s %s -nt %s > %s', executable, model, dnafile, treefile)
	system(cmd)

	read.tree(treefile)
}
