library(optparse)
library(ape)

assert <- function(assertion, assertion.message) {
	if (!all(assertion))
		stop(assertion.message)
}

check.tree <- function(tree.file) {
	# can read
	tree <- read.tree(tree.file)
	
	# ensure binary tree
	assert(is.binary(tree), "tree must be binary")
	
	# ensure tip labels have correct length
	tip.labs <- tree$tip.labels
	assert(sum(!is.na(tree$tip.label)) == tree$Nnode + 1, "tree tip labels not complete")
	
	# ensure no negative edge lengths
	assert(tree$edge.length > 0, "tree must have only positive edge lengths")
	
	# ensure no duplicate names
	assert(length(unique(tip.labels)) == length(tip.labels), "duplicate names in tree")
	
	f
}

check.info <- function(info.file) {
	# can read
	info <- read.csv(info.file, stringsAsFactors=F)
	
	# correct columns present
	assert(c("FULLSEQID", "COLDATE", "TYPE", "CENSORED") %in% names(info), "info file missing column")
	
	# correct column types
	assert(is.character(info$FULLSEQID) | is.numeric(info$FULLSEQID), "FULLSEQID in info file parsed incorectly")
	assert(is.character(info$TYPE), "TYPE in info file parsed incorectly")
	assert(is.numeric(info$CENSORED), "CENSORED in info must be numeric")
		
	info
}

check.settings <- function(settings.file) {
	# can read
	settings <- readLines(settings.file)
	
	# filter comments
	settings <- settings[grepl("^[^#]", settings)]
	
	# ensure only useable options present
	op <- OptionParser()
op <- add_option(op, "--cartoon", type='logical', action='store_true', default=F)
op <- add_option(op, "--cutoff", type='character', default=NA)
op <- add_option(op, "--distby", type='double', default=NA)
op <- add_option(op, "--distmax", type='double', default=NA)
op <- add_option(op, "--distmin", type='double', default=NA)
op <- add_option(op, "--dnashapescale", type='numeric', default=1)
op <- add_option(op, "--dupshift", type='numeric', default=0.01)
op <- add_option(op, "--weight", type='character', default=NA)
op <- add_option(op, "--histbymonth", type='logical', action='store_true', default=F)
op <- add_option(op, "--histfreqby", type='numeric', default=2)
op <- add_option(op, "--histheight", type='numeric', default=1.7)
op <- add_option(op, "--info", type='character', default=NA)
op <- add_option(op, "--liktol", type='numeric', default=1e-3)
op <- add_option(op, "--marklatent", type='logical', action='store_true', default=F)
op <- add_option(op, "--maxcoltime", type='numeric', default=NA)
op <- add_option(op, "--maxvl", type='numeric', default=NA)
op <- add_option(op, "--method", type='character', default="correlation")
op <- add_option(op, "--mincoltime", type='numeric', default=NA)
op <- add_option(op, "--nsteps", type='numeric', default=1000)
op <- add_option(op, "--ogrname", type='character', default="REFERENCE")
op <- add_option(op, "--ogr", type='logical', action='store_true', default=F)
op <- add_option(op, "--raxml", type='character', default="raxml")
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--seed", type='numeric', default=1989)
op <- add_option(op, "--therapy2", type='character', default=NA)
op <- add_option(op, "--therapy3end", type='character', default=NA)
op <- add_option(op, "--therapy3", type='character', default=NA)
op <- add_option(op, "--therapy", type='character', default=NA)
op <- add_option(op, "--threads", type='numeric', default=2)
op <- add_option(op, "--tmpdir", type='character', default=NULL)
op <- add_option(op, "--useall", type='logical', action='store_true', default=F)
op <- add_option(op, "--usedups", type='logical', action='store_true', default=F)
op <- add_option(op, "--vlfile", type='character', default=NA)
op <- add_option(op, "--xtitle", type='character', default="Collection Year")
op <- add_option(op, "--yearby", type='double', default=NA)
op <- add_option(op, "--yearend", type='character', default=NA)
op <- add_option(op, "--yearstart", type='character', default=NA)
op <- add_option(op, "--outputfolder", type='character', default="plots")
op <- add_option(op, "--types", type='character', default="PLASMA,PBMC")
op <- add_option(op, "--typevalues", type='character', default="16,1,18,5")
op <- add_option(op, "--rainbow", type='logical', action='store_true', default=T)
op <- add_option(op, "--black", type='logical', action='store_false', dest="rainbow")
op <- add_option(op, "--plotdups", type='logical', action='store_true', default=F)
op <- add_option(op, "--colourvalues", type='character', default=NA)
op <- add_option(op, "--histdupfreqby", type='numeric', default=10)
op <- add_option(op, "--training", type='numeric', default=0)
op <- add_option(op, "--alltraining", type='logical', action='store_true', default=F)
op <- add_option(op, "--plotgroups", type='logical', action='store_true', default=F)
op <- add_option(op, "--latentedges", type='logical', action='store_true', default=F)
op <- add_option(op, "--textheight", type='character', default="2016-01-01")
op <- add_option(op, "--model", type='character', default="GTRGAMMA")
op <- add_option(op, "--vlyearby", type='numeric', default=2)
	args <- parse_args(op, args=settings)
	
	args
}

compare.tree.info <- function(tree, info) {
	TRUE
}

compare.tree.settings <- function(tree, info) {
	# check that outgroup is in tree file
	if (settings$ogr)
		assert(settings$ogrname %in% tree$tip.label, "outgroup not in tree")
		
	TRUE
}

compare.info.settings <- function(info, settings) {
	# check that real is set properly
	if (settings$real)
		assert(is.character(info$COLDATE) & !is.na(as.Date(info$COLDATE)), "use --real flag with calendar dates")
	else
		assert(is.numeric(info$COLDATE), "--real flag is not used with numeric dates")
	
	if (settings$usedups) {
		assert("DUPLICATE" %in% names(info), "info file missing columns for usedup option")
		assert(is.character(info$DUPLICATE) | is.numeric(info$DUPLICATE), "DUPLICATE in info file parsed incorrectly")
	}
	
	if (!is.na(settings$weight)) {
		assert(settings$weight %in% names(info), "info file missing columns for usedup option")
		assert(is.numeric(info[, settings$weight]), "weight in info must be numeric")
	}
	
	TRUE
}

compare.tree.info.settings <- function(f, info, settings) {
	# sequence names in info file
	seq.ids <- if (settings$ogr) c(settings$ogrname, info$FULLSEQID) else info$FULLSEQID
	assert(tree$tip.label %in% seq.ids, "info file does not contain all sequence ids in the tree")
	
	TRUE
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--settings", type='character')
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
settings.file <- args$settings
	
# check file integrity
f <- check.tree(tree.file)
info <- check.info(info.file)
settings <- check.settings(settings.file)

# check file compatibility
cont <- compare.tree.info(tree, info)
cont <- compare.tree.settings(tree, settings)
cont <- compare.info.settings(info, settings)
cont <- compare.tree.info.settings(tree, info, settings)
