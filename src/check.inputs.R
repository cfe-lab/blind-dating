library(seqinr)
library(optparse)

assert <- function(assertion, assertion.message) {
	if (!all(assertion))
		stop(assertion.message)
}

check.fasta <- function(fasta.file) {
	# can read
	f <- read.fasta(fasta.file)
	
	# ensure alignment
	f.len <- sapply(f, length)
	assert(f.len == f.len[1], "fasta file must contain an alignment")
	
	# ensure no duplicate names
	assert(length(unique(names(f))) == length(f), "duplicate names in fasta file")
	
	# ensure no duplicates
	assert(length(unique(unname(f))) == length(f), "fasta file contains duplicates")
	
	f
}

check.info <- function(info.file) {
	# can read
	info <- read.csv(info.file, stringsAsFactors=F)
	
	# correct columns present
	assert(all(c("FULLSEQID", "COLDATE", "TYPE", "CENSORED") %in% names(info)), "info file missing column")
	
	# correct column types
	assert(is.character(info$FULLSEQID) | is.numeric(info$FULLSEQID), "FULLSEQID in info file parsed incorectly")
	assert(is.character(info$TYPE), "TYPE in info file parsed incorectly")
	assert(is.numeric(info$CENSORED), "CENSORED in info must be numeric")
	assert(is.character(info$DUPLICATE) | is.numeric(info$DUPLICATE), "DUPLICATE in info file parsed incorrectly")
	assert(!is.na(as.Date(info$COLDATE)) | is.numeric(info$COLDATE), "COLDATE must be either in date format or numeric")
		
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
op <- add_option(op, "--freqweights", type='logical', action='store_true', default=F)
op <- add_option(op, "--histbymonth", type='logical', action='store_true', default=F)
op <- add_option(op, "--histfreqby", type='numeric', default=2)
op <- add_option(op, "--histheight", type='numeric', default=1.7)
op <- add_option(op, "--info", type='character', default=NA)
op <- add_option(op, "--liktol", type='numeric', default=1e-3)
op <- add_option(op, "--marklatent", type='logical', action='store_true', default=F)
op <- add_option(op, "--maxcoltime", type='numeric', default=NA)
op <- add_option(op, "--maxvl", type='numeric', default=100000)
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
	args <- parse_args(op, args=settings)
	
	args
}

compare.fasta.info <- function(f, info) {
	TRUE
}

compare.info.settings <- function(info, settings) {
	# check that real is set properly
	assert(!settings$real | !is.na(as.Date(info$COLDATE)), "--real flag is to be used calendar dates")
	assert(settings$real | is.numeric(info$COLDATE), "--real flag is not used with numeric dates")
	
	if (settings$usedups) {
		assert(all(c("DUPLICATE", "COUNT") %in% names(info)), "info file missing columns for usedups option")
		assert(is.numeric(info$COUNT), "COUNT in info must be numeric")
	}
	
	TRUE
}

compare.fasta.info.settings <- function(f, info, settings) {
	# sequence names in info file
	names.f <- if (settings.ogr) c(settings.ogrname, names(f)) else names(f)
	assert(names.f %in% info$FULLSEQID, "info file does not contain all sequence ids in the fasta file")
}

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--settings", type='character')
args <- parse_args(op)

fasta.file <- args$fasta
info.file <- args$info
settings.file <- args$settings
	
# check file integrity
f <- check.fasta(fasta.file)
info <- check.info(info.file)
settings <- check.settings(settings.file)

# check file compatibility
cont <- compare.fasta.info(f, info)
cont <- compare.info.settings(info, settings)
