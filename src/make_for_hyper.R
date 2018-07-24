library(seqinr)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "read.info.R"), chdir=T)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--hyperfasta", type='character')
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_optoon(op, "--useall", type='logical', action='store_true', default=F)
args <- parse_args(op)

fasta.file <- args$fasta
info.file <- args$info
hyper.fasta.file <- args$hyperfasta
use.date <- args$real
use.all <- args$use.all

f <- read.fasta(fasta.file)

if (!use.all) {
	info <- read.info(info.file, names(f))

	if (use.date)
		info$COLDATE <- as.Date(info$COLDATE)

	f.clade <- f[info$COLDATE == min(info$COLDATE)]
} else {
	f.clade <- f
}

ref <- list(consensus(do.call(rbind, f.clade)))

write.fasta(c(ref, f), c("REFERENCE", names(f)), hyper.fasta.file)