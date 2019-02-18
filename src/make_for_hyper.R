library(seqinr)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "read.info.R"), chdir=T)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--hyperfasta", type='character')
op <- add_option(op, "--rename", type='character', default=NA)
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--useall", type='logical', action='store_true', default=F)
args <- parse_args(op)

fasta.file <- args$fasta
info.file <- args$info
hyper.fasta.file <- args$hyperfasta
rename.file <- args$rename
use.date <- args$real
use.all <- args$useall

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

if (is.na(rename.file)) {
	f.names <- names(f)
} else {
	f.names <- paste0("S", 1:length(f))
	
	write.csv(data.frame(label=names(f), newlabel=f.names), rename.file, row.names=F)
}

write.fasta(c(ref, f), c("REFERENCE", f.names), hyper.fasta.file)