library(seqinr)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "read.info.R"), chdir=T)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--hyperfasta", type='character')
args <- parse_args(op)

fasta.file <- args$fasta
info.file <- args$info
hyper.fasta.file <- args$hyperfasta

f <- read.fasta(fasta.file)
info <- read.info(info.file, names(f))

info$COLDATE <- as.Date(info$COLDATE)

ref <- list(consensus(do.call(rbind, f[info$COLDATE == min(info$COLDATE)])))

write.fasta(c(ref, f), c("REFERENCE", names(f)), hyper.fasta.file)