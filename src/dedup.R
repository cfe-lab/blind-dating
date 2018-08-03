library(seqinr)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, 'read.info.R'), chdir=T)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--newfasta", type='character')
op <- add_option(op, "--newinfo", type='character')
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--weight", type='character', default="COLDATE")
op <- add_option(op, "--lowest", type="logical", action='store_true', default=T)
op <- add_option(op, "--highest", type="logical", dest="lowest", action='store_false')
args <- parse_args(op)

fasta.file <- args$fasta
info.file <- args$info
new.fasta.file <- args$newfasta
new.info.file <- args$newinfo
use.real <- args$real
weight <- args$weight
use.lowest <- args$lowest

f <- read.fasta(fasta.file)

same <- lapply(f, function(x) which(unlist(lapply(f, function(y) all(x == y)))))
same <- same[unlist(lapply(same, length)) > 1]
same <- unique(lapply(same, function(x) {names(x) <- NULL; x}))

info <- read.info(info.file, names(f))

o <- order(if (use.real) as.numeric(as.Date(info[, weight])) else info[, weight], decreasing=!use.lowest)
info <- info[o, ]
o <-  order(info$CENSORED)
info <- info[o, ]
o <- match(info$FULLSEQID, names(f))
info <- info[match(names(f), info$FULLSEQID), ]

same <- lapply(same, function(x) o[o %in% x])

dont.keep <- unlist(lapply(same, function(x) x[-1]))
all.rows <- unlist(same)
nall.rows <- if (is.null(all.rows)) 1:nrow(info) else -all.rows
dups <- do.call(c, lapply(same, function(x) lapply(1:length(x), function(y) info$FULLSEQID[x[-y]])))
rep.dup <- unlist(lapply(same, function(x) rep(info$FULLSEQID[x[1]], length(x))))

info$KEPT <- 1
info$KEPT[dont.keep] <- 0
info$DUPLICATES[all.rows] <- unlist(lapply(dups, paste, collapse=";"))
info$DUPLICATE[all.rows] <- rep.dup
info$DUPLICATE[nall.rows] <- info$FULLSEQID[nall.rows]

write.csv(info, new.info.file, row.names=F, na="")

f <- f[info$KEPT == 1]
write.fasta(f, names(f), new.fasta.file)