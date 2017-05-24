library(seqinr)

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.param <- grep("f", args.all)
	source.dir <- dirname(args.all[file.param + 1])
}

source(file.path(source.dir, 'read.info.R'), chdir=T)

args <- commandArgs(trailingOnly = T)

fasta.file <- args[1]
info.file <- args[2]

f <- read.fasta(fasta.file)

same <- lapply(f, function(x) which(unlist(lapply(f, function(y) all(x == y)))))
same <- same[unlist(lapply(same, length)) > 1]
same <- unique(lapply(same, function(x) {names(x) <- NULL; x}))

info <- read.info(info.file, names(f)[names(f) != "REFERENCE"])

o <- match(info$FULLSEQID, names(f))[order(as.numeric(as.Date(info$COLDATE)))]

same <- lapply(same, function(x) o[o %in% x])

dont.keep <- unlist(lapply(same, function(x) x[-1]))
all.rows <- unlist(same)
dups <- do.call(c, lapply(same, function(x) lapply(1:length(x), function(y) info$FULLSEQID[x[-y]])))

info$KEPT <- 1
info$KEPT[dont.keep] <- 0
info$DUPLICATE[all.rows] <- unlist(lapply(dups, paste, collapse=";"))

info <- info[!is.na(info$PATID), ]
write.csv(info, info.file, row.names=F, na="")