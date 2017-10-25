#!/usr/bin/Rscript

library(seqinr)

args.all <- commandArgs(F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

source(file.path(source.dir, 'read.info.R'), chdir=T)

get.mode <- function(x, multi=NA) {
	tab <- table(x)
	m <- names(tab)[tab == max(tab)]
	if (length(m) > 1)
		multi
	else
		m
}

args <- commandArgs(trailing=T)

fasta.file <- args[1]
info.file <- args[2]

hyper.fasta.file <- gsub(".fasta", ".hyper.fasta", fasta.file)

fasta <- read.fasta(fasta.file)
info <- read.info(info.file, names(fasta))

info$COLDATE <- as.Date(info$COLDATE)
first <- fasta[info$COLDATE == min(info$COLDATE)]

f.hyper <- list(REF=sapply(1:length(first[[1]]), function(i) get.mode(sapply(first, function(x) x[[i]]), multi='N')))

write.fasta(c(f.hyper, fasta), c("REF", names(fasta)), hyper.fasta.file)