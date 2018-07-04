#!/usr/bin Rscript

filter.and.move <- function(seq.row, files, input.folder, output.folder) {
	filt <- grepl(seq.row[1], files)
	
	lapply(files[filt], function(x) {cat(x, '\n', sep=""); file.copy(file.path(input.folder, x), file.path(output.folder, seq.row[3], x))})
}

args <- commandArgs(trailingOnly = T)

align.csv.folder <- args[1] # %FOLDER%/align_csv
sequence.key.file <- args[2] # SequenceKey.csv
output.folder <- args[3] # %FOLDER%

sequence.key <- read.csv(sequence.key.file, stringsAsFactors=F)

align.csv.files <- dir(align.csv.folder, ".+csv")

sup <- apply(sequence.key, 1, filter.and.move, align.csv.files, align.csv.folder, output.folder)