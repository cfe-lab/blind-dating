library(pegas)
library(optparse)

op <- OptionParser()
#op <- add_option(op, "--name", type='character')
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--output", type='character')
args <- parse_args(op)

#name <- args$name
fasta.file <- args$fasta
info.file <- args$info
output.file <- args$output

f <- read.FASTA(fasta.file)
info <- read.csv(info.file)

info <- info[match(labels(f), info$FULLSEQID), ]

dist <- dist.dna(f)
types <- factor(info$TYPE)

g <- amova(dist ~ types)

write(capture.output(print(g)), output.file)