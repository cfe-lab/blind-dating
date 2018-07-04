library(seqinr)
library(optparse)

echo <- function(...) if (verbose) cat(..., '\n', sep="")

op <- OptionParser()
op <- add_option(op, c("-c", "--align-csv"), type='character')
op <- add_option(op, c("-f", "--align-fasta"), type='character')
op <- add_option(op, c("-q", "--q-cutoff"), type='numeric', default=15)
op <- add_option(op, c("-n", "--count-cutoff"), type="numeric", default=1)
op <- add_option(op, c("-p", "--seq-id-prefix"), type='character', default=NA)
op <- add_option(op, c("-s", "--seq-id-suffix"), type='character', default=NA)
op <- add_option(op, c("-v", "--verbose"), type='logical', action="store_true", default=F)
args <- parse_args(op)

align.csv.file <- args[["align-csv"]]
align.fasta.file <- args[["align-fasta"]]
q.cutoff <- args[["q-cutoff"]]
count.cutoff <- args[["count-cutoff"]]
seq.id.prefix <- args[["seq-id-prefix"]]
seq.id.suffix <- args[["seq-id-suffix"]]
verbose <- args[["verbose"]]

echo("extract_miseq.R")
echo("Params:")
echo("align-csv: ", align.csv.file)
echo("align-fasta: ", align.fasta.file)
echo("q-cutoff: ", q.cutoff)
echo("count-cutoff: ", count.cutoff)
echo("seq-id-prefix: ", seq.id.prefix)
echo("seq-id-suffix: ", seq.id.suffix)
echo("verbose: ", verbose)
echo()

echo("Reading csv...")
align.csv <- read.csv(align.csv.file, stringsAsFactors=F, col.names=c("refnames", "qcut", "rank", "count", "offset", "seq"), colClasses=c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'character'))

echo("Parsing csv...")
align.csv.filter <- subset(align.csv, qcut >= q.cutoff & count > count.cutoff)
align.fasta <- strsplit(align.csv.filter$seq, "")

names(align.fasta) <- paste(align.csv.filter$rank, align.csv.filter$count, sep="_")
if (!is.na(seq.id.prefix))
	names(align.fasta)  <- paste(seq.id.prefix, names(align.fasta), sep="_")
if (!is.na(seq.id.suffix))
	names(align.fasta)  <- paste(names(align.fasta), seq.id.suffix, sep="_")

echo("Writing fasta...")
write.fasta(align.fasta, names(align.fasta), align.fasta.file)

echo("Done")