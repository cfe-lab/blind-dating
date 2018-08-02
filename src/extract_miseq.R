library(seqinr)
library(optparse)

echo <- function(...) if (verbose) cat(..., '\n', sep="")

get.info <- function(name, sequence.key) {
	short.name <- gsub("-.+", "", name)
	subset(sequence.key, enum == short.name)
}

op <- OptionParser()
op <- add_option(op, c("-c", "--align-csv-folder"), type='character')
op <- add_option(op, c("-k", "--sequence-key"), type='character')
op <- add_option(op, c("-f", "--align-fasta"), type='character')
op <- add_option(op, c("-i", "--info"), type='character')
op <- add_option(op, c("-q", "--q-cutoff"), type='numeric', default=15)
op <- add_option(op, c("-n", "--count-cutoff"), type="numeric", default=100)
op <- add_option(op, c("-p", "--percent-cutoff"), type="numeric", default=1)
op <- add_option(op, c("-t", "--top"), type="numeric", default=NA)
op <- add_option(op, c("-v", "--verbose"), type='logical', action="store_true", default=F)
args <- parse_args(op)

align.csv.folder <- args[["align-csv-folder"]]
sequence.key.file <- args[["sequence-key"]]
align.fasta.file <- args[["align-fasta"]]
info.file <- args[["info"]]
q.cutoff <- args[["q-cutoff"]]
count.cutoff <- args[["count-cutoff"]]
percent.cutoff <- args[["percent-cutoff"]]
top <- args[["top"]]
verbose <- args[["verbose"]]

echo("extract_miseq.R")
echo("Params:")
echo("align-csv-folder: ", align.csv.folder)
echo("sequence-key: ", sequence.key.file)
echo("align-fasta: ", align.fasta.file)
echo("info: ", info.file)
echo("q-cutoff: ", q.cutoff)
echo("count-cutoff: ", count.cutoff)
echo("percent-cutoff: ", percent.cutoff)
echo("top: ", top)
echo("verbose: ", verbose)
echo()

echo("Reading key...")
sequence.key <- read.csv(sequence.key.file,  stringsAsFactors=F, col.names=c("enum", "id", "patient", "date", "type", "censored", "note"), colClasses=c('character', 'character', 'character', 'Date', 'character', 'numeric', "character"))

echo("Reading csv...")
csv.files <- dir(align.csv.folder)
align.csv <- lapply(file.path(align.csv.folder, csv.files), read.csv, stringsAsFactors=F, col.names=c("refnames", "qcut", "rank", "count", "offset", "seq"), colClasses=c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'character'))

echo("Parsing csv...")
name <- gsub(".+[- ](.+-V3.+V3LOOP.+)_align.csv", "\\1", csv.files, perl=T)
align.csv.all <- lapply(1:length(csv.files), function(i) 
	cbind(name=name[i], align.csv[[i]], get.info(name[i], sequence.key), percent=with(align.csv[[i]], count / sum(count)))
)

align.csv.filter <- do.call(rbind, lapply(align.csv.all, function(x) {
	d <- subset(x, qcut >= q.cutoff & count > count.cutoff & count > sum(count) * percent.cutoff * 0.01)
	if (nrow(d) > 0 && !is.na(top))
		d[1:min(c(nrow(d), top)), ]
	else
		d
}))

align.fasta <- strsplit(align.csv.filter$seq, "")

align.csv.filter$full.name <- with(align.csv.filter, sprintf("%s_%d_%4f_%d", name, rank, percent, date))
names(align.fasta) <- align.csv.filter$full.name

warnings()

echo("Writing fasta...")
write.fasta(align.fasta, names(align.fasta), align.fasta.file)

echo("Writing info...")
write.csv(align.csv.filter, info.file, row.names=F)

echo("Done")