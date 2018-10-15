library(seqinr)
library(optparse)
library(digest)
library(parallel)

echo <- function(...) if (verbose) cat(..., '\n', sep="")

get.info <- function(name, sequence.key) {
	subset(sequence.key, enum == name)
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
op <- add_option(op, c("--threads"), type="numeric", default=1)
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
threads <- args[["threads"]]
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
echo("threads: ", threads)
echo("verbose: ", verbose)
echo()

my.apply <- if (threads > 1) function(...) mclapply(..., mc.cores=threads) else lapply

echo("Reading key...")
sequence.key <- read.csv(
	sequence.key.file,
	stringsAsFactors=F,
	col.names=c("enum", "id", "patient", "date", "type", "censored", "note"),
	colClasses=c('character', 'character', 'character', 'Date', 'character', 'numeric', "character")
)

echo("Reading csv...")
csv.files <- dir(align.csv.folder)
align.csv <- lapply(
	file.path(align.csv.folder, csv.files),
	read.csv,
	stringsAsFactors=F,
	col.names=c("refnames", "qcut", "rank", "count", "offset", "seq"),
	colClasses=c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'character')
)

echo("Parsing csv...")
name <- gsub(".+[- ](.+-V3.+).-(V3LOOP.+)_align.csv", "\\1\\2", csv.files, perl=T)
run <- gsub(".+[- ].+-V3.+(.)-V3LOOP.+_align.csv", "\\1", csv.files, perl=T)
name <- gsub("-.+", "", name)
u.name <- unique(name)

align.csv <- lapply(
	1:length(align.csv),
	function(i) {
		align.csv[[i]]$run <- run[i]
		align.csv[[i]]$file <- csv.files[i]
		align.csv[[i]]
	}
)

# combine replicate runs
align.csv.all <- my.apply(
	u.name,
	function(x) {
		echo(x)
		
		same.csv <- subset(
			as.data.frame(do.call(rbind, align.csv[name == x])),
			qcut >= q.cutoff
		)
		csv.filter <- split(same.csv, tolower(same.csv$seq))
		
		total.count <- sum(same.csv$count)
		
		as.data.frame(
			cbind(
				do.call(
					rbind,
					lapply(
						csv.filter,
						function(y) {
							with(
								y,
								data.frame(
									name=x,
									file=paste0(file, collapse=":"),
									refnames=refnames[1],
									qcut=min(qcut),
									rank=paste0(rank, collapse='-'),
									count=sum(count),
									offset=offset[1],
									seq=seq[1],
									run=paste0(run, collapse='-'),
									percent=sum(count) / total.count
								)
							)
						}
					)
				),
				get.info(x, sequence.key),
				row.names = NULL
			)
		)		
	}
)

echo("Filtering...")
align.csv.filter <- do.call(
	rbind,
	my.apply(
		align.csv.all,
		function(x) {
			d <- subset(
				x,
				count > count.cutoff & 
					count > sum(count) * percent.cutoff * 0.01
			)
			if (nrow(d) > 0 && !is.na(top))
				d[1:min(c(nrow(d), top)), ]
			else
				d
		}
	)
)

align.fasta <- strsplit(as.character(align.csv.filter$seq), "")

align.csv.filter$full.name <- with(
	align.csv.filter,
	sprintf(
		"%s_%s_%4f_%d",
		name,
		unname(substr(sapply(paste(run, rank, sep="_"), digest), 1, 6)),
		percent,
		date
	)
)
names(align.fasta) <- align.csv.filter$full.name

echo("Writing fasta...")
write.fasta(align.fasta, names(align.fasta), align.fasta.file)

echo("Writing info...")
write.csv(align.csv.filter, info.file, row.names=F)

if (verbose)
	warnings()

echo("Done")