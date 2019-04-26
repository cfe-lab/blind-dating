library(seqinr)
library(optparse)

op <- OptionParser()
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--types", type='character')
op <- add_option(op, "--dates", type='character')
args <- parse_args(op)

info.file <- args$info
fasta.file <- args$fasta
data.file <- args$data
types.file <- args$types
dates.file <- args$dates

info <- read.csv(info.file)
f <- read.fasta(fasta.file)
data <- read.csv(data.file)

info <- info[match(names(f), info$FULLSEQID), ]
row.names(info) <- NULL

m <- merge(data, info, by.x='ID', by.y='DUPLICATE', all=T)
m <- m[match(names(f), m$FULLSEQID), ]

write.table(
	with(m, data.frame(traits=FULLSEQID, subset=TYPE)),
	types.file,
	row.names=F,
	sep="\t",
	quote=F
)
write.table(
	with(m, data.frame(name=FULLSEQID, date=Estimated.Date)),
	dates.file,
	row.names=F,
	col.names=F,
	sep="\t",
	quote=F
)
				   	   