read.info <- function(f, label) {
	d <- read.csv(f, stringsAsFactors=F)
	m <- match(label, d$FULLSEQID)
	d[m, ]
}