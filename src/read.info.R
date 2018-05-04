read.info <- function(f, label) {
	d <- read.csv(f, stringsAsFactors=F)
	m <- match(label, d$FULLSEQID)
	d.f <- d[m, ]
	rownames(d.f) <- 1:nrow(d.f)
	d.f
}