library(optparse)

op <- OptionParser()
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--data", type='character')
args <- parse_args(op)

info.file <- args$info
data.file <- args$data

info <- read.csv(info.file, stringsAsFactors=F)
data <- read.csv(data.file, stringsAsFactors=F)

info <- subset(info, CENSORED > 0 & DUPLICATE %in% data$ID)
info.split <- split(info, info$TYPE)
info.split <- lapply(info.split, function(x) {
	x[sapply(unique(x$DUPLICATE), function(y) which(x$DUPLICATE == y)[1]), ]
})
info.uniq <- do.call(rbind, info.split)

info.uniq$est.date <- data[match(info.uniq$DUPLICATE, data$ID), "Estimated.Date"]

sup <- lapply(names(info.split), function(x) {
	res <- t.test(subset(info.uniq, TYPE == x)$est.date, subset(info.uniq, TYPE != x)$est.date)
	
	cat(x, ": ", res$statistic, " (p = ", res$p.value, ")\n", sep="")
})