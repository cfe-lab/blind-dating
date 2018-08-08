library(optparse)
library(ggplot2)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

get.n <- function(x) data.frame(y=max(x)+offset, label=sprintf("n = %d", length(x)))

op <- OptionParser()
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--plot", type='character')
op <- add_option(op, "--offset", type='numeric', default=100)
args <- parse_args(op)

info.file <- args$info
data.file <- args$data
plot.file <- args$plot

info <- read.csv(info.file, stringsAsFactors=F)
data <- read.csv(data.file, stringsAsFactors=F)

info <- subset(info, CENSORED > 0 & DUPLICATE %in% data$ID)
info.split <- split(info, info$TYPE)
info.split <- lapply(info.split, function(x) {
	x[sapply(unique(x$DUPLICATE), function(y) which(x$DUPLICATE == y)[1]), ]
})
info.uniq <- as.data.frame(do.call(rbind, info.split))

info.uniq$est.date <- data[match(info.uniq$DUPLICATE, data$ID), "Estimated.Date"]

pdf(plot.file)
ggplot(info.uniq, aes(x=TYPE, y=as.Date(est.date, origin='1970-01-01'), fill=TYPE)) +
	geom_violin(show.legend=F) +
	stat_summary(fun.data=get.n, geom='text')
	scale_y_date(name="Estimate Integration Date") +
	scale_x_discrete(name="Cell Type", breaks=c("CD4-CM", "CD4-EM", "CD4-TM", "CD4-NAIVE"), limits=c("CD4-CM", "CD4-EM", "CD4-TM", "CD4-NAIVE"), labels=c("Central Memory", "Effector Memory", "Transitory Memory", "Naive")) +
	theme_classic() +
	theme(text=element_text(size=15))
dev.off()

sup <- lapply(names(info.split), function(x) {
	res <- t.test(subset(info.uniq, TYPE == x)$est.date, subset(info.uniq, TYPE != x)$est.date)
	
	cat(x, ": ", res$statistic, " (p = ", res$p.value, ")\n", sep="")
})