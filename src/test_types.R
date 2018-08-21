# test_types.R
# Plots the distributions of ages of the different CD4 subsets and tests for differences between the distributions

library(optparse)
library(ggplot2)
library(ggpubr)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

get.n <- function(x) data.frame(y=text.height, label=sprintf("n = %d", length(x)))

do.test <- function(x, y, name) {
	res <- wilcox.test(x$est.date, y$est.date)
	res <- list(stat=res$statistic[[1]], p=res$p.value)
	names(res) <- paste(name, names(res), sep=".")
		
	res
}

op <- OptionParser()
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--out", type='character')
op <- add_option(op, "--plot", type='character')
op <- add_option(op, "--textheight", type='character', default="2016-01-01")
op <- add_option(op, "--bymonth", type='logical', action='store_true', default=F)
args <- parse_args(op)

info.file <- args$info
data.file <- args$data
out.file <- args$out
plot.file <- args$plot
by.month <- args$bymonth
text.height <- as.Date(args$textheight)

info <- read.csv(info.file, stringsAsFactors=F)
data <- read.csv(data.file, stringsAsFactors=F)

info <- subset(info, CENSORED > 0 & DUPLICATE %in% data$ID)
info.split <- split(info, info$TYPE)
info.split <- lapply(info.split, function(x) {
	x[sapply(unique(x$DUPLICATE), function(y) which(x$DUPLICATE == y)[1]), ]
})
info.uniq <- as.data.frame(do.call(rbind, info.split))

info.uniq$est.date <- data[match(info.uniq$DUPLICATE, data$ID), "Estimated.Date"]

test.res <- as.data.frame(do.call(rbind, lapply(c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"), function(y) 
	unlist(lapply(c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"), function(x) {
		if (x == y) {
			res <- list(stat=NA, p=NA)	
			names(res) <- paste(x, names(res), sep=".")
			res
		} else {
			do.test(subset(info.uniq, TYPE == x), subset(info.uniq, TYPE == y), x)
		}
	}))
)), stringsAsFactors=F)

res.for.plot <- unname(test.res[, c(2, 4, 6, 8)])
res.for.plot <- data.frame(
	p=res.for.plot[upper.tri(res.for.plot, diag=F)],
	x=c("CD4-NAIVE", "CD4-NAIVE", "CD4-CM", "CD4-NAIVE", "CD4-CM", "CD4-TM"),
	y=c("CD4-CM", rep("CD4-TM", 2), rep("CD4-EM", 3)),
	stringsAsFactors=F
)
my.comparisons <- lapply(which(res.for.plot$p < 0.05), function(r) unlist(unname(res.for.plot[r, c("x", "y")])))

row.names(test.res) <- c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM")

write.table(test.res, file=out.file, na="NA", col.names=names(test.res), row.names=row.names(test.res), sep=",")

pdf(plot.file)
ggplot(
	info.uniq,
	aes(
		x=factor(TYPE, c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM")),
		y=as.Date(est.date, origin='1970-01-01')
	)
) +
	geom_boxplot(
		aes(colour=factor(TYPE, c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"))),
		varwidth=T,
		outlier.alpha=0,
		fill=NA,
		show.legend=F
	) +
	geom_jitter(
		aes(colour=factor(TYPE, c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"))),
		height=0,
		width=0.25
	) +
	stat_compare_means(comparisons=my.comparisons, label='p.signif') +
	stat_compare_means(label.y=text.height, size=5) +
	scale_y_date(
		name="Estimate Integration Date",
		labels=if (by.month) function(x) as.character(x, "%b %Y") else waiver()
	) +
	scale_x_discrete(
		name="Cell Type",
		breaks=c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"),
		limits=c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"),
		labels=c("N", "CM", "TM", "EM")
	) +
	scale_colour_manual(
		name="Cell Type",
		breaks=c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"),
		limits=c("CD4-NAIVE", "CD4-CM", "CD4-TM", "CD4-EM"),
		values=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
	) +
	theme_classic() +
	theme(text=element_text(size=20), axis.text = element_text(colour='black'), legend.position='none')
dev.off()