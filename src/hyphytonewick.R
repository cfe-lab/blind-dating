library(optparse)
library(ape)
library(magrittr)
library(circlize)

op <- OptionParser()
op <- add_option(op, "--hyphy", type='character')
op <- add_option(op, "--plot", type='character')
args <- parse_args(op)

hyphy.file <- args$hyphy
plot.file <- args$plot

hyphy <- readLines(hyphy.file)

tree <- read.tree(text=gsub(".+=", "", hyphy[1]))

hyphy.states <- hyphy[c(-1, -2)] %>%
	gsub("\\{\\{", "list(", .) %>%
	gsub("\\}\\}", ")", .) %>%
	gsub("\\[", "\\[\\[", .) %>%
	gsub("\\]", "\\]\\]", .) %>%
	gsub("^_", "", .)
	
branchClasses <- list()
eval(parse(text=hyphy.states))

trans <- branchClasses[grep("-->", names(branchClasses))]
trans %<>% sapply(length)
trans.from <- gsub("(.+) --> (.+)", "\\1", names(trans))
trans.to <-gsub("(.+) --> (.+)", "\\2", names(trans))

subsets <- c("N", "CM", "TM", "EM")
cols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

pdf(plot.file, height=7, width=7)
circos.initialize(factor=subsets, xlim=c(0, 6))
circos.track(ylim=c(0, 1), bg.col=cols)

for (i in 1:length(trans)) {
	w <- trans[i] / 24 * c(-1, 1)
	id.diff <- (which(subsets == trans.to[i]) - which(subsets == trans.from[i]) + 4) %% 4
	from.x <- 7.5 - 2 * id.diff
	to.x <- -1.5 +  2 * id.diff
	circos.link(
		trans.from[i],
		from.x + w,
		trans.to[i],
		to.x + w,
		col=cols[which(subsets == trans.from[i])],
		arr.type='triangle'
	)
}

dev.off()