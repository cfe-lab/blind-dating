library(optparse)
library(ape)
library(magrittr)

op <- OptionParser()
op <- add_option(op, "--hyphy", type='character')
op <- add_option(op, "--tree", type='character')
args <- parse_args(op)

hyphy.file <- args$hyphy
tree.file <- args$tree

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

trans <- branchClasses[1:(length(branchClasses) - 4)]
trans %<>% sapply(length)
trans.from <- gsub("(.+) --> (.+)", "\\1", names(trans))
trans.to <-gsub("(.+) --> (.+)", "\\2", names(trans))