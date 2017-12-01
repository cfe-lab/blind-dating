

args <- commandArgs(trailingOnly=T)

info.path <- args[1]
censor <- as.numeric(args[2])
seed <- as.numeric(args[3])

set.seed(seed)

info.files <- dir(info.path, full.names=T)
info <- lapply(info.files, read.csv)
info <- lapply(info, function(x) {y <- x; y[sample.int(nrow(x), censor), 'CENSORED'] <- 1; y})

suppress <- lapply(1:length(info), function(i) write.csv(info[[i]], gsub(".csv", ".cens.csv", info.files[[i]]), row.names=F))