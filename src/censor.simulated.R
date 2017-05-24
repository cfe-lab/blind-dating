

args <- commandArgs(trailingOnly=T)

info.path <- args[1]
censor <- as.integer(args[2])
seed <- as.numeric(args[3])

set.seed(seed)

files <- dir(info.path)

info <- lapply(paste(info.path, files, sep="/"), read.csv, stringsAsFactors=T)

info <- lapply(info, function(x) {x[sample(nrow(x), censor), 'CENSORED'] <- 1; x})

suppress <- lapply(1:length(info), function(i) write.csv(info[[i]], paste0("info/", files[i]), row.names=F))