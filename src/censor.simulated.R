

args <- commandArgs(trailingOnly=T)

info.path <- args[1]
count <- as.integer(args[2])
censor <- as.integer(args[3])
seed <- as.numeric(args[3])

set.seed(seed)

info.ori <- read.csv(info.path, stringsAsFactors=T)
info.ori <- info.ori[with(info.ori, CENSORED == 0 & KEPT == 1), ]

info <- lapply(1:count, function(x) {y <- info.ori; y[sample.int(nrow(info.ori), censor, replace=F), 'CENSORED'] <- 1; y})

suppress <- lapply(1:length(info), function(i) write.csv(info[[i]], paste0("info/SIM_", i, ".csv"), row.names=F))