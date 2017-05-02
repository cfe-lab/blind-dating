library(seqinr)

args <- commandArgs(trailingOnly = T)

fasta.file <- args[1]

f <- read.fasta(fasta.file)

same <- lapply(f, function(x) which(unlist(lapply(f, function(y) all(x == y)))))

same <- same[unlist(lapply(same, length)) > 1]

print(unique(same))