library(seqinr)
nef.start <- 8149
nef.end <- 8795

f <- read.fasta("seqs/mullens-aligned-trimmed.fasta")
#!/usr/local/bin/Rscript

df <- do.call(rbind, f)

cons <- apply(df[grep("Oct-1998", row.names(df)), ], 2, get.mode, dup="-")
filter <- cons != "-"

all.ids <- (1:length(cons))[filter]
nef.ids <- all.ids[all.ids >= nef.start & all.ids <= nef.end]
notnef.ids <- all.ids[all.ids < nef.start | all.ids > nef.end]

all.diffs <- colSums(t(df[, all.ids]) != cons[all.ids]) / length(all.ids)
nef.diffs <- colSums(t(df[, nef.ids]) != cons[nef.ids]) / length(nef.ids)
notnef.diffs <- colSums(t(df[, notnef.ids]) != cons[notnef.ids]) / length(notnef.ids)


f.cleaned <- read.fasta("seqs/mullens-aligned-stripped.cleaned.fasta")
df <- df[row.names(df) %in% names(f.cleaned), ]

filter <- apply(df, 2, function(x) !any(x == '-'))

cons <- apply(df[grep("Oct-1998", row.names(df)), ], 2, get.mode, dup="-")
filter <- filter | cons != "-"

all.ids <- (1:length(cons))[filter]
nef.ids <- all.ids[all.ids >= nef.start & all.ids <= nef.end]
notnef.ids <- all.ids[all.ids < nef.start | all.ids > nef.end]

all.diffs <- colSums(t(df[, all.ids]) != cons[all.ids]) / length(all.ids)
nef.diffs <- colSums(t(df[, nef.ids]) != cons[nef.ids]) / length(nef.ids)
notnef.diffs <- colSums(t(df[, notnef.ids]) != cons[notnef.ids]) / length(notnef.ids)

cor(all.diffs, nef.diffs)
cor(notnef.diffs, nef.diffs)
