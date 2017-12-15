#!/usr/bin/Rscript

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

compare.dates <- function(data, ori.data, censor, f, ...) {
	if (censor) {
		data <- subset(data, censored > 0)
		ori.data <- subset(ori.data, censored > 0)
	}

	m <- match(data$label, ori.data$label)

	f(data[, 'est.date'], ori.data[m, 'est.date'], ...)
}

concord <- function(x, y) {
	mu.x <- sum(x) / length(x)
	mu.y <- sum(y) / length(y)
	s.x <- sum((x - mu.x)^2) / length(x) 
	s.y <-  sum((y - mu.y)^2) / length(y)
	s.xy <- sum((x - mu.x) * (y - mu.y)) / length(y)
	
	2 * s.xy / (s.x + s.y + (mu.x - mu.y)^2)
}

args <- commandArgs(T)

stats.dir <- args[1]
ori.data.file <- args[2]
output.stats <- args[3]
filt <- args[4]

files <- dir(stats.dir, filt, full.names=T)
stats <- do.call(rbind, lapply(files, read.csv, stringsAsFactors=F))
data <- lapply(gsub("stats\\.", "data.", files), read.csv, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

ori.data <- read.csv(ori.data.file, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

#rmses <- tryCatch(unlist(lapply(data, compare.dates, ori.data, F, function(x, y) sqrt(sum((x - y)^2)/length(x)))), exception=function(e) NA)
#maes <- tryCatch(unlist(lapply(data, compare.dates, ori.data, F, function(x, y) sum(abs(x - y))/length(x))), exception=function(e) NA)
#cors <- tryCatch(unlist(lapply(data, compare.dates, ori.data, F, cor, method='pearson')), exception=function(e) NA)
#cors.s <- tryCatch(unlist(lapply(data, compare.dates, ori.data, F, cor, method='spearman')), exception=function(e) NA)
#concord <- tryCatch(unlist(lapply(data, compare.dates, ori.data, F, concord)), exception=function(e) NA)
rmses.cens <- unlist(lapply(data, compare.dates, ori.data, T, function(x, y) sqrt(sum((x - y)^2)/length(x))))
maes.cens <- unlist(lapply(data, compare.dates, ori.data, T, function(x, y) sum(abs(x - y))/length(x)))
cors.cens <- unlist(lapply(data, compare.dates, ori.data, T, cor, method='pearson'))
cors.s.cens <- unlist(lapply(data, compare.dates, ori.data, T, cor, method='spearman'))
concord.cens <- unlist(lapply(data, compare.dates, ori.data, T, concord))

stats <- as.data.frame(cbind(stats, full.rmse.cens=rmses.cens, full.mae.cens=maes.cens, full.cor.cens=cors.cens, spearman.cor.cens=cors.s.cens, concord.cens=concord.cens), stringsAsFactors=F)
stats.col.names <- c(
	"Patient",
	"Training Samples",
	"Censored Samples",
	"Total Samples",
	"Training Time Points",
	"Censored Time Points",
	"Total Time Points",
	"Minimum Training Time Point",
	"Maximum Training Time Point",
	"Minimum Censored Time Point",
	"Maximum Censored Time Point",
	"Minimum Time Point",
	"Maximum Time Point",
	"AIC",
	"null AIC",
	"p-value",
	"Model Intercept",
	"Model Slope",
	"Model Error",
	"Estimated Root Date",
	"Model Fit",
	"Training RMSE",
	"Censored RMSD",
	"Total RMSD",
	"Training MAE",
	"Censored MAE",
	"Total MAE",
	"Total Concordance",
	"Bin Test",
	"Original RMSE",
	"Original MAE",
	"Original Correlation",
	"Original Correlation",
	"Original Concordance",
)
write.table(stats, output.stats, row.names=F, col.names=stats.col.names, sep=",")