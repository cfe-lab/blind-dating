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

compare.dates.rmse <- function(data, ori.data) {
	m <- match(data$label, ori.data$label)
	dna <- which(data$censored == 1)
	ori.dna <- m[dna]
	
	sqrt(sum((data[dna, 'est.date'] - ori.data[ori.dna, 'est.date'])^2 / length(dna)))
}

compare.dates.cor <- function(data, ori.data) {
	m <- match(data$label, ori.data$label)
	dna <- which(data$censored == 1)
	ori.dna <- m[dna]
	
	cor(data[dna, 'est.date'], ori.data[ori.dna, 'est.date'])
}

avg.var <- function(data) {
	dna <- data[[1]][which(data[[1]]$censored == 1), 'label']
		
	est.date <- do.call(rbind, lapply(data, function(x) x[match(dna, x$label), 'est.date']))
		
	mean(apply(est.date, 1, sd))
}

avg.good.var <- function(data, stats) {
	dna <- data[[1]][which(data[[1]]$censored == 1), 'label']
	good.data <- which(stats$null.AIC - stats$AIC > 10 & stats$Model.Slope > 0)
		
	if (length(good.data) > 0) {
		est.date <- do.call(rbind, lapply(data[good.data], function(x) x[match(dna, x$label), 'est.date']))
		
		mean(apply(est.date, 1, sd))
	} else
		0
}


args <- commandArgs(T)

stats.dir <- args[1]
ori.data.file <- args[2]
output.stats <- args[3]

stats <- do.call(rbind, lapply(dir(stats.dir, "*stats*", full.names=T), read.csv, stringsAsFactors=F))
data <- lapply(dir(stats.dir, "*data*", full.names=T), read.csv, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff", "ci.low", "ci.high"), stringsAsFactors=F)

ori.data <- read.csv(ori.data.file, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff", "ci.low", "ci.high"), stringsAsFactors=F)

rmses <- unlist(lapply(data, compare.dates, ori.data, ))
cors <- unlist(lapply(data, compare.dates, ori.data, cor, method='pearson'))
cors.s <- unlist(lapply(data, compare.dates, ori.data, cor, method='spearman'))
avg.vars <- avg.var(data)
avg.good.vars <- avg.good.var(data, stats)

stats <- as.data.frame(cbind(stats, full.rmse=rmses, full.cor=cors, avg.var=avg.vars, avg.good.var=avg.good.vars))
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
	"Training RMSE",
	"Censored RMSD",
	"Original RMSE",
	"Original Correlation (Pearson)",
	"Original Correlation (Spearman)",
	"Average Estimated Deviation",
	"Average Estimated Deviation (Convergent)"
)
write.table(stats, output.stats, row.names=F, col.names=stats.col.names, sep=",")