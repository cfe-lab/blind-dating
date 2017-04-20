

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

compare.dates <- function(data, ori.data) {
	m <- match(data$label, ori.data$label)
	dna <- which(data$type == "PBMC")
	ori.dna <- m[dna]
	
	sqrt(sum((data[dna, 'est.date'] - ori.data[ori.dna, 'est.date'])^2 / sum(dna)))
}

args <- commandArgs(T)

stats.dir <- args[1]
ori.data.file <- args[2]
output.stats <- args[3]

stats <- do.call(rbind, lapply(dir(stats.dir, "*stats*", full.names=T), read.csv, stringsAsFactors=F))
data <- lapply(dir(stats.dir, "*data*", full.names=T), read.csv, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

ori.data <- read.csv(ori.data.file, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

good <- stats[, "p.value"] < 0.001

rmses <- unlist(lapply(data, compare.dates, ori.data))
stats <- as.data.frame(cbind(stats, full.rmse=rmses))
stats.col.names <- c(
	"Patient",
	"RNA Samples",
	"DNA Samples",
	"Total Samples",
	"RNA Time Points",
	"DNA Time Points",
	"Total Time Points",
	"Minimum RNA Time Point",
	"Maximum RNA Time Point",
	"Minimum DNA Time Point",
	"Maximum DNA Time Point",
	"Minimum Time Point",
	"Maximum Time Point",
	"AIC",
	"null AIC",
	"p-value",
	"Model Intercept",
	"Model Slope",
	"Estimated Mutation Rate",
	"Training RMSE",
	"Censored RMSD",
	"Full RMSE"
)
write.table(stats, output.stats, row.names=F, col.names=stats.col.names, sep=",")

good.rmses <- rmses[good]
output.data <- data.frame(good=good/length(data), mean.rmse=mean(good.rmses), median.rmse=median(good.rmses), sd.rmse=sd(good.rmses), mean.AIC.diff=mean(stats[good, "null.AIC"] - stats[good, "AIC"]))
write.csv(stats, stats.file, row.names=F)