

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
	dna <- which(data$censored == 0)
	ori.dna <- m[dna]
	
	sqrt(sum((data[dna, 'est.date'] - ori.data[ori.dna, 'est.date'])^2 / length(dna)))
}

compare.dates.cor <- function(data, ori.data) {
	m <- match(data$label, ori.data$label)
	dna <- which(data$censored == 0)
	ori.dna <- m[dna]
	
	cor(data[dna, 'est.date'], ori.data[ori.dna, 'est.date'])
}

args <- commandArgs(T)

stats.dir <- args[1]
ori.data.file <- args[2]
output.stats <- args[3]

stats <- do.call(rbind, lapply(dir(stats.dir, "*stats*", full.names=T), read.csv, stringsAsFactors=F))
data <- lapply(dir(stats.dir, "*data*", full.names=T), read.csv, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

ori.data <- read.csv(ori.data.file, col.names=c("label", "type", "censored", "date", "dist", "est.date", "date.diff"), stringsAsFactors=F)

rmses <- unlist(lapply(data, compare.dates.rmse, ori.data))
cors <- unlist(lapply(data, compare.dates.cor, ori.data))
stats <- as.data.frame(cbind(stats, full.rmse=rmses, full.cor=cors))
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
	"Estimated Mutation Rate",
	"Training RMSE",
	"Censored RMSD",
	"Full RMSE",
	"Full Correlation" 
)
write.table(stats, output.stats, row.names=F, col.names=stats.col.names, sep=",")