library(ape)

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

source(file.path(source.dir, 'read.info.R'), chdir=T)

args <- commandArgs(trailingOnly = T)

tree.file <- args[1]
info.file <- args[2]
pat.id <- args[3]
use.date <- if (length(args) >= 4) as.integer(args[4]) else 1

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")

tree <- read.tree(tree.file)
info <- read.info(info.file, tree$tip.label)
n <- length(tree$tip.label)

data <- data.frame(label=tree$tip.label, type=info$TYPE, censored=info$CENSORED, date=if (use.date == 1) as.numeric(as.Date(info$COLDATE)) else info$COLDATE, dist=node.depth.edgelength(tree)[1:n], stringsAsFactors=F)

g <- lm(dist ~ date, data=data, subset=censored == 0)
g.null <- lm(dist ~ 1, data=data, subset=censored == 0)

a <- coef(g)[[1]]
b <- coef(g)[[2]]
p <- predict(g, data, interval="confidence")

data <- as.data.frame(cbind(data, est.date=data$dist/b-a/b, date.diff=data$dist/b-a/b-data$date, cf.low=p[,2], cf.high=p[,3]))
write.table(data, data.file, col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Estimated Date", "Date Difference", "CI low", "CI high"), row.names=F, sep=",")

stats <- data.frame(
	pat=pat.id,
	samples.rna=sum(data$censored == 0),
	samples.dna=sum(data$censored == 1),
	samples.total=n,
	rna.time.points=length(unique(data[data$censored == 0, 'date'])),
	dna.time.points=length(unique(data[data$censored == 1, 'date'])),
	total.time.points=length(unique(data$date)),
	min.rna.time.point=min(data[data$censored == 0, 'date']),
	max.rna.time.point=max(data[data$censored == 0, 'date']),
	min.dna.time.point=min(data[data$censored == 1, 'date']),
	max.dna.time.point=max(data[data$censored == 1, 'date']),
	min.time.point=min(data$date),
	max.time.point=max(data$date),
	AIC=AIC(g),
	null.AIC=AIC(g.null),
	p=1-pchisq(AIC(g.null) - AIC(g) + 2, 1),
	a=a,
	b=b,
	error=0,
	mu=-a/b,
	train.RMSE=sqrt(sum(data$date.diff[data$censored == 0]^2)/sum(data$censored == 0)),
	cens.RMSD=sqrt(sum(data$date.diff[data$censored == 1]^2)/sum(data$censored == 1)),
	train.MAE=sum(abs(data$date.diff[data$censored == 0]))/sum(data$censored == 0),
	cens.MAE=sum(abs(data$date.diff[data$censored == 1]))/sum(data$censored == 1)
)
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
	"Training MAE",
	"Censored MAE"
)
write.table(stats, stats.file, col.names=stats.col.names, row.names=F, sep=",")
