library(ape)
library(optparse)

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

concord <- function(x, y) {
	mu.x <- sum(x) / length(x)
	mu.y <- sum(y) / length(y)
	s.x <- sum((x - mu.x)^2) / length(x) 
	s.y <-  sum((y - mu.y)^2) / length(y)
	s.xy <- sum((x - mu.x) * (y - mu.y)) / length(y)
	
	2 * s.xy / (s.x + s.y + (mu.x - mu.y)^2)
}

op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--patid", type='character')
op <- add_option(op, "--real", type='logical', action='store_true', default=F)
op <- add_option(op, "--usedups", type='logical', action='store_true', default=F)
op <- add_option(op, "--seed", type="numeric", default=1989)
op <- add_option(op, "--cutoff", type="character", default=NA)
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info
pat.id <- args$patid
use.date <- args$real
use.all <- args$usedups
cutoff <- args$cutoff

set.seed(args$seed)

data.file <- paste0("stats/", pat.id, ".data.csv")
stats.file <- paste0("stats/", pat.id, ".stats.csv")
regression.file <- paste0("stats/", pat.id, ".regression.rds")

tree <- read.tree(tree.file)
info <- if (use.all) read.csv(info.file, stringsAsFactors=F) else read.info(info.file, tree$tip.label)
n <- length(tree$tip.label)

data <- data.frame(label=info$FULLSEQID, type=info$TYPE, censored=info$CENSORED, date=if (use.date == 1) as.numeric(as.Date(info$COLDATE)) else as.numeric(info$COLDATE), dist=node.depth.edgelength(tree)[if (use.all) match(info$DUPLICATE, tree$tip.label) else 1:n], stringsAsFactors=F)

g <- lm(dist ~ date, data=data, subset=censored == 0)
g.null <- lm(dist ~ 1, data=data, subset=censored == 0)

a <- coef(g)[[1]]
b <- coef(g)[[2]]

data <- as.data.frame(cbind(data, est.date=data$dist/b-a/b, date.diff=data$dist/b-a/b-data$date))
write.table(data, data.file, col.names=c("ID", "Type", "Censored", "Collection Date", "Divergence", "Estimated Date", "Date Difference"), row.names=F, sep=",")


cutoff  <- if (is.na(cutoff)) {
	min(data$date)
} else {
	if (use.date)
		as.numeric(as.Date(cutoff))
	else
		as.numeric(cutoff)
}

ci <- predict(g, newdata=data.frame(date=-a / b), interval='confidence')
rownames(ci) <- NULL
ci <- as.data.frame(ci)	
ci$lwr <- ci$lwr / b - a / b
ci$upr <- ci$upr / b - a / b

bin.test=binom.test(sum(data[data$censored == 1, "date.diff"] < 0), sum(data$censored == 1), alternative='two.sided')
T.test <- t.test(data[data$censored == 1, "date.diff"], alternative='two.sided')

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
	root.date=-a / b,
	ci.lwr=ci$lwr,
	ci.upr=ci$upr,
	fit=as.numeric(AIC(g.null) - AIC(g) > 10 && ci$lwr < cutoff && b > 0),
	train.RMSE=sqrt(sum(data$date.diff[data$censored == 0]^2)/sum(data$censored == 0)),
	cens.RMSD=sqrt(sum(data$date.diff[data$censored == 1]^2)/sum(data$censored == 1)),
	cens.RMSD=sqrt(sum(data$date.diff^2)/nrow(data)),
	train.MAE=sum(abs(data$date.diff[data$censored == 0]))/sum(data$censored == 0),
	cens.MAE=sum(abs(data$date.diff[data$censored == 1]))/sum(data$censored == 1),
	tot.MAE=sum(abs(data$date.diff))/nrow(data),
	tot.concord=concord(data$date, data$est.date),
	bin.test=bin.test$p.value,
	bin.mean=bin.test$estimate,
	t.test=T.test$p.value,
	t.mean=T.test$estimate
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
	"ERD CI low",
	"ERD CI high",
	"Model Fit",
	"Training RMSE",
	"Censored RMSD",
	"Total RMSD",
	"Training MAE",
	"Censored MAE",
	"Total MAE",
	"Total Concordance",
	"Bin Test (p)",
	"Bin Test (mean)",
	"T-Test (p)",
	"T-Test (mean)"
)
write.table(stats, stats.file, col.names=stats.col.names, row.names=F, sep=",")

saveRDS(g, regression.file)
