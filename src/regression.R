library(ape)
library(optparse)
library(chemCal)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "read.info.R"), chdir=T)

get.val <- function(x, default) if (is.null(x)) default else x

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
op <- add_option(op, "--real", type='logical', action='store_true')
op <- add_option(op, "--usedups", type='logical', action='store_true')
op <- add_option(op, "--cutoff", type="character")
op <- add_option(op, "--data", type='character')
op <- add_option(op, "--stats", type='character')
op <- add_option(op, "--regression", type='character')
op <- add_option(op, "--freqweights", type='logical', action='store_true')
op <- add_option(op, "--settings", type='character', default=NA)
args <- parse_args(op)

settings.file <- args$settings
if (!is.na(settings.file)) {
	settings <- readLines(settings.file)
	settings.filter <- unlist(lapply(op@options, function(x) settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]))
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}


tree.file <- args$tree
info.file <- args$info
pat.id <- get.val(args$patid, NA)
use.date <- get.val(args$real, F)
use.all <- get.val(args$usedups, F)
cutoff <- get.val(args$cutoff, NA)
data.file <- get.val(args$data, NA)
stats.file <- get.val(args$stats, NA)
regression.file <- get.val(args$regression, NA)
freq.weights <- get.val(args$freqweights, F)

set.seed(args$seed)

if (is.na(data.file)) data.file <- paste0("stats/", pat.id, ".data.csv")
if (is.na(stats.file)) stats.file <- paste0("stats/", pat.id, ".stats.csv")
if (is.na(regression.file)) regression.file <- paste0("stats/", pat.id, ".regression.rds")

tree <- read.tree(tree.file)
info <- if (use.all) read.csv(info.file, stringsAsFactors=F) else read.info(info.file, tree$tip.label)
n <- length(tree$tip.label)

data <- data.frame(label=info$FULLSEQID, type=info$TYPE, censored=info$CENSORED, date=if (use.date == 1) as.numeric(as.Date(info$COLDATE)) else as.numeric(info$COLDATE), dist=node.depth.edgelength(tree)[if (use.all) match(info$DUPLICATE, tree$tip.label) else 1:n], stringsAsFactors=F)

weights <- if (freq.weights) {
	info$COUNT[if (use.all) match(info$DUPLICATE, tree$tip.label) else 1:n, ]
} else {
	NULL
}

g <- lm(dist ~ date, data=data, weights=weights, subset=censored == 0)
g.null <- lm(dist ~ 1, data=data, weights=weights, subset=censored == 0)

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

ci <- inverse.predict(g, 0)

bin.test <- tryCatch(binom.test(sum(data[data$censored == 1, "date.diff"] < 0), sum(data$censored == 1), alternative='two.sided'), error=function(x) list(p.value=NA, estimate=NA))
T.test <- tryCatch(t.test(data[data$censored == 1, "date.diff"], alternative='two.sided'), error=function(x) list(p.value=NA, estimate=NA))

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
	rsquared=summary(g)$r.squared,
	root.date=-a / b,
	ci.lwr=ci[['Confidence Limits']][1],
	ci.upr=ci[['Confidence Limits']][2],
	fit=as.numeric(AIC(g.null) - AIC(g) > 10 && (-a / b < cutoff || (!is.na(ci[['Confidence Limits']][1]) && ci[['Confidence Limits']][1] < cutoff)) && b > 0),
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
	"Rsquared",
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
