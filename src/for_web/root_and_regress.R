library(ape)
library(optparse)
library(chemCal)

DATE_FMT <- "%Y-%m-%d"

# utility function
to.Date <- function(...) as.Date(..., origin = "1970-01-01")

# parse command line arguments
op <- OptionParser()
op <- add_option(op, "--runid", type = 'character')
op <- add_option(op, "--tree", type = 'character')
op <- add_option(op, "--info", type = 'character')
op <- add_option(op, "--rootedtree", type = 'character')
op <- add_option(op, "--data", type = 'character')
op <- add_option(op, "--stats", type = 'character')
args <- parse_args(op)

run.id <- args$runid
tree.file <- args$tree
info.file <- args$info
rooted.tree.file <- args$rootedtree
data.file <- args$data
stats.file <- args$stats

### TODO check arguments

# read tree and info
tree <- read.tree(tree.file)
info <- read.csv(info.file, stringsAsFactors = FALSE)

### TODO check info has the correct columns

info <- info[match(tree$tip.label, info$ID), ]
info$Date <- as.numeric(as.Date(info$Date, format = DATE_FMT))

# root tree
dates <- info$Date
dates[info$Query == 1] <- NA

rooted.tree <- rtt(tree, dates, opt.tol = 1e-16)

# linear regression
info$Divergence <- node.depth.edgelength(rooted.tree)[1:Ntip(rooted.tree)]

model <- lm(Divergence ~ Date, data = info, subset = Query == 0)
null.model <- lm(Divergence ~ 1, data = info, subset = Query == 0)

est.date <- as.data.frame(do.call(
	rbind,
	lapply(info$Divergence, function(x)
		unlist(inverse.predict(model, x))
	)
))
root.date <- inverse.predict(model, 0)

# make output data frames
data <- data.frame(
	ID = info$ID,
	EstimatedDate = as.character(to.Date(est.date$Prediction), format = DATE_FMT),
	EstimatedDate95Low = as.character(to.Date(est.date$`Confidence Limits1`), format = DATE_FMT),
	EstimatedDate95High = as.character(to.Date(est.date$`Confidence Limits2`), format = DATE_FMT),
	stringsAsFactors = FALSE
)
data.censored <- data[info$Query == 1, ]

stats <- data.frame(
	RunID = run.id,
	dAIC = AIC(null.model) - AIC(model),
	EstimatedRootDate = as.character(to.Date(root.date$Prediction), format = DATE_FMT),
	EstimatedRootDate95Low = as.character(to.Date(root.date$`Confidence Limits`[1]), format = DATE_FMT),
	EstimatedRootDate95High = as.character(to.Date(root.date$`Confidence Limits`[2]), format = DATE_FMT),
	EstimatedEvolutionaryRate = coef(model)[[2]],
	Fit = as.numeric(((AIC(null.model) - AIC(model)) > 10) &&
					 	(root.date$`Confidence Limits`[1] < min(info$Date))),
	stringsAsFactors = FALSE
)

# write output (rooted tree, data and stats)
write.tree(rooted.tree, rooted.tree.file)
write.csv(data.censored, data.file, row.names = FALSE)
write.csv(stats, stats.file, row.names = FALSE)