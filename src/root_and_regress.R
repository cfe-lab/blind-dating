library(ape)
library(optparse)
library(dplyr)
library(chemCal)

DATE_FMT <- "%Y-%m-%d"

# utility function
to.Date <- function(...) as.Date(..., origin = "1970-01-01")

get.child.lengths <- function(node) {
	sum(tree$edge.length[tree$edge[, 1] == node])
}

is.cherry <- function(node) {
	all(tree$edge[tree$edge[, 1] == node, 2] <= Ntip(tree))
}

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

if (!all(c("ID", "Date", "Query") %in% names(info))) {
	stop("Info file column names are incorrect")
}

info <- select(info, ID, Date, Query)
info <- info[match(tree$tip.label, info$ID), ]

if (any(is.na(info))) {
	stop("Info file missing data")
}

info$Date <- as.numeric(as.Date(info$Date, format = DATE_FMT))

if (any(is.na(info$Date))) {
	stop("Date format incorrect (should be yyyy-mm-dd)")
}

#browser()

# remove zero length branches if they exist
if (any(tree$edge.length < 1e-7)) {
	warning("Tiny length branches detected. This could be caused by duplicate sequences.")
}

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
	Query=info$Query,
	Date=info$Date,
	stringsAsFactors = FALSE
) %>%
	arrange(desc(Query), Date) %>%
	select(-Date)

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
write.csv(data, data.file, row.names = FALSE)
write.csv(stats, stats.file, row.names = FALSE)
