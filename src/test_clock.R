library(ape)

# read info file code
args.all <- commandArgs(trailingOnly = F)
source.dir <- "../../src"
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

# compute the diffrence of each pair of elememts of each vector
cross.diff <- function(x, y) {
	unlist(lapply(x, subtract, y))
}

# compute relative rate function
compute.relative.rate <- function(seq.A, seq.B, tree) {
	node.depth <- node.depth.edgelength(tree)
	
	T.A <- seq.A$NDATE
	T.B <- seq.B$NDATE
	d.A <- node.depth[as.numeric(rownames(seq.A))]
	d.B <- node.depth[as.numeric(rownames(seq.B))]
	
	cross.diff(d.A, d.B) / cross.diff(T.A, T.B)
}

# store distance
get.dists <- function(data, tree) {
	node.depth <- node.depth.edgelength(tree)
	
	data.frame(date=data$COLDATE, dist=node.depth[as.numeric(rownames(data))])
}

# read parameters
op <- OptionParser()
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
args <- parse_args(op)

tree.file <- args$tree
info.file <- args$info

# read inputs
tree <- ape::read.tree(tree.file)
info <- read.info(info.file, tree$tip.label)

#split info by time
info <- subset(info, CENSORED == 0)
info$NDATE <- as.numeric(as.Date(info$COLDATE))
info.s <- split(info, info$COLDATE)
n.t <- length(info.s)

# compute relative rates
relative.rates <- unlist(lapply(1:(n.t - 1), function(i) lapply((i+1):n.t, function(j)  compute.relative.rate(info.s[[i]], info.s[[j]], tree))))

dists <- get.dists(info, tree)








cross.self <- function(x, foo, ...) unlist(lapply(x, foo, x, ...))

data <- read.csv(data.file)
data <- subset(data, Censored == 0)

rr <- cross.self(data$Divergence, subtract) / cross.self(data$Collection.Date, subtract)
names(rr) <- cross.self(data$ID, paste, sep=".")
rr <- rr[is.finite(rr) & !is.na(rr)]

top <- mean(rr) + 2 * sd(rr)
bot <- mean(rr) - 2 * sd(rr)

bad <- rr[rr < bot | rr > top] %>% names %>% strsplit(split="\\.") %>% unlist %>% table %>% sort

sort(bad / sapply(names(bad), function(x) sum(data[x == data$ID, 'Collection.Date'] != data$Collection.Date)))
