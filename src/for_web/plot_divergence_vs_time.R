library(optparse)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidytree)
library(lubridate)

# Placeholders
STEPS <- 100
TREE_COLOUR <- "darkgray"
TREE_ALPHA <- 0.5
TREE_SIZE <- 1
TIP_SIZE <- 5
REGRESSION_SIZE <- 2
REGRESSION_COLOUR <- "blue"
REGRESSION_LINETYPE <- "dashed"
DIV_NAME <- "Divergence from Root"
DATE_NAME <- "Collection Date"
DATE_FMT <- "%Y-%m-%d"
EXPAND <- 0.025

COLOUR_BREAKS <- c("0", "1")
COLOUR_VALUES <- c("black", "red")
SHAPE_BREAKS <- c("0", "1")
SHAPE_VALUES <- c(16, 18)

# utility function
to.Date <- function(...) as.Date(..., origin = "1970-01-01")

# plot theme
my.theme <- theme_bw() +
	theme(
		panel.grid = element_blank(),
		legend.position = 'none',
		text = element_text(size=20),
		axis.ticks = element_line(colour = "black"),
		axis.text = element_text(colour = "black"),
	)

# from node.dating (https://github.com/brj1/node.dating/blob/random/R/node.dating.R)
# to process children before parents
get.node.order <- function(t) {
	n.tips <- length(t$tip.label)
	
	nodes <- n.tips + 1
	
	for (i in 1:(t$Nnode + n.tips)) {
		to.add <- t$edge[t$edge[, 1] == nodes[i], 2]
		
		nodes <- c(nodes, to.add[to.add > 0])
		
		i <- i + 1
	}
	
	nodes <- rev(nodes)
	
	nodes
}

# from node.dating (https://github.com/brj1/node.dating/blob/random/R/node.dating.R)
estimate.dates <- function(
	t,
	node.dates,
	mu = estimate.mu(t, node.dates, output.type='numeric'),
	node.mask = 1:length(tree$tip.label),
	node.order = get.node.order(t),
	min.date = -.Machine$double.xmax,
	max.date = .Machine$double.xmax,
	show.steps = 0,
	opt.tol = 1e-8,
	nsteps = 1000,
	lik.tol = 0,
	is.binary = is.binary.phylo(t),
	output.type = 'vector')
{
	# check parameters
	if (any(mu < 0))
		stop(paste0("mu (", mu, ") less than 0"))
	
	if (!(output.type %in% c('vector', 'list', 'phylo4d')))
		stop(paste0("unknown output type: ", output.type))
	if (output.type == 'phylo4d' && !require(phylobase))
		stop(paste0("library phylobase required for phylo4d output"))
	
	# init vars
	mu <- if (length(mu) == 1) {
		rep(mu, length(t$edge.length))
	} else if (length(mu) == nrow(t$edge)) {
		mu
	} else {
		stop(paste0("mu must be a vector with length equal 1 or equal to the number of edges"))
	}
	n.tips <- length(t$tip.label)
	dates <- if (length(node.dates) == n.tips) {
		c(node.dates, rep(NA, t$Nnode))
	} else if (length(node.dates) == n.tips + t$Nnode) {
		node.dates
	} else {
		stop(paste0(
			"node.dates must be a vector with length equal to the number of tips or ",
			"equal to the number of nodes plus the number of tips"
		))
	}
	
	lik.sens <- if (lik.tol == 0) opt.tol else lik.tol
	
	# Don't count initial step if all values are seeded
	iter.step <- if (any(is.na(dates))) 0 else 1
	
	children <- lapply(
		1:(t$Nnode + n.tips),
		function(x) {
			which(t$edge[,1] == x)
		}
	)
	parent <- lapply(
		1:(t$Nnode + n.tips),
		function(x) {
			which(t$edge[,2] == x)
		}
	)
	
	nodes <- node.order[!(node.order %in% node.mask)]
	
	min.dates <- dates
	min.dates[-node.mask] <- min.date
	
	for (n in rev(nodes)) {
		par <- t$edge[parent[[n]], 1]
		
		min.dates[n] <- max(min.dates[par], min.date)
	}
	
	max.dates <- dates
	max.dates[-node.mask] <- max.date
	
	for (n in rev(nodes)) {
		ch <- t$edge[children[[n]], 2]
		
		max.dates[n] <- min(max.dates[ch], max.date)
	}
	
	# calculate likelihood functions
	scale.lik <- sum(-lgamma(t$edge.length+1)+(t$edge.length+1)*log(mu))
	
	calc.Like <- function(ch.node, ch.edge, x) {
		tim <- ch.node - x
		
		t$edge.length[ch.edge]*log(tim)-mu[ch.edge]*tim
	}
	
	opt.fun <- function(x, ch, p, ch.edge, p.edge, use.parent=T) {
		sum(
			if (!use.parent || length(dates[p]) == 0 || is.na(dates[p])) {		
				calc.Like(dates[ch], t$edge.length[ch.edge], x)
			} else {
				calc.Like(
					c(dates[ch], x),
					c(t$edge.length[ch.edge],
					  t$edge.length[p.edge]),
					c(rep(x, length(dates[ch])), dates[p])
				)
			},
			na.rm=T
		)
	}
	
	get.bounds <- function(bounds) {
		x <- c(bounds[1] + opt.tol, bounds[2] - opt.tol)
		if (x[2] <= x[1])
			x <- mean(bounds)
		x
	}
	
	solve.lin2 <- function(bounds, p.times, p.edge) {	
		y <- p.times + t$edge.length[p.edge] / mu[p.edge]
		x <- get.bounds(bounds)
		if (bounds[1] < y && y < bounds[2])
			x <- c(y, x)
		
		x[which.max(unlist(lapply(x, function(z) sum(calc.Like(z, p.edge, p.times)))))]
	}
	
	solve.lin <- function(bounds, ch.times, ch.edge) {	
		y <- ch.times - t$edge.length[ch.edge] / mu[ch.edge]
		x <- get.bounds(bounds)
		if (bounds[1] < y && y < bounds[2])
			x <- c(y, x)
		
		x[which.max(unlist(lapply(x, function(z) sum(calc.Like(ch.times, ch.edge, z)))))]
	}
	
	solve.poly2 <- function(bounds, a, b, c.0) {
		x <- get.bounds(bounds)
		
		if (b ^ 2 - 4 * a * c.0 >= 0) {
			if (a == 0) {
				y <- -c.0 / b
				
				if (bounds[1] < y && y < bounds[2])
					x <- c(y, x)
			} else {
				x.1 <- (-b + sqrt(b ^ 2 - 4 * a * c.0)) / (2 * a)
				x.2 <- (-b - sqrt(b ^ 2 - 4 * a * c.0)) / (2 * a)
				
				if (bounds[1] < x.1 && x.1 < bounds[2])
					x <- c(x.1, x)
				if (bounds[1] < x.2 && x.2 < bounds[2])
					x <- c(x.2, x)
			}
		}
		
		x
	}
	
	solve.bin <- function(bounds, ch.times, ch.edge) {
		ch.edge.length <- t$edge.length[ch.edge]
		a <- sum(mu[ch.edge])
		b <- ch.edge.length[1] + ch.edge.length[2] - 
			a * (ch.times[1] + ch.times[2])
		c.0 <- a*ch.times[1] * ch.times[2] - 
			ch.times[1] * ch.edge.length[2] - 
			ch.times[2] * ch.edge.length[1]
		
		x <- solve.poly2(bounds, a, b, c.0)					
		
		x[which.max(unlist(lapply(x, function(y) sum(calc.Like(ch.times, ch.edge, y)))))]
	}
	
	
	solve.bin2 <- function(bounds, ch.times, ch.edge, par.time, par.edge) {
		ch.edge.length <- t$edge.length[ch.edge]
		par.edge.length <- t$edge.length[par.edge]
		a <- mu[ch.edge] - mu[par.edge]
		b <- ch.edge.length + par.edge.length - 
			a * (ch.times + par.time)
		c.0 <- a*ch.times * par.time - 
			ch.times * par.edge.length - 
			par.time * ch.edge.length
		
		x <- solve.poly2(bounds, a, b, c.0)					
		
		x[which.max(unlist(lapply(x, function(y)
			sum(calc.Like(c(ch.times, y), c(ch.edge, par.edge), c(y, par.time)))
		)))]
	}
	
	solve.poly3 <- function(bounds, a, b, c.0, d) {
		x <- get.bounds(bounds)
		
		if (a == 0)
			x <- c(x, solve.poly2(bounds, b, c.0, d))
		else {
			delta.0 <- complex(real=b^2 - 3 * a * c.0)
			delta.1 <- complex(real=2 * b^3 - 9 * a * b * c.0 + 27 * a^2 * d)
			C <- ((delta.1 + sqrt(delta.1^2 - 4 * delta.0^3)) / 2)^(1/3)
			
			x.1 <- Re(-1 / (3 * a) *
					  	(b + complex(real=1) * C + delta.0 / (complex(real=1) * C)))
			x.2 <- Re(-1 / (3 * a) *
					  	(b + complex(real=-1/2, imaginary=sqrt(3)/2) * C +
					  	 	delta.0 / (complex(real=-1/2, imaginary=sqrt(3)/2) * C)))
			x.3 <- Re(-1 / (3 * a) *
					  	(b + complex(real=-1/2, imaginary=-sqrt(3)/2) * C +
					  	 	delta.0 / (complex(real=-1/2, imaginary=-sqrt(3)/2) * C)))
			
			if (bounds[1] < x.1 && x.1 < bounds[2])
				x <- c(x.1, x)
			if (bounds[1] < x.2 && x.2 < bounds[2])
				x <- c(x.2, x)
			if (bounds[1] < x.3 && x.3 < bounds[2])
				x <- c(x.3, x)
		}
		
		x
	}
	
	solve.cube <- function(bounds, ch.times, ch.edge, par.time, par.edge) {
		ch.edge.length <- t$edge.length[ch.edge]
		par.edge.length <- t$edge.length[par.edge]
		
		a <- sum(mu[ch.edge]) - mu[par.edge]
		b <- sum(ch.edge.length) + par.edge.length - a * (sum(ch.times) + par.time)
		c.0 <- a * (ch.times[1] * ch.times[2] + ch.times[1] * par.time + ch.times[2] * par.time) -
			(ch.times[1] + ch.times[2]) * par.edge.length -
			(ch.times[1] + par.time) * ch.edge.length[2] -
			(ch.times[2] + par.time) * ch.edge.length[1]
		d <- ch.edge.length[1] * ch.times[2] * par.time +
			ch.edge.length[2] * ch.times[1] * par.time + 
			par.edge.length * ch.times[1] * ch.times[2] -
			a * prod(ch.times) * par.time
		
		x <- solve.poly3(bounds, a, b, c.0, d)
		
		x[which.max(unlist(lapply(x, function(y)
			sum(calc.Like(c(ch.times, y), c(ch.edge, par.edge), c(y, y, par.time)))
		)))]
	}
	
	estimate <- function(node) {
		ch.edge <- children[[node]]
		ch <- t$edge[ch.edge, 2]
		
		p.edge <- parent[[node]]
		p <- t$edge[p.edge, 1]
		
		m <- if (length(p) == 0 || is.na(dates[p])) {
			min.dates[node]
		} else {
			dates[p]
		}
		
		M <- min(max.dates[node], dates[ch], na.rm=T)
		
		good.dates <- !is.na(dates[ch])
		n.ch <- sum(good.dates)
		
		if (is.binary) {
			if (m + 2 * opt.tol >= M) {
				mean(c(m, M))
			} else {
				if (length(dates[p]) == 0 || is.na(dates[p])) {
					if (n.ch == 2)
						solve.bin(c(m, M), dates[ch], ch.edge)
					else
						solve.lin(c(m, M), dates[ch][good.dates], ch.edge[good.dates])
				} else {
					if (n.ch == 2)
						solve.cube(c(m, M), dates[ch], ch.edge, dates[p], p.edge)
					else if (n.ch == 1)
						solve.bin2(c(m, M), dates[ch][good.dates], ch.edge[good.dates], dates[p], p.edge)
					else
						solve.lin2(c(m, M), dates[p], p.edge)
				}
			}
		} else {				
			res <- optimize(opt.fun, c(m, M), ch, p, ch.edge, p.edge, maximum=T)
			
			res$maximum
		}
	}
	
	# iterate to estimate dates
	lik <- NA
	
	repeat
	{
		for (n in nodes) {		
			dates[n] <- estimate(n)
		}
		
		all.lik <- calc.Like(dates[t$edge[,2]], 1:length(t$edge.length), dates[t$edge[,1]]) + scale.lik
		new.lik <- sum(all.lik)
		
		if (show.steps > 0 && ((iter.step %% show.steps) == 0)) {
			cat(paste("Step: ", iter.step, ", Likelihood: ", new.lik, "\n", sep=""))
		}
		
		if (
			(lik.tol > 0 &&
			 (!is.na(lik) && (is.infinite(lik) || is.infinite(new.lik) || new.lik - lik < lik.tol))) ||
			(nsteps > 0 && iter.step >= nsteps) ||
			(lik.tol <= 0 && nsteps <= 0)
		) {
			if (is.infinite(lik) || is.infinite(new.lik)) {
				warning("Likelihood infinite")
			}
			else if (!is.na(lik) && new.lik + lik.sens < lik) {			
				warning("Likelihood less than previous estimate")
			}
			
			break
		} else {
			lik <- new.lik
		}
		
		iter.step <- iter.step + 1
	}
	
	if (show.steps > 0) {
		cat(paste("Step: ", iter.step, ", Likelihood: ", new.lik, "\n", sep=""))
	}
	
	if (output.type == 'vector')
		dates
	else if (output.type == 'phylo') {
		time.t <- t
		time.t$edge.length <- dates[t$edge[, 2]] - dates[t$edge[, 1]]
		
		time.t
	} else if (output.type == 'list') {
		time.t <- t
		time.t$edge.length <- dates[t$edge[, 2]] - dates[t$edge[, 1]]
		
		list(tree=t, time.tree=time.t, node.date=dates, mu=mu, log.lik=new.lik, edge.lik=all.lik)
	} else if (output.type == 'phylo4d') {
		from.edge <- unlist(lapply(1:(n.tips + t$Nnode), function(x) {
			if (any(t$edge[,2] == x)) which(t$edge[,2] == x) else NA
		}))
		parent <- t$edge[from.edge, 1]
		
		df <- data.frame(date=dates,
						 ancestor.date=dates[parent],
						 edge.time=dates-dates[parent],
						 edge.lik=all.lik[from.edge]
		)
		
		if (output.type == 'phylo4d')
			phylo4d(t, all.data=df, metadata=list(mu=mu, log.lik=new.lik))
	}
}

# parse command line arguments
op <- OptionParser()
op <- add_option(op, "--info", type = 'character')
op <- add_option(op, "--rootedtree", type = 'character')
op <- add_option(op, "--stats", type = 'character')
op <- add_option(op, "--plotprefix", type = 'character')
args <- parse_args(op)

info.file <- args$info
rooted.tree.file <- args$rootedtree
stats.file <- args$stats
plot.prefix <- args$plotprefix

### TODO check arguments

# read tree, info and stats
tree <- read.tree(rooted.tree.file)
info <- read.csv(info.file, stringsAsFactors = FALSE)
stats <- read.csv(stats.file, stringsAsFactors = FALSE)

info <- info[match(tree$tip.label, info$ID), ]
info$Date <- as.numeric(as.Date(info$Date, format = DATE_FMT))

# estimate node dates
root.date <- as.numeric(as.Date(stats$EstimatedRootDate, format = DATE_FMT))
mu <- stats$EstimatedEvolutionaryRate

node.dates <- tryCatch(
	{
		estimate.dates(
			tree,
			c(info$Date, root.date, rep(NA, Nnode(tree) - 1)),
			mu,
			node.mask = 1:(Ntip(tree) + 1),
			lik.tol=0,
			nsteps=STEPS,
			show.steps=0,
			opt.tol=1e-16
		)
	},
	error=function(e) {
		cat("Warning: Error in estimate.dates: ")
		message(e)
		cat("\n")
		estimate.dates(
			tree,
			info$Date,
			mu,
			lik.tol=0,
			nsteps=STEPS,
			show.steps=0,
			opt.tol=1e-16
		)
	}
)

info.combined <- mutate(info, node = 1:Ntip(tree)) %>% 
	full_join(
		data.frame(node = 1:(Ntip(tree) + Nnode(tree)), date = node.dates),
		by = "node"
	) %>%
	transmute(
		node = node,
		date = date,
		Censored = as.factor(Query)
	 )

# calculate date breaks
dates <- to.Date(info.combined$date)
date.range <- year(max(dates)) - year(min(dates))
sep <- max(1, round(date.range / 5))
date.labels <- dates %>%
	min() %>%
	as.character("%Y") %>%
	seq(by = sep, length.out = 7)
date.breaks <- date.labels %>%
	as.character() %>%
	as.Date(format = "%Y") %>%
	as.numeric()

# make plot
fort.tree <- as_tibble(tree) %>%
	full_join(info.combined, by = "node") %>%
	as.treedata() %>%
	fortify(layout = 'slanted', yscale = "date")

p <- ggplot(fort.tree, aes(x = x, y = y)) +
	geom_tree(layout = 'slanted', colour = TREE_COLOUR, size = TREE_SIZE, alpha = TREE_ALPHA) +
	geom_tippoint(aes(colour = Censored, shape = Censored), size = TIP_SIZE) +
	scale_x_continuous(name = DIV_NAME, expand = expansion(mult = EXPAND)) +
	scale_y_continuous(
		name = DATE_NAME,
		breaks = date.breaks,
		labels = date.labels,
		expand = expansion(mult = EXPAND)
	) +
	scale_colour_manual(breaks = COLOUR_BREAKS, limits = COLOUR_BREAKS, values = COLOUR_VALUES) +
	scale_shape_manual(breaks = SHAPE_BREAKS, limits = SHAPE_BREAKS, values = SHAPE_VALUES) +
	coord_flip() +
	my.theme

p$layers <- c(
	geom_abline(
		intercept = root.date,
		slope = 1 / mu,
		colour = REGRESSION_COLOUR,
		size = REGRESSION_SIZE,
		linetype = REGRESSION_LINETYPE
	),
	p$layer
)

# save plots
ggsave(
	paste0(plot.prefix, ".pdf"),
	p,
	'pdf',
	colormode = 'srgb',
	useDingbats = FALSE
)

ggsave(
	paste0(plot.prefix, ".png"),
	p,
	'png'
)