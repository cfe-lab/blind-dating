## rtt.R (2014-06-16)

##   Root a tree by root-to-tip regression

## Copyright (c) 2014, Rosemary McCloskey, BC Centre for Excellence in HIV/AIDS

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rtt <- function(t, tip.dates, ncpu = 1, objective = "correlation", opt.tol = .Machine$double.eps^0.25) 
{
  
  ## These are objective functions which can be used to evaluate the "goodness" of
  ## a regression fit.
  if (objective == "correlation")
    objective <- function(x, y) cor.test(y, x)$estimate
  else if (objective == "rsquared")
    objective <- function(x, y) summary(lm(y ~ x))$r.squared
  else if (objective == "rms")
    objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
  else
    stop('objective must be one of "correlation", "rsquared", or "rms"')
  
  # Find the tips of the tree that're missing dates
  missing.indices <- which(is.na(tip.dates))
  valid.indices <- which(!is.na(tip.dates))
  
  tip.lengths <- node.depth.edgelength(t)
  
  # Computes distance-to-tips, indexing [i,] gives the distances from the i'th node to all tips
  if(is.rooted(t))  # NOTE: The last version unrooted the tree if it was rooted, this threw away information if the tree was rooted
                    # now we don't unroot it, but pull the distances off matrix as if it were unrooted.
    dist <- dist.nodes(t)[, 1:(t$Nnode + 1)] 
  else
    dist <- dist.nodes(t)[, 1:(t$Nnode + 2)]
  
  ## Do root-to-tip regressions for every possible choice of root,
  ## that is, compute the objective function for every choice of root.
  choice.f <- function(row) {
    if(row %in% missing.indices) # don't compute the objective function for tips that're missing data
      return (-Inf)
    else
      return (objective(tip.dates[valid.indices], dist[row, valid.indices])) # Only do the regression over data that exist
  }
  
  # Apply the objective function
  fits <- if (ncpu > 1)
    unlist(parallel::mclapply(1:nrow(dist), choice.f, mc.cores = ncpu))
  else 
    unlist(lapply(1:nrow(dist), choice.f))

  ## Find the best one (highest value of objective function).
  fit.edge <- apply(t$edge, 2, function(e) fits[e])
  obj.edge <- apply(fit.edge, 1, mean)
  
  ## Compatibility with Path-O-Gen: find the last maximum, not the first.
  best.edge <- length(obj.edge) - which.max(rev(obj.edge)) + 1
  best.edge.parent <- t$edge[best.edge, 1]
  best.edge.child <- t$edge[best.edge, 2]
  best.edge.length <- t$edge.length[best.edge]
  
  ## Find the best location on that edge.
  f <- function(x) {
    dist <- x * dist[best.edge.parent, ] + (1 - x) * dist[best.edge.child, ]
    objective(tip.dates, dist)
  }
  best.pos <- optimize(f, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum
  
  ## Reroot the tree at the optimal location
  new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", edge.length = 1, Nnode = 1L, root.edge = 1)
  class(new.root) <- "phylo"
  t <- bind.tree(t, new.root, where = best.edge.child, position = best.pos * best.edge.length)
  t <- collapse.singles(t)
  t <- root(t, "new.root")
  t <- drop.tip(t, "new.root")

  # Construct the linear model of the data
  tip.lengths <- node.depth.edgelength(t)
  
  distances <- tip.lengths[valid.indices]
  times <- tip.dates[valid.indices]
  model <- lm(times ~ distances)
  
  # Predict the missing dates
  tip.dates[missing.indices] <- predict(model, data.frame(distances=tip.lengths[missing.indices]))
  t$tip.dates <- tip.dates
  return (t)
}
