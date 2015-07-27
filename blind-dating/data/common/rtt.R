# Copyright (c) 2014, Rosemary McCloskey, BC Centre for Excellence in HIV/AIDS
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the BC Centre for Excellence in HIV/AIDS nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL The BC Centre for Excellence in HIV/AIDS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

library(ape)

if (R.Version()$major < 3 | R.Version()$minor < 1) {
  library(multicore)
} else {
  library(parallel)
}

# (Re)root a tree by root-to-tip regression.
#
# Rambaut, Andrew. "Estimating the rate of molecular evolution: incorporating 
# non-contemporaneous sequences into maximum likelihood phylogenies." 
# Bioinformatics 16.4 (2000): 395-399.
#
# t:         tree to root
#
# tip.dates: vector of dates for the tips, in the same order as t$tip.label
#
# ncpu:      number of cores to use
#
# objective: how to measure the "goodness" of the regression fit. Should be a
#            function taking two parameters (x, y), and returning a value which
#            is *larger* if the regression is better. Three objective functions,
#            which are the choices Path-O-Gen offers, are given below.
#
# opt.tol:   tolerance for optimization precision. By default, the optimize()
#            function uses a tolerance of .Machine$double.eps^0.25 (see ?optimize).
#            On my system, this is about 0.0001, which is probably not precise
#            enough for branch lengths, so you may want to set it to a smaller value.
rtt <- function (t, tip.dates, ncpu = 1, objective = "correlation", 
                 opt.tol = .Machine$double.eps^0.25)  {
  if (objective == "correlation") 
    objective <- function(x, y) cor.test(y, x)$estimate
  else if (objective == "rsquared") 
    objective <- function(x, y) summary(lm(y ~ x))$r.squared
  else if (objective == "rms") 
    objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
  else stop("objective must be one of \"correlation\", \"rsquared\", or \"rms\"")
  
  missing.indices <- which(is.na(tip.dates))
  valid.indices <- which(!is.na(tip.dates))
  
  ut <- unroot(t)
  dist <- dist.nodes(ut)[, 1:(ut$Nnode + 2)]
  
  # Save the tip labels, they might get reordered during the re-root
  # which would muss up th ordering of tip.dates
  saved.tips <- t$tip.label
  ut$tip.label <- unlist(lapply(1:length(ut$tip.label), toString)) # Give them 1..n, as strings
  
  f <- function (x, parent, child) {
    if(child %in% missing.indices) {
      return (-.Machine$double.max)
    }
    edge.dist <- x * dist[parent, valid.indices] + (1 - x) * dist[child, valid.indices]
    objective(tip.dates[valid.indices], edge.dist)
  }
  
  obj.edge <- if (ncpu > 1)
    unlist(parallel::mclapply(1:nrow(ut$edge), function (e) {
      opt.fun <- function (x) f(x, ut$edge[e,1], ut$edge[e,2])
      optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
    }, mc.cores=ncpu))
  else apply(ut$edge, 1, function (e) {
    opt.fun <- function (x) f(x, e[1], e[2])
    optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
  })
  
  best.edge <- which.max(obj.edge)
  
  best.edge.parent <- ut$edge[best.edge, 1]
  best.edge.child <- ut$edge[best.edge, 2]
  best.edge.length <- ut$edge.length[best.edge]
  
  opt.fun <- function (x) f(x, best.edge.parent, best.edge.child)
  best.pos <- optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum
  
  new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", 
                   edge.length = 1, Nnode = 1L, root.edge = 1)
  class(new.root) <- "phylo"
  ut <- bind.tree(ut, new.root, where = best.edge.child, position = best.pos * 
                    best.edge.length)
  ut <- collapse.singles(ut)
  ut <- root(ut, "new.root")
  rt <- drop.tip(ut, "new.root")
  
  # Reorder the tip dates, and restore the labels
  permutation <- as.integer(rt$tip.label)
  tip.dates <- tip.dates[permutation]
  rt$tip.label <- saved.tips[permutation]
  
  # Due to the potential reordering, we need to 
  # re-define these
  valid.indices <- which(!is.na(tip.dates))
  missing.indices <- which(is.na(tip.dates))
  
  # Construct the linear model of the data
  tip.lengths <- node.depth.edgelength(rt)
  
  distances <- tip.lengths[valid.indices]
  times <- tip.dates[valid.indices]
  model <- lm(distances ~ times)
  a<-model$coefficients[[1]]
  b<-model$coefficients[[2]]
  
  # Predict the missing dates
  tip.dates[missing.indices] <- (tip.lengths[missing.indices]/b - a/b)
  rt$tip.dates <- tip.dates
  rt
}
