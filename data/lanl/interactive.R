#!/usr/bin/Rscript

library(ape)
library(gplots)
#library(lmtest)
library(MASS)
library(RColorBrewer)
#library(BSDA)

setwd('/Users/art/git/cfe-papers/blind-dating/data/lanl/')
source('../common/rtt.R')

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)


#options <- commandArgs(trailingOnly = TRUE)
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))


# Parse command line args
tree.file <- 'trees.good/patient_10138.tre'
do.rtt <- FALSE


# Read in tree
tr <- read.tree(tree.file)
tr <- drop.tip(tr, "REFERENCE")

# Start pdf output


# Extract the data from the tip labels
tip.dates <- extract_dates(tr$tip.label)
tip.types <- types(tr$tip.label)

# Mark the PBMC and PLASMA cells
tip.pbmc <- tip.types == "PBMC"
tip.plasma <- tip.types == "PLASMA"

# Re-root the tree
if(do.rtt) {
	plasma.dates <- tip.dates
	plasma.dates[tip.pbmc] <- NA
	tr <- rtt(tr, plasma.dates)
}

# Create a regression line solely on the plasma data
plasma.dates <- tip.dates[tip.plasma]
pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates

distances <- node.depth.edgelength(tr)[1:length(tip.dates)]
plasma.dists <- distances[tip.plasma]
pbmc.dists <- distances[tip.pbmc]

# Build a linear model that relates distance(in) to dates(out)
# # model <- glm(plasma.dates ~ plasma.dists)
# model <- glm(plasma.dists ~ plasma.dates, family = "Gamma")
model <- glm(plasma.dists ~ plasma.dates)

# try fitting a model to first 3 time points of plasma only
model2 <- glm(plasma.dists[plasma.dates < 2000] ~ plasma.dates[plasma.dates < 2000])

#ratio.test <- lrtest(model)
#prb <- ratio.test["Pr(>Chisq)"][[1]][2]
a<-model2$coefficients[[1]]
b<-model2$coefficients[[2]]

pbmc.p.dates <-  pbmc.dists/b - a/b

# # Plot only the plasma data and fit
# k <- kde2d(plasma.dists, plasma.dates, n=300)
# image(k, col=r, main="Distance vs Time (No PBMC)",  xlab="Expected Substitutions", ylab="Time")
# points(plasma.dists, plasma.dates, main="Distance vs Time (No PBMC)",  xlab="Expected Substitutions", ylab="Time")
# # abline(model)
# abline(-a/b,1/b)

# Same but swap axis
# k <- kde2d(plasma.dates, plasma.dists, n=100)
# image(k, col=r, main="Distance vs Time (No PBMC)",  ylab="Expected Substitutions", xlab="Time")
# points(plasma.dates, plasma.dists, main="Distance vs Time (No PBMC)",  ylab="Expected Substitutions", xlab="Time")
# # abline(-a/b,1/b)
# abline(model)

# # Plot both PBMC and Plasma
# k <- kde2d(plasma.dists, plasma.dates, n=300)
# plot(distances, tip.dates, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time")
# image(k, col=r, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time", add=T)
# points(distances, tip.dates, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time")
# # abline(model)
# abline(-a/b,1/b)
# points(pbmc.dists,pbmc.s.dates, col="red")

# Plot both PBMC and Plasma (swapped)
#k <- kde2d(plasma.dates, plasma.dists, n=100)
par(mar=c(5,5,2,2))
plot(tip.dates, distances, ylab="Expected Substitutions", xlab="Time (days)", type='n', cex.lab=1.2, xlim=c(0, 2500))
#image(k, col=r, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time", add=T)
points(plasma.dates, plasma.dists)
# abline(-a/b,1/b)
#abline(model, lty=2)
points(pbmc.s.dates, pbmc.dists, col="red", pch=2, cex=0.5)
for (i in 1:length(pbmc.s.dates)) {
	lines(x=c(pbmc.s.dates[i], pbmc.p.dates[i]), y=rep(pbmc.dists[i], 2), col=rgb(1,0,0,0.2), lwd=1)
}
abline(model2, lty=2)
legend(x=1800, y=0.03, legend=c('Plasma', 'PBMC'), pch=c(1,2), col=c('black', 'red'))



par(mar=c(5,1,1,1))
hist(pbmc.p.dates-pbmc.s.dates, breaks=15, border=NA, col=rgb(0,0,1,0.3), main=NA, xlab='PBMC date residuals (days)', cex.axis=1.1, yaxt='n', ylab='', xlim=c(-1000, 1000), cex.lab=1.2)
abline(v=median(pbmc.p.dates-pbmc.s.dates), col='white', lwd=3, lend=0)

# Plot a histogram of the difference between both the predicted and actual PBMC dates
difference <- pbmc.p.dates-pbmc.s.dates

# m <- mean(difference)
# hist(difference,  main="Histogram of Difference", sub="Time delta")
# abline(v = m, col = "red", lwd = 2)
# abline(v = 0, col = "blue", lwd = 2)

# m <- mean(abs(difference))
# hist(abs(difference), main="Histogram of Absolute Difference", sub="Time delta")
# abline(v = m, col = "red", lwd = 2)

# Tree and test results
#plot(tr, cex = 0.5)

# textplot(sprintf("LR Test; slope = 0 for null:\n
	# %f < %f ? (%s)\n
	# %f < %f ? (%s)\n
	# %f < %f ? (%s)\n",
	# prb, 0.01, if(prb < 0.01) "Reject H_0; Good" else "Accept H_0; ",
	# prb, 0.001, if(prb < 0.001) "Reject H_0; Good" else "Accept H_0; ",
	# prb, 0.0001, if(prb < 0.0001) "Reject H_0; Good" else "Accept H_0; ")
# )

#EDA(residuals(model)) # where fit <- glm(...)

# Close the PDF
dev.off()

# if(prb < 0.001) {
	# file.rename(output.pdf, sprintf('good/%s', output.pdf))
# } else {
	# file.rename(output.pdf, sprintf('bad/%s', output.pdf))
# }

