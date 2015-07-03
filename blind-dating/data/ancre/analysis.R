#!/usr/bin/Rscript

library(ape)
library(gplots)
library(lmtest)
library(MASS)
library(RColorBrewer)
library(BSDA)

source('include/rtt.R')

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)


options <- commandArgs(trailingOnly = TRUE)
types <- function(x) gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

if(length(options) != 3) {
	stop("Invalid usage!\nExample usage:\n    /analyze.tree.R [some newick tree] [root to tip 0=no, 1=root based on rtt] [output]\n");
	return;
}


# Parse command line args
tree.file <- options[1]
output.pdf <- options[3]
do.rtt <- F
if(options[2] == "1") do.rtt <- T

print(tree.file)

# Read in tree
tr <- read.tree(tree.file)
tr <- drop.tip(tr, "REFERENCE")

# Start pdf output
pdf(output.pdf, width=11.5, height=8.5)
plot(tr, cex=0.5)

# Extract the data from the tip labels
tip.dates <- extract_dates(tr$tip.label)
n <- length(tr$tip.label)

spercent <- floor(n*0.5)

tip.plasma <- rep(T, n)
tip.plasma[sample(1:n, spercent)] <- F
tip.pbmc <- !tip.plasma

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

ratio.test <- lrtest(model)
prb <- ratio.test["Pr(>Chisq)"][[1]][2]
a<-model$coefficients[[1]]
b<-model$coefficients[[2]]

pbmc.p.dates <-  pbmc.dists/b - a/b

# # Plot only the plasma data and fit
# k <- kde2d(plasma.dists, plasma.dates, n=300)
# image(k, col=r, main="Distance vs Time (No PBMC)",  xlab="Expected Substitutions", ylab="Time")
# points(plasma.dists, plasma.dates, main="Distance vs Time (No PBMC)",  xlab="Expected Substitutions", ylab="Time")
# # abline(model)
# abline(-a/b,1/b)

# Same but swap axis
k <- kde2d(plasma.dates, plasma.dists, n=100)
image(k, col=r, main="Distance vs Time (No PBMC)",  ylab="Expected Substitutions", xlab="Time")
points(plasma.dates, plasma.dists, main="Distance vs Time (No PBMC)",  ylab="Expected Substitutions", xlab="Time")
# abline(-a/b,1/b)
abline(model)

# # Plot both PBMC and Plasma
# k <- kde2d(plasma.dists, plasma.dates, n=300)
# plot(distances, tip.dates, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time")
# image(k, col=r, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time", add=T)
# points(distances, tip.dates, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time")
# # abline(model)
# abline(-a/b,1/b)
# points(pbmc.dists,pbmc.s.dates, col="red")

# Plot both PBMC and Plasma (swapped)
k <- kde2d(plasma.dates, plasma.dists, n=100)
plot(tip.dates, distances, main="Distance vs Time",  ylab="Expected Substitutions", xlab="Time")
image(k, col=r, main="Time vs Distance",  xlab="Expected Substitutions", ylab="Time", add=T)
points(tip.dates, distances, main="Distance vs Time",  ylab="Expected Substitutions", xlab="Time")
# abline(-a/b,1/b)
abline(model)
points(pbmc.s.dates, pbmc.dists, col="red")

# Plot a histogram of the difference between both the predicted and actual PBMC dates
difference <- pbmc.p.dates-pbmc.s.dates

m <- mean(difference)
hist(difference,  main="Histogram of Difference", sub="Time delta")
abline(v = m, col = "red", lwd = 2)
abline(v = 0, col = "blue", lwd = 2)

m <- mean(abs(difference))
hist(abs(difference), main="Histogram of Absolute Difference", sub="Time delta")
abline(v = m, col = "red", lwd = 2)

# Tree and test results
plot(tr, cex = 0.5)

textplot(sprintf("LR Test; slope = 0 for null:\n
	%f < %f ? (%s)\n
	%f < %f ? (%s)\n
	%f < %f ? (%s)\n",
	prb, 0.01, if(prb < 0.01) "Reject H_0; Good" else "Accept H_0; ",
	prb, 0.001, if(prb < 0.001) "Reject H_0; Good" else "Accept H_0; ",
	prb, 0.0001, if(prb < 0.0001) "Reject H_0; Good" else "Accept H_0; ")
)

EDA(residuals(model)) # where fit <- glm(...)

# Close the PDF
dev.off()

if(prb < 0.001) {
	file.rename(output.pdf, sprintf('good/%s', output.pdf))
} else {
	file.rename(output.pdf, sprintf('bad/%s', output.pdf))
}

