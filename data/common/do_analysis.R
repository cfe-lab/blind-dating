#!/usr/bin/Rscript

library(nlme)
library(BSDA)

# get args
args <- commandArgs(trailing=T)

suffix <- args[1]
legend.posx <- as.double(args[2])
legend.posy <- as.double(args[3])
legend.ncols <- as.integer(args[4])

# read files
df <- read.csv(paste("data", args[1], ".csv", sep=""), header=T)

good <- read.csv("good_fits.txt", header=F)
names(good) <- c("id", "name")

# find good data
df.good <- df[df$Patient %in% good$id,]

df.good$Is.Latent <- df.good$Error > 0
df.good$Trans.Error <- sign(df.good$Error)*sqrt(abs(df.good$Error))

# graphing
n <- length(good$id)
colors <- unlist(lapply(seq(1, n), function(x) {
	if (x < n / 3) {
		return(rgb(1 - x * 3 / n, x * 3 / n, 0, .2))
	} else if (x < 2 * n / 3) {
		return(rgb(0, 2 - x * 3 / n, x * 3 / n - 1, .2))
	} else {
		return(rgb(x * 3 / n - 2, 0, 3 - x * 3 / n, .2))
	}
}))
names(colors) <- good$id

c <- NA

plot.graph <- function(y) {
	with(df.good, plot(Distance, Error, pch=16, col=colors[as.character(Patient)], ylab="Difference from Date Reconstruction (days)", xlab="Genetic Distance from Root (no. substitutions)"))
	x <- seq(0, 1, by=0.0001)
	lines(x, unlist(lapply(x, y)))
	abline(0, 0, col=rgb(0, 0, 0, .5))
	for (i in 1:2)
		legend(legend.posx, legend.posy, good$name, col=colors[as.character(good$id)], pch=16, ncol=legend.ncols)
}

do.analysis <- function(y, glm.function, ...) {
	model <- glm.function(..., data=df.good)

	plot.new()
	summary.output <- capture.output(print(summary(model)))
	for (i in 1:(length(summary.output))) {
		text(.5, 1 - i/length(summary.output), summary.output[i])
	}
	
	EDA(resid(model))
		
	assign("c", unlist(lapply(coefficients(model), mean)) * sd(df.good$Error), envir=.GlobalEnv)
	plot.graph(y)
}

do.analysis.trans <- function(y, glm.function, ...) {
	model <- glm.function(..., data=df.good)
	
	plot.new()
	summary.output <- capture.output(print(summary(model)))
	for (i in 1:(length(summary.output))) {
		text(.5, 1 - i/length(summary.output), summary.output[i])
	}
	
	EDA(resid(model))
	
	assign("c", unlist(lapply(coefficients(model), mean)) * sd(df.good$Trans.Error), envir=.GlobalEnv)
	plot.graph(y)
}

# binomial test
bin.test <- binom.test(sum(df.good$Is.Latent), length(df.good$Is.Latent), alternative="greater")
print(bin.test)

# scale factors
sd(df.good$Error)
sd(df.good$Trans.Error)

# GLM fitting
pdf(paste("lme.fixed.null", args[1], ".pdf", sep=""))
try(do.analysis(function(x){c[[1]]}, glm, formula=Error/sd(Error) ~ 1))
dev.off()

pdf(paste("lme.fixed.dist", args[1], ".pdf", sep=""))
try(do.analysis(function(x){c[[1]]+c[[2]]*x}, glm, formula=Error/sd(Error) ~ Distance))
dev.off()

pdf(paste("lme.mixed.null", args[1], ".pdf", sep=""))
try(do.analysis(function(x){c[[1]]}, lme, fixed=Error/sd(Error) ~ 1, random=~1 | Patient))
dev.off()

pdf(paste("lme.mixed.dist", args[1], ".pdf", sep=""))
try(do.analysis(function(x){c[[1]]+c[[2]]*x}, lme, fixed=Error/sd(Error) ~ Distance, random=~1 | Patient))
dev.off()

pdf(paste("lme.mixed.rand.dist", args[1], ".pdf", sep=""))
try(do.analysis(function(x){c[[1]]+c[[2]]*x}, lme, fixed=Error/sd(Error) ~ Distance, random=~Distance | Patient))
dev.off()

# GLM fitting (with transform)
pdf(paste("lme.fixed.null.sqrt", args[1], ".pdf", sep=""))
try(do.analysis.trans(function(x){sign(c[[1]])*(c[[1]])^2}, glm, formula=Trans.Error/sd(Trans.Error) ~ 1))
dev.off()

pdf(paste("lme.fixed.dist.sqrt", args[1], ".pdf", sep=""))
try(do.analysis.trans(function(x){sign(c[[1]]+c[[2]]*x)*(c[[1]]+c[[2]]*x)^2}, glm, formula=Trans.Error/sd(Trans.Error) ~ Distance))
dev.off()

pdf(paste("lme.mixed.null.sqrt", args[1], ".pdf", sep=""))
try(do.analysis.trans(function(x){sign(c[[1]])*(c[[1]])^2}, lme, fixed=Trans.Error/sd(Trans.Error) ~ 1, random=~1 | Patient))
dev.off()

pdf(paste("lme.mixed.dist.sqrt", args[1], ".pdf", sep=""))
try(do.analysis.trans(function(x){sign(c[[1]]+c[[2]]*x)*(c[[1]]+c[[2]]*x)^2}, lme, fixed=Trans.Error/sd(Trans.Error) ~ Distance, random=~1 | Patient))
dev.off()

pdf(paste("lme.mixed.rand.dist.sqrt", args[1], ".pdf", sep=""))
try(do.analysis.trans(function(x){sign(c[[1]]+c[[2]]*x)*(c[[1]]+c[[2]]*x)^2}, lme, fixed=Trans.Error/sd(Trans.Error) ~ Distance, random=~Distance | Patient))
dev.off()
