#!/usr/bin/Rscript

#library(nlme)
#library(BSDA)
#library(Rmpfr)
library(lme4)
#library(BSDA)

# bug fix
dnorm <- function (x, mean = 0, sd = 1, log = FALSE) {
    if (is.numeric(x) && is.numeric(mean) && is.numeric(sd)) 
        stats__dnorm(x, mean, sd, log = log)
    else if ((x.mp <- is(x, "mpfr")) || is(mean, "mpfr") || (s.mp <- is(sd, 
        "mpfr"))) {
        s.mp <- is(sd, "mpfr")
        prec <- pmax(53, getPrec(x), getPrec(mean), getPrec(sd))
        if (!x.mp) 
            x <- mpfr(x, prec)
        x <- (x - mean)/sd
        twopi <- 2 * Const("pi", prec)
        if (!s.mp) 
            sd <- mpfr(sd, prec)
        if (log) 
            -(log(sd) + (log(twopi) + x * x)/2)
        else exp(-x^2/2)/(sd * sqrt(twopi))
    }
    else stop("invalid arguments (x,mean,sd)")
}

pnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
    if (is.numeric(q) && is.numeric(mean) && is.numeric(sd)) 
        stats__pnorm(q, mean, sd, lower.tail = lower.tail, log.p = log.p)
    else if ((q.mp <- is(q, "mpfr")) || is(mean, "mpfr") || is(sd, 
        "mpfr")) {
        stopifnot(length(lower.tail) == 1, length(log.p) == 1)
        rr <- q <- ((if (q.mp) 
            q
        else as(q, "mpfr")) - mean)/sd
        if (any(neg <- (q < 0)))
            rr[neg] <- pnorm(-q[neg], lower.tail = !lower.tail, 
                log.p = log.p)
        if (any(pos <- !neg)) {
            q <- q[pos]
            prec.q <- max(.getPrec(q))
            rt2 <- sqrt(mpfr(2, prec.q))
            rr[pos] <- if (lower.tail) {
                eq2 <- erf(q/rt2)
                if (log.p && any(sml <- abs(eq2) < 0.5)) {
                  r <- q
                  r[sml] <- log1p(eq2[sml]) - log(2)
                  r[!sml] <- log((1 + eq2[!sml])/2)
                  r
                }
                else {
                  r <- (1 + eq2)/2
                  if (log.p) 
                    log(r)
                  else r
                }
            }
            else {
                r <- erfc(q/rt2)/2
                if (log.p) 
                  log(r)
                else r
            }
        }
        rr
    }
    else stop("(q,mean,sd) must be numeric or \"mpfr\"")
}

paste.0 <- function(...) {paste(..., sep='')}

# get args
args <- commandArgs(trailing=T)

suffix <- args[1]
data.name <- args[2]
#legend.posx <- as.double(args[2])
#legend.posy <- as.double(args[3])
#legend.ncols <- as.integer(args[4])

# read files
df <- read.csv(paste0("data", suffix, ".csv"), header=T)

good <- read.csv("good_fits.txt", header=F)
names(good) <- c("id", "name")

# find good data
df.good <- df[df$Patient %in% good$id,]

df.good$Is.Latent <- df.good$Error > 0
#df.good$Trans.Error <- sign(df.good$Error)*sqrt(abs(df.good$Error))

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

plot.hist <- function(err, colour, main) {
	d <- density(err)
	
	xlim <- c(min(err) * 1.1, max(err) * 1.1)
	ylim <- c(0, max(d$y) * 1.1)
	
	plot(c(-1001,-1100), xlim=xlim, ylim=ylim, xlab="Difference", ylab="Density", axes=T, main=main)
	
	polygon(d, col=colour, border=rgb(1,1,1,0))
	
	xlim
}

qqplotdist <- function(x, dist) {
	xu <- unique(sort(x))
	xq <- unlist(lapply(xu, function(i) sum(x == i))) / (length(x) + 1)
	for (i in 1:(length(xu)-1)) {
		xq[i + 1] <- xq[i + 1] + xq[i]
	}
	
	plot(xu, dist(xq))
}

# set up pdf printing
#pdf(paste("analysis", suffix, ".pdf", sep=""))

print.to.plot <- function(x) {
	plot.new()
	output <- capture.output(print(x))
	for (i in 1:(length(output))) {
		text(.5, 1 - i/length(output), output[i])
	}
}

print.lines.to.plot <- function(x) {
	plot.new()
	for (i in 1:(length(x))) {
		output <- capture.output(cat(x[i]))
		if (length(output) > 0)
			text(.5, 1 - i/length(x), output)
	}
}

do.test <- function(df) {
	b <- binom.test(sum(df$Is.Latent), length(df$Is.Latent))
	b$patient <- df$Patient[1]
	b
}

get.est <- function(x) 1/(1+exp(-x))

#print.to.plot(data.name)

# binomial test
bin.test <- do.test(df.good)
#print.to.plot(bin.test)

bin.glme.test <- glmer(Is.Latent ~ (1 | Patient), data=df.good, family=binomial)
bin.glm.test <- glm(Is.Latent ~ 1, data=df.good, family=binomial)
#s <- summary(bin.glme.test)
#print.to.plot(s)
#print.to.plot(c(paste0("AIC: ", AIC(bin.glme.test), ", null AIC: ", AIC(glm(Is.Latent ~ 1, data=df.good, family=binomial))), paste0("p-value: ", 1-pchisq(AIC(glm(Is.Latent ~ 1, data=df.good, family=binomial))- AIC(bin.glme.test) + 2, df=1)), paste0("mean: ", 1/(1 + exp(-coef(s)[, "Estimate"])))))
#EDA(resid(bin.glme.test))
#dev.off()

bin.tests <- lapply(split(df.good, df.good$Patient), do.test)

results <- t(rbind(sapply(bin.tests, function(x) data.frame(patient=c(x$patient), value=c(x$estimate), conf.low=c(x$conf.int[1]), conf.high=c(x$conf.int[2]), p.value=(x$p.value), stringsAsFactors=F))))

results <- rbind(data.frame(patient=c("all"), value=c(bin.test$estimate), conf.low=c(bin.test$conf.int[1]), conf.high=c(bin.test$conf.int[2]), p.value=c(bin.test$p.value), stringsAsFactors=F), results)

write.csv(data.frame(lapply(results, as.character), stringsAsFactors=FALSE), file=paste0("analysis.bin", suffix, ".csv"), row.names=F)

write.csv(data.frame(patient=c("all (glme)", "all (null)"), value=get.est(c(mean(unlist(coef(bin.glme.test))), unlist(coef(bin.glm.test)))), AIC=c(AIC(bin.glme.test), AIC(bin.glm.test)), stringsAsFactors=F), file=paste0("analysis.glme", suffix, ".csv"), row.names=F)

write.table(data.frame(patient=row.names(coef(bin.glme.test)[[1]]), value=get.est(unlist(coef(bin.glme.test)))), file=paste0("analysis.glme", suffix, ".csv"), row.names=F, append=T, sep=",", col.names=F)

if (0) {
# norm exp fit
logdnormexp <- function(x, lambda, sigma) {log(abs(lambda))+abs(lambda)/2*(abs(lambda)*sigma^2-2*sign(lambda)*x)+pnorm(-(abs(lambda)*sigma^2-sign(lambda)*x)/sigma, log=T)}

analysis.frame <- data.frame(Patient=c(good$name, "Overall"), lambda=NA, sigma=NA, lik=NA, norm.mu=NA, norm.sigma=NA, norm.lik=NA)

suppress <- apply(good, 1, function (x) {
		pat.error <- mpfr(unlist(df.good[df.good$Patient == x[1], "Error"]), 20)
		i <- which(x[1] == good[,1])
				
		lambda <- 1/mean(pat.error)
		analysis.frame[i,"lambda"] <<- as.double(lambda)
		sigma <- sqrt(sd(pat.error)^2-mean(pat.error)^2)
		analysis.frame[i,"sigma"] <<- as.double(sigma)
		lik <- if (is.na(sigma) || sigma <= 0) {
				-Inf
			} else {
				sum(logdnormexp(pat.error, lambda, sigma))
			}
		analysis.frame[i,"lik"] <<- as.double(lik)
				
		norm.mu <- mean(pat.error)
		analysis.frame[i,"norm.mu"] <<- as.double(norm.mu)
		norm.sigma <- sd(pat.error)
		analysis.frame[i,"norm.sigma"] <<- as.double(norm.sigma)
		norm.lik <- sum(dnorm(pat.error, norm.mu, norm.sigma, log=T))
		analysis.frame[i,"norm.lik"] <<- as.double(norm.lik)
		
		lims <- plot.hist(as.double(pat.error), "green", x[2])
		b <- seq(lims[1], lims[2], length.out=100)
		if (!is.infinite(lik))
			lines(b, exp(logdnormexp(b, lambda, sigma)), col="#FF000080")
		lines(b, dnorm(b, norm.mu, norm.sigma), col="#0000FF80")
				
#		print.lines.to.plot(c(paste.0("Patient ", x[2]), "", "Normal Fit: ", paste.0("mu: ", as.double(norm.mu)), paste.0("sigma: ", as.double(norm.sigma)), paste.0("likelihood: ", as.double(norm.lik)), "", "Normexp Fit: ", paste.0("lambda: ", as.double(lambda)), paste.0("sigma: ", as.double(sigma)), paste.0("likelihood: ", as.double(lik))))
	})

pat.error <- mpfr(unlist(df.good$Error), 20)
i <- length(good$id) + 1

lambda <- 1/mean(pat.error)
analysis.frame[i,"lambda"] <- as.double(lambda)
sigma <- sqrt(sd(pat.error)^2-mean(pat.error)^2)
analysis.frame[i,"sigma"] <- as.double(sigma)
lik <- if (is.na(sigma) || sigma == 0) {
		-Inf
	} else {
		sum(logdnormexp(pat.error, lambda, sigma))
	}
analysis.frame[i,"lik"] <- as.double(lik)
				
norm.mu <- mean(pat.error)
analysis.frame[i,"norm.mu"] <- as.double(norm.mu)
norm.sigma <- sd(pat.error)
analysis.frame[i,"norm.sigma"] <- as.double(norm.sigma)
norm.lik <- sum(dnorm(pat.error, norm.mu, norm.sigma, log=T))
analysis.frame[i,"norm.lik"] <- as.double(norm.lik)

lims <- plot.hist(as.double(pat.error), "green", "Overall")
b <- seq(lims[1], lims[2], length.out=100)
if (!is.infinite(lik))
	lines(b, exp(logdnormexp(b, lambda, sigma)), col="#FF000080")
lines(b, dnorm(b, norm.mu, norm.sigma), col="#0000FF80")

write.csv(analysis.frame, paste("analysis", suffix, ".csv", sep=""))
		
#print.lines.to.plot(c("Overall", "", "Normal Fit: ", paste.0("mu: ", as.double(norm.mu)), paste.0("sigma: ", as.double(norm.sigma)), paste.0("likelihood: ", as.double(norm.lik)), "", "Normexp Fit: ", paste.0("lambda: ", as.double(lambda)), paste.0("sigma: ", as.double(sigma)), paste.0("likelihood: ", as.double(lik))))


dev.off()
}
if (0) {
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
}