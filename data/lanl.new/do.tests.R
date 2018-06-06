#!/usr/bin/Rscript

library(chemCal)
library(Barnard)

get.deviant <- function(data, g, region) {
	ci <- do.call(rbind, lapply(data$Divergence, function(x) unlist(inverse.predict(g, x))))[, 4:5]
	region <- tolower(substr(region, 1, 1))
	
	subset(data, ((region == 'l' | region == 'o') & Collection.Date < ci[, 1]) | ((region == 'r' | region == 'o') & Collection.Date > ci[, 2]) | region == 'i')
}

do.all <- function(i) {
	fast.data <- get.deviant(all.data[[i]], all.g[[i]], 'LEFT')
	latent.data <- get.deviant(all.data[[i]], all.g[[i]], 'RIGHT')
	good.data <- get.deviant(all.data[[i]], all.g[[i]], 'INSIDE')
	
	cat(all.pat.ids[i], "\n")
	
	print(matrix(c(sum(latent.data$Type == 'PBMC'), sum(good.data$Type == 'PBMC'), sum(fast.data$Type == 'PBMC'), sum(latent.data$Type == 'PLASMA'), sum(good.data$Type == 'PLASMA'), sum(fast.data$Type == 'PLASMA')), ncol=2))
	
	cat("DNA v RNA contigency\n")
	barnard.test(sum(fast.data$Type == 'PBMC') + sum(good.data$Type == 'PBMC'), sum(fast.data$Type == 'PLASMA') + sum(good.data$Type == 'PLASMA'), sum(latent.data$Type == 'PBMC'), sum(latent.data$Type == 'PLASMA'))
	cat("Latent DNA v RNA\n")
	if (nrow(latent.data) > 0)
		print(binom.test(sum(latent.data$Type == 'PBMC'), nrow(latent.data)))
	else
		cat("N.A.\n")
	cat("Slow v Fast\n")
	if (nrow(latent.data) + nrow(fast.data) > 0)
		print(binom.test(c(nrow(latent.data), nrow(fast.data))))
	else
		cat("N.A.\n")
	cat("Slow v Fast (only DNA)\n")
	if (sum(latent.data$Type == 'PBMC') + sum(fast.data$Type == 'PBMC') > 0)
		print(binom.test(c(sum(latent.data$Type == 'PBMC'), sum(fast.data$Type == 'PBMC'))))
	else
		cat("N.A.\n")
}

skew.pat.ids <- paste0("patient_", c(821, 13889))

skew.data <- lapply(skew.pat.ids, function(x) read.csv(paste0("stats/", x, ".data.csv")))
skew.g <- lapply(skew.pat.ids, function(x) readRDS(paste0("stats/", x, ".regression.rds")))

skew.fast.data <- do.call(rbind, lapply(1:length(skew.pat.ids), function(i) get.deviant(skew.data[[i]], skew.g[[i]], 'LEFT')))
skew.latent.data <- do.call(rbind, lapply(1:length(skew.pat.ids), function(i) get.deviant(skew.data[[i]], skew.g[[i]], 'RIGHT')))
skew.good.data <- do.call(rbind, lapply(1:length(skew.pat.ids), function(i) get.deviant(skew.data[[i]], skew.g[[i]], 'INSIDE')))

cat("Only good\n")

print(matrix(c(sum(skew.latent.data$Type == 'PBMC'), sum(skew.good.data$Type == 'PBMC'), sum(skew.fast.data$Type == 'PBMC'), sum(skew.latent.data$Type == 'PLASMA'), sum(skew.good.data$Type == 'PLASMA'), sum(skew.fast.data$Type == 'PLASMA')), ncol=2))

cat("DNA v RNA contigency\n")
barnard.test(sum(skew.fast.data$Type == 'PBMC') + sum(skew.good.data$Type == 'PBMC'), sum(skew.fast.data$Type == 'PLASMA') + sum(skew.good.data$Type == 'PLASMA'), sum(skew.latent.data$Type == 'PBMC'), sum(skew.latent.data$Type == 'PLASMA'))
cat("Latent DNA v RNA\n")
binom.test(sum(skew.latent.data$Type == 'PBMC'), nrow(skew.latent.data))
cat("Slow v Fast\n")
binom.test(c(nrow(skew.latent.data), nrow(skew.fast.data)))
cat("Slow v Fast (only DNA)\n")
binom.test(c(sum(skew.latent.data$Type == 'PBMC'), sum(skew.fast.data$Type == 'PBMC')))

all.pat.ids <- paste0("patient_", c(821, 822, 824, 13889, 10769, 34391, 34411))

all.data <- lapply(all.pat.ids, function(x) read.csv(paste0("stats/", x, ".data.csv")))
all.g <- lapply(all.pat.ids, function(x) readRDS(paste0("stats/", x, ".regression.rds")))

all.fast.data <- do.call(rbind, lapply(1:length(all.pat.ids), function(i) get.deviant(all.data[[i]], all.g[[i]], 'LEFT')))
all.latent.data <- do.call(rbind, lapply(1:length(all.pat.ids), function(i) get.deviant(all.data[[i]], all.g[[i]], 'RIGHT')))
all.good.data <- do.call(rbind, lapply(1:length(all.pat.ids), function(i) get.deviant(all.data[[i]], all.g[[i]], 'INSIDE')))

sup <- lapply(1:length(all.pat.ids), do.all)

cat("All\n")

print(matrix(c(sum(all.latent.data$Type == 'PBMC'), sum(all.good.data$Type == 'PBMC'), sum(all.fast.data$Type == 'PBMC'), sum(all.latent.data$Type == 'PLASMA'), sum(all.good.data$Type == 'PLASMA'), sum(all.fast.data$Type == 'PLASMA')), ncol=2))

cat("DNA v RNA contigency\n")
barnard.test(sum(all.fast.data$Type == 'PBMC') + sum(all.good.data$Type == 'PBMC'), sum(all.fast.data$Type == 'PLASMA') + sum(all.good.data$Type == 'PLASMA'), sum(all.latent.data$Type == 'PBMC'), sum(all.latent.data$Type == 'PLASMA'))
cat("Latent DNA v RNA\n")
binom.test(sum(all.latent.data$Type == 'PBMC'), nrow(all.latent.data))
cat("Slow v Fast\n")
binom.test(c(nrow(all.latent.data), nrow(all.fast.data)))
cat("Slow v Fast (only DNA)\n")
binom.test(c(sum(all.latent.data$Type == 'PBMC'), sum(all.fast.data$Type == 'PBMC')))