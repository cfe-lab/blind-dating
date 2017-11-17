#!/usr/bin/Rscript

library(ggplot2)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

data.rtt <- lapply(dir("stats", "[0-9].cens.data", full.names=T), read.csv)
stats.rtt <- do.call(rbind, lapply(dir("stats", "[0-9].cens.stats", full.names=T), read.csv))

good.rtt <- with(stats.rtt, null.AIC - AIC > 10 & Estimated.Root.Date < 0)
data.rtt <- data.rtt[good.rtt]

data.rtt <- lapply(data.rtt, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); x})

g <- ggplot() + theme_bw() + theme(panel.grid=element_blank())

for (x in data.rtt) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', alpha=1/6, color="#00000000")
}

pdf("ancre.density.rtt.pdf")
g + scale_x_continuous(name="Scaled Date Differnce") + scale_y_continuous(name="Density")
dev.off()

data.ogr <- lapply(dir("stats", "cens-with_ref.data", full.names=T), read.csv)
stats.ogr <- do.call(rbind, lapply(dir("stats", "cens-with_ref.stats", full.names=T), read.csv))

good.ogr <- with(stats.ogr, null.AIC - AIC > 10 & Estimated.Root.Date < 0)
data.ogr <- data.ogr[good.ogr]

data.ogr <- lapply(data.ogr, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); x})

g <- ggplot() + theme_bw() + theme(panel.grid=element_blank())

for (x in data.ogr) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', alpha=1/5, color="#00000000")
}

pdf("ancre.density.ogr.pdf")
g + scale_x_continuous(name="Scaled Date Differnce") + scale_y_continuous(name="Density")
dev.off()