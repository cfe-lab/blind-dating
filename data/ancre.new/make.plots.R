#!/usr/bin/Rscript

library(ggplot2)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

my.theme <- theme(
	text=element_text(size=20),
	axis.text=element_text(size=15, colour='black'),
	legend.text=element_text(size=15),
	legend.position=0,
	axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
	legend.justification=c(0, 1),
	legend.spacing=unit(0, 'cm'),
	legend.margin=margin(0, 0, 0, 0, 'cm'),
	legend.key.size=unit(1.2, 'cm'),
	legend.key=element_rect(fill="#00000000", colour="#00000000"),
	legend.background=element_blank(),
	legend.box.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank()
) 

my.theme2 <- theme(
 text=element_text(size=35),
 axis.text=element_text(size=30, colour='black'),
 legend.text=element_text(size=25),
 legend.position=c(.99, .99),
 legend.justification=c(1, 1),
 legend.spacing=unit(0, 'cm'),
 legend.box='horizontal',
 legend.margin=margin(0.01, 0.01, 0.01, 0.01, 'cm'),
 legend.key.size=unit(1.2, 'cm'),
 legend.key=element_rect(fill="#00000000", colour="#00000000"),
 legend.background=element_blank(),
 legend.box.background=element_blank(),
 panel.grid.major=element_blank(),
 panel.grid.minor=element_blank()
)

concord <- function(x, y) {
	mu.x <- sum(x) / length(x)
	mu.y <- sum(y) / length(y)
	s.x <- sum((x - mu.x)^2) / length(x) 
	s.y <-  sum((y - mu.y)^2) / length(y)
	s.xy <- sum((x - mu.x) * (y - mu.y)) / length(y)
	
	2 * s.xy / (s.x + s.y + (mu.x - mu.y)^2)
}

my.theme.hist <- theme(
text=element_text(size=35),
axis.text=element_text(size=30, colour='black'),
axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=20),
legend.text=element_text(size=30),
legend.title=element_text(size=30),
legend.position=c(.98, .98),
legend.justification=c(1, 1),
legend.spacing=unit(0, 'cm'),
legend.margin=margin(0, 0, 0, 0, 'cm'),
legend.key.size=unit(1.2, 'cm'),
legend.background=element_blank(),
legend.box.background=element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank()
)

THERAPY_COLOUR <- "#a0a0a0"
THERAPY <- 1795
data <- read.csv("stats/patient_13334.cens.data.csv")
data <- subset(data, Censored == 1)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')
pdf("plots/patient_13334.cens.hist.pdf")
ggplot(data, aes(x=Estimated.Date)) +
#	geom_rect(xmin=THERAPY, xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR) +
	geom_histogram(breaks=seq(min(floor(data$Estimated.Date / 365.25)) * 365.25, max(floor(data$Estimated.Date / 365.25) + 1) * 365.25, by=365.25), fill='red') +
	geom_segment(x=0, xend=0, y=8, yend=7, arrow=arrow(length = unit(0.5, "cm")))  +
	theme_bw() +
	my.theme.hist +
	scale_x_continuous(name="Years since first collection", breaks=seq(min(floor(data$Estimated.Date / 365.25)) * 365.25, max(floor(data$Estimated.Date / 365.25) + 1) * 365.25, by=365.25), labels=seq(-1, 3)) +
	scale_y_continuous(name="Frequency")
dev.off()

files <- dir("info", "cens.csv")

data.files.rtt <- paste0("stats/", gsub(".csv", ".data.csv", files))
stats.files.rtt <- paste0("stats/", gsub(".csv", ".stats.csv", files))
data.rtt <- lapply(data.files.rtt, read.csv)
stats.rtt <- do.call(rbind, lapply(stats.files.rtt, read.csv))

good.rtt <- with(stats.rtt, Model.Fit == 1)
data.rtt <- data.rtt[good.rtt]

data.rtt <- lapply(data.rtt, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); subset(x, x$Censored == 1)})

g <- ggplot() + theme_bw() + my.theme2

for (x in data.rtt) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', alpha=1/length(data.rtt), color="#00000000")
}

pdf("ancre.density.rtt.pdf")
g + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 3))
dev.off()

data.files.ogr <- paste0("stats/", gsub(".csv", "-with_ref.data.csv", files))
stats.files.ogr <- paste0("stats/", gsub(".csv", "-with_ref.stats.csv", files))
data.ogr <- lapply(data.files.ogr, read.csv)
stats.ogr <- do.call(rbind, lapply(stats.files.ogr, read.csv))

good.ogr <- with(stats.ogr, Model.Fit == 1)
data.ogr <- data.ogr[good.ogr]

data.ogr <- lapply(data.ogr, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); subset(x, x$Censored == 1)})

g <- ggplot() + theme_bw() + my.theme2
for (x in data.ogr) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', alpha=1/length(data.ogr), color="#00000000")
}

pdf("ancre.density.ogr.pdf")
g + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 3))
dev.off()

data.rtt <- lapply(data.files.rtt, read.csv)
data.rtt <- data.rtt[good.ogr & good.rtt]
data.ogr <- lapply(data.files.ogr, read.csv)
data.ogr <- data.ogr[good.ogr & good.rtt]

data.all <- do.call(rbind, lapply(1:length(data.rtt), function(i) cbind(Patient=stats.rtt[good.ogr & good.rtt, "Patient"][i], merge(data.rtt[[i]], data.ogr[[i]], by=c("ID", "Type", "Censored", "Collection.Date"), suffixes=c(".rtt", ".ogr")))))
data.all <- subset(data.all, Censored == 1)

pdf("ancre.comp.pdf")
ggplot(data.all) + geom_abline() + geom_point(aes(x=Estimated.Date.rtt, y=Estimated.Date.ogr, colour=Patient), show.legend=F) + scale_colour_brewer(name="", palette='Dark2') + scale_x_continuous(name="Estimated Date RTT (days after first sampling)") + scale_y_continuous(name="Estimated Date OGR (days after first sampling)") + theme_bw() + my.theme
dev.off()

cat("Concordance:\n")
sup <- lapply(unique(data.all$Patient), function(x) cat(x, "; ", with(subset(data.all, Patient == x), concord(Estimated.Date.rtt, Estimated.Date.ogr)), "\n", sep=""))
cat("Total: ", with(data.all, concord(Estimated.Date.rtt, Estimated.Date.ogr)), "\n", sep="")