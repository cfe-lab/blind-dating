#!/usr/bin/Rscript
library(ggplot2)

my.theme <- theme(
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
data <- read.csv("stats/cartoon.data.csv")
data <- subset(data, Censored == 1)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')
pdf("plots/cartoon.hist.pdf")
ggplot(data, aes(x=Estimated.Date)) +
	geom_rect(xmin=2.1, xmax=Inf, ymin=0, ymax=2, fill=THERAPY_COLOUR) +
	geom_histogram(breaks=c(-1, 0, 1, 2, 3), fill='red') +
	geom_segment(x=0, xend=0, y=2, yend=1.5, arrow=arrow(length = unit(0.5, "cm")))  +
	theme_bw() +
	my.theme +
	scale_x_continuous(name="Years since first collection") +
	scale_y_continuous(name="Frequency", breaks=c(0, 1, 2))
dev.off()

my.theme <- theme(
legend.justification=c(1, 1),
legend.position=0,
legend.spacing=unit(0, 'cm'),
legend.margin=margin(.0, .0, .0, .0, 'cm'),
legend.box.background=element_blank(),
text=element_text(size=35),
legend.key.size=unit(1.2, 'cm'),
axis.text=element_text(size=30, colour='black'),
axis.text.x=element_text(angle=60, hjust=1),
legend.text=element_text(size=30),
panel.grid=element_blank(),
legend.background=element_blank(),
legend.key=element_rect(fill="#00000000", colour="#00000000")
)


data <- read.csv("VL.cartoon.csv")

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')
pdf("plots/VL.cartoon.pdf")
ggplot(data, aes(x=Date, y=VL)) +
	geom_rect(xmin=2.1, xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR) +
	geom_line(size=1) +
	theme_bw() +
	my.theme +
	coord_cartesian(ylim=c(10, 120000), xlim=c(0, 4)) +
	scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*''))) +
	scale_x_continuous(name="Years since first collection", breaks=c(0, 1, 2, 3, 4)) +
	geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=8) +
	scale_colour_manual(name="", breaks=c(0, 1), limits=c(0, 1), labels=c("Training", "Censored"), values=c('black', 'red')) + scale_shape_manual(name="", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 5))
dev.off()