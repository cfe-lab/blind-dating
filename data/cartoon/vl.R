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
THERAPY <- 4.4666666


THERAPY_COLOUR <- "#111111"
THERAPY_LTY <- 3
THERAPY_COLOUR2 <- "#111111"
THERAPY_LTY2 <- 6

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

data$my.colour <- as.numeric(data$Date)
data$my.colour[data$Censored > 0] <- "censored"
my.colour.break <- unique(data$my.colour)
my.colour.filter <- suppressWarnings(!is.na(as.numeric(my.colour.break)))
my.colour.value <- rep('black', length(my.colour.break))
my.colour.scale <- (as.numeric(my.colour.break[my.colour.filter]) - 0) / (4.166666666 - 0)
my.colour.value[my.colour.filter] <- hsv(my.colour.scale * .75, 0.5, 0.5)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')
pdf("plots/VL.cartoon.pdf")
ggplot(data, aes(x=Date, y=VL)) +
	geom_vline(xintercept=as.numeric(THERAPY), linetype=THERAPY_LTY, colour=THERAPY_COLOUR) +
	geom_line(size=1) +
	theme_bw() +
	my.theme +
	coord_cartesian(ylim=c(10, 120000), xlim=c(0, 8)) +
	scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*''))) +
	scale_x_continuous(name="Years since first collection", breaks=seq(0, 8)) +
	geom_point(aes(colour=factor(my.colour), shape=Type), data=subset(data, Used != ""), size=8) +
	scale_colour_manual(name="", breaks=my.colour.break, limits=my.colour.break, values=my.colour.value) + scale_shape_manual(name="", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 5))
dev.off()
