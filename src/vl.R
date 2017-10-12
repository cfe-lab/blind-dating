library(ggplot2)

data <- read.csv("VL.***REMOVED***.csv", stringsAsFactors=F)
data$Date <- as.Date(data$Date)
data$Used <- gsub("(V3| & )", "", data$Used)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')

THERAPY_COLOUR <- "#a0a0a0"
THERAPY_COLOUR2 <- "#c0c0c0"

p <- ggplot(data, aes(x=Date, y=VL))
p <- p + geom_rect(xmin=as.numeric(as.Date("1997-08-01"), origin="1970-01-01"), xmax=as.numeric(as.Date("1997-11-01"), origin="1970-01-01"), ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR)
p <- p + geom_rect(xmin=as.numeric(as.Date("2006-08-01"), origin="1970-01-01"), xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR)
p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*'')))
p <- p + scale_x_date(name="Collection Year", breaks=as.Date(paste0(seq(1996, 2017, by=2), "-01-01")), labels=seq(1996, 2017, by=2))
p <- p + coord_cartesian(ylim=c(10, 100000), xlim=as.Date(c("1996-01-01", "2017-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=6)
p <- p + theme(
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
p <- p + scale_colour_manual(name="", breaks=c(0, 1), limits=c(0, 1), labels=c("Training", "Censored"), values=c('black', 'red'))
p <- p + scale_shape_manual(name="", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 5))
p <- p + guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1))

pdf("VL.***REMOVED***.pdf")
p
dev.off()

library(ggplot2)

data <- read.csv("VL.***REMOVED***.csv", stringsAsFactors=F)
data$Date <- as.Date(data$Date)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')

THERAPY_COLOUR <- "#a0a0a0"
THERAPY_COLOUR2 <- "#c0c0c0"

p <- ggplot(data, aes(x=Date, y=VL))
p <- p + geom_rect(xmin=as.numeric(as.Date("2000-08-01"), origin="1970-01-01"), xmax=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR2)
p <- p + geom_rect(xmin=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), xmax=Inf, ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR)
p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*'')))
p <- p + scale_x_date(name="Collection Year", breaks=as.Date(paste0(seq(1997, 2017, by=2), "-01-01")), labels=seq(1997, 2017, by=2))
p <- p + coord_cartesian(ylim=c(10, 1000000), xlim=as.Date(c("1997-01-01", "2017-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=6)
p <- p + theme(
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
p <- p + scale_colour_manual(name="", breaks=c(0, -1, 1), limits=c(0, -1, 1), labels=c("Training", "Training 2", "Censored"), values=c('black', 'darkblue', 'red'))
p <- p + scale_shape_manual(name="", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 5))
p <- p + guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1))

pdf("VL.***REMOVED***.pdf")
p
dev.off()



library(ggplot2)

THERAPY_COLOUR <- "#a0a0a0"
THERAPY_COLOUR2 <- "#c0c0c0"

#data <- read.csv("stats/patient_***REMOVED***.data.csv", stringsAsFactors=F)
data <- read.csv("stats/patient_***REMOVED***-with_ref.data.csv", stringsAsFactors=F)

data <- subset(data, Censored == 1)

m <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(min(data$Estimated.Date), origin="1970-01-01")))
M <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(max(data$Estimated.Date), origin="1970-01-01"))) + 1

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=5, colormodel='rgb')

p <- ggplot(data, aes(x=Estimated.Date, fill=factor(as.Date(Collection.Date, origin="1970-01-01"), levels=c("2007-09-20", "2011-07-19", "2015-09-30", "2016-06-29")))) +
	geom_rect(xmin=as.numeric(as.Date("1997-08-01"), origin="1970-01-01"), xmax=as.numeric(as.Date("1997-11-01"), origin="1970-01-01"), ymin=0, ymax=7, fill=THERAPY_COLOUR) +
	geom_rect(xmin=as.numeric(as.Date("2006-08-01"), origin="1970-01-01"), xmax=Inf, ymin=0, ymax=7, fill=THERAPY_COLOUR) +
	geom_segment(x=as.numeric(as.Date("1996-08-12"), origin="1970-01-01"), xend=as.numeric(as.Date("1996-08-12"), origin="1970-01-01"), y=7, yend=5, arrow=arrow(length = unit(0.5, "cm"))) +
	geom_histogram(breaks=as.numeric(as.Date(paste0(seq(m, M), "-01-01")))) +
	scale_fill_manual(name="Collection Date", values=c("#FF9999", "#ff4a4a", "#DD0000", "#a30000"), labels=c("Sep. 2007", "Jul. 2011", "Sep. 2015", "Jun. 2016")) +
	scale_x_continuous(name="Estimated integration year", breaks=as.numeric(as.Date(paste0(seq(m, M), "-01-01"))), labels=seq(m, M)) +
	scale_y_continuous(name="Frequency", breaks=c(0, 2, 4, 6, 8), limits=c(0, 12)) +
	guides(fill=guide_legend(override.aes=list(size=8, colour="#00000000"), ncol=2)) +
	theme_bw() +
	theme(
		text=element_text(size=35),
		axis.text=element_text(size=30, colour='black'),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
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
	
#pdf("plots/patient_***REMOVED***.hist.pdf")
pdf("plots/patient_***REMOVED***-with_ref.hist.pdf")

p
dev.off()


library(ggplot2)

THERAPY_COLOUR <- "#a0a0a0"
THERAPY_COLOUR2 <- "#c0c0c0"

#data <- read.csv("stats/patient_***REMOVED***.comb.data.csv", stringsAsFactors=F)
data <- read.csv("stats/patient_***REMOVED***-with_ref.comb.data.csv", stringsAsFactors=F)

data <- subset(data, Censored > 0)

m <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(min(data$Estimated.Date), origin="1970-01-01")))
M <- as.numeric(gsub("(.+)-.+-.+", "\\1", as.Date(max(data$Estimated.Date), origin="1970-01-01"))) + 1

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=5, colormodel='rgb')

p <- ggplot(data, aes(x=Estimated.Date, fill=factor(as.Date(Collection.Date, origin="1970-01-01"), levels=c("2013-03-15", "2016-08-15")))) +
	geom_rect(xmin=as.numeric(as.Date("2000-08-01"), origin="1970-01-01"), xmax=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), ymin=0, ymax=4, fill=THERAPY_COLOUR2) +
	geom_rect(xmin=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), xmax=Inf, ymin=0, ymax=4, fill=THERAPY_COLOUR) +
	geom_segment(x=as.numeric(as.Date("1997-02-11"), origin="1970-01-01"), xend=as.numeric(as.Date("1997-02-11"), origin="1970-01-01"), y=4, yend=3, arrow=arrow(length = unit(0.5, "cm"))) +
	geom_histogram(breaks=as.numeric(as.Date(paste0(seq(m, M), "-01-01")))) +
	scale_fill_manual(name="Collection Date", values=c("#ff4a4a", "#a30000"), labels=c("Mar. 2013", "Aug. 2016")) +
	scale_x_continuous(name="Estimated integration year", breaks=as.numeric(as.Date(paste0(seq(m, M), "-01-01"))), labels=seq(m, M)) +
	scale_y_continuous(name="Frequency", breaks=c(0, 1, 2, 3, 4), limits=c(0, 6.5)) +
	guides(fill=guide_legend(override.aes=list(size=8, colour="#00000000"))) +
	theme_bw() +
	theme(
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
	
#pdf("plots/patient_***REMOVED***.hist.pdf")
pdf("plots/patient_***REMOVED***-with_ref.hist.pdf")

p
dev.off()

