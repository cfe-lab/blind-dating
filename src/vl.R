library(ggplot2)
library(optparse)

op <- OptionParser()
op <- add_option(op, "--vlfile", type='character')
op <- add_option(op, "--patid", type='character')
op <- add_option(op, "--yearstart", type='character', default=NA)
op <- add_option(op, "--yearend", type='character', default=NA)
op <- add_option(op, "--yearby", type='double', default=2)
op <- add_option(op, "--vlstart", type='double', default=10)
op <- add_option(op, "--vlend", type='double', default=NA)
op <- add_option(op, "--therapy", type='character', default=NA)
op <- add_option(op, "--therapy2", type='character', default=NA)
op <- add_option(op, "--two", type='logical', action='store_true', default=F)
args <- parse_args(op)

vl.file <- args$vlfile
pat.id <- args$patid
therapy2 <- args$therapy2
therapy <- args$therapy
year.start <- args$yearstart
year.end <- args$yearend
vl.start <- args$vlstart
vl.end <- args$vlend
two <- args$two

therapy <- as.numeric(as.Date(strsplit(therapy, ' ')[[1]]))

if (length(therapy) %% 2 == 1)
	therapy <- c(therapy, Inf)
	
therapy <- t(matrix(therapy, nrow=2))

if (!is.na(therapy2)) {
	therapy2 <- as.numeric(as.Date(strsplit(therapy2, ' ')[[1]]))

	if (length(therapy2) %% 2 == 1)
		therapy2 <- c(therapy2, Inf)
	
	therapy2 <- t(matrix(therapy2, nrow=2))
}

year.by <- args$yearby

data <- read.csv(vl.file, stringsAsFactors=F)
data$Date <- as.Date(data$Date)
data$Used <- gsub("(V3| & )", "", data$Used)

if (is.na(year.start)) {
	year.start <- as.numeric(format(min(data$Date), "%Y"))
}

if (is.na(year.end)) {
	year.end <- as.numeric(format(max(data$Date), "%Y")) + 1
}

if (is.na(vl.end)) {
	vl.end <- max(data$VL)
}

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=4.2, colormodel='rgb')

THERAPY_COLOUR <- "#a0a0a0"
THERAPY_COLOUR2 <- "#c0c0c0"

p <- ggplot(data, aes(x=Date, y=VL))

for (i in 1:nrow(therapy))
	p <- p + geom_rect(xmin=therapy[i, 1], xmax=therapy[i, 2], ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR)

if (!any(is.na(therapy2))) {
	for (i in 1:nrow(therapy2))
		p <- p + geom_rect(xmin=therapy2[i, 1], xmax=therapy2[i, 2], ymin=-Inf, ymax=Inf, fill=THERAPY_COLOUR2)
}

p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load", breaks=10^c(1, 3, 5), labels=sapply(c(1, 3, 5), function(x) bquote(''*10^{.(x)}*'')))
p <- p + scale_x_date(name="Collection Year", breaks=as.Date(paste0(seq(year.start, year.end, by=year.by), "-01-01")), labels=seq(year.start, year.end, by=year.by))
p <- p + coord_cartesian(ylim=c(vl.start, vl.end), xlim=as.Date(paste0(c(year.start, year.end), "-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=6)

if (two)
	p <- p + geom_point(aes(y=VL * 4, shape=Type), data=subset(data, Used != "" & Date > as.Date("2008-01-01")), colour="#FF8800", size=6)

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

if (two) {
	p <- p + scale_colour_manual(name="", breaks=c(0, -1, 1), limits=c(0, -1, 1), labels=c("Training", "Training 2", "Censored"), values=c('black', 'darkblue', 'red'))
} else {
	p <- p + scale_colour_manual(name="", breaks=c(0, 1), limits=c(0, 1), labels=c("Training", "Censored"), values=c('black', 'red'))
}
	
p <- p + scale_shape_manual(name="", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 5))
p <- p + guides(shape=guide_legend(order=2), colour=guide_legend(override.aes=list(shape=15, size=8), order=1))

pdf(paste0("VL.", pat.id, ".pdf"))
p
dev.off()