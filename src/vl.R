pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=3.15, colormodel='rgb')

p <- ggplot(data, aes(x=Date, y=VL))
p <- p + geom_rect(xmin=as.numeric(as.Date("2006-08-01"), origin="1970-01-01"), xmax=Inf, ymin=-Inf, ymax=Inf, fill="#80808010")
p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load", breaks=10^(1:5), labels=sapply(1:5, function(x) bquote(''*10^{.(x)}*'')))
p <- p + scale_x_date(name="Year", breaks=as.Date(paste0(seq(1992, 2017, by=4), "-01-01")), labels=seq(1992, 2017, by=4))
p <- p + coord_cartesian(ylim=c(10, 100000), xlim=as.Date(c("1997-01-01", "2017-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=2.5)
p <- p + theme(
	legend.justification=c(1, 1),
	legend.position=c(.98,.98),
	legend.spacing=unit(.25, 'cm'),
	legend.margin=margin(.05, .05, .05, .05, 'cm'),
	legend.key.size=unit(.25, 'cm'),
	legend.box.background=element_rect(fill="white", colour="white"),	
	text=element_text(size=10),
	axis.text=element_text(size=9),
	panel.grid=element_blank()
)
p <- p + scale_colour_manual(name=c("Data set"), breaks=c(0, 1), limits=c(0, 1), labels=c("Training", "Censored"), values=c('black', 'red'))
p <- p + scale_shape_manual(name="Type", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 18))
p



pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=3.15, colormodel='rgb')

p <- ggplot(data, aes(x=Date, y=VL))
p <- p + geom_rect(xmin=as.numeric(as.Date("2000-08-01"), origin="1970-01-01"), xmax=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), ymin=-Inf, ymax=Inf, fill="#a0a0a010")
p <- p + geom_rect(xmin=as.numeric(as.Date("2006-09-01"), origin="1970-01-01"), xmax=Inf, ymin=-Inf, ymax=Inf, fill="#80808010")
p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load", breaks=10^(1:5), labels=sapply(1:5, function(x) bquote(''*10^{.(x)}*'')))
p <- p + scale_x_date(name="Year", breaks=as.Date(paste0(seq(1996, 2017, by=4), "-01-01")), labels=seq(1996, 2017, by=4))
p <- p + coord_cartesian(ylim=c(10, 1000000), xlim=as.Date(c("1997-01-01", "2017-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=factor(Censored), shape=Type), data=subset(data, Used != ""), size=2.5)
p <- p + theme(
	legend.justification=c(1, 1),
	legend.position=c(.98,.98),
	legend.spacing=unit(.25, 'cm'),
	legend.margin=margin(.05, .05, .05, .05, 'cm'),
	legend.key.size=unit(.25, 'cm'),
	legend.box.background=element_rect(fill="white", colour="white"),	
	text=element_text(size=10),
	axis.text=element_text(size=9),
	panel.grid=element_blank()
)
p <- p + scale_colour_manual(name=c("Data set"), breaks=c(0, 1), limits=c(0, 1), labels=c("Training", "Censored"), values=c('black', 'red'))
p <- p + scale_shape_manual(name="Type", breaks=c("RNA", "DNA"), limits=c("RNA", "DNA"), labels=c("RNA", "DNA"), values=c(16, 18))
p






p <- ggplot(data, aes(x=Estimated.Date, fill=factor(as.Date(Collection.Date, origin="1970-01-01")))) + geom_histogram(bins=10) + scale_fill_manual(name="Collection Date", values=c("#FF9999", "#ff4a4a", "#DD0000", "#a30000"), labels=c("Sep. 2007", "July 2011", "Sep. 2015", "June 2016")) + scale_x_continuous(name="Estimated integration date (year)", breaks=as.numeric(as.Date(paste0(seq(1992, 2017, by=4), "-01-01"))), labels=seq(1992, 2017, by=4)) + scale_y_continuous(name="Frequency", breaks=c(0, 4, 8)) + theme_bw() + theme(
			text=element_text(size=10),
			axis.text=element_text(size=10),
			legend.text=element_text(size=10),
			legend.position=c(.02, .98),
			legend.justification=c(0, 1),
			legend.spacing=unit(0, 'cm'),
			legend.margin=margin(0, 0, 0, 0, 'cm'),
			legend.key.size=unit(.35, 'cm'),
			legend.background=element_blank(),
			legend.box.background=element_blank(),
			panel.grid.major=element_blank(),
		    panel.grid.minor=element_blank()
		)