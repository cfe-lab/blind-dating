p <- ggplot(data, aes(x=Date, y=VL))
p <- p + geom_rect(xmin=as.numeric(as.Date("2006-08-01"), origin="1970-01-01"), xmax=Inf, ymin=-Inf, ymax=Inf, fill="#aaaaaa")
p <- p + geom_line(size=.5) 
p <- p + scale_y_log10(name="Viral load")
p <- p + scale_x_date(name="Year")
p <- p + coord_cartesian(ylim=c(10, 100000), xlim=as.Date(c("1996-01-01", "2017-01-01")))
p <- p + theme_bw()
p <- p + geom_point(aes(colour=Used, shape=Type), data=subset(data, Used != ""), size=2.5)
p <- p + theme(
	legend.justification=c(1, 1),
	legend.position=c(.98,.98),
	legend.spacing=unit(.25, 'cm'),
	legend.margin=margin(.05, .05, .05, .05, 'cm'),
	legend.key.size=unit(.25, 'cm'),
	legend.box.background=element_rect(fill="white", colour="white"),	
	text=element_text(size=10),
	panel.grid=element_blank()
)
p <- p + scale_colour_manual(name=c("Data set"), breaks=c("NEF", "V3", "V3 & NEF"), limits=c("NEF", "V3", "V3 & NEF"), labels=c("NEF", "V3", "Both"), values=c("#bb6600", "#1111ff", "#336633"))
p <- p + scale_shape_manual(name="Type", breaks=c("PLASMA", "PBMC","WHOLE BLOOD"), limits=c("PLASMA", "PBMC","WHOLE BLOOD"), labels=c("Plasma RNA", "PBMC DNA", "Whole Blood Draw DNA"), values=c(16, 17, 15))
p