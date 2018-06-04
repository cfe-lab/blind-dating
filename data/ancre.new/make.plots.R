#!/usr/bin/Rscript

library(ggplot2)

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

my.theme <- theme(
			text=element_text(size=35),
			axis.text=element_text(size=25, colour='black'),
			legend.text=element_text(size=25),
			legend.position=c(0.01, 1),
			axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
			legend.justification=c(0, 1),
			legend.spacing=unit(0, 'cm'),
			legend.margin=margin(0, 0, 0, 0, 'cm'),
			legend.key.size=unit(1, 'cm'),
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

pdf.options(family="Helvetica", fonts="Helvetica", width=7, height=7, colormodel='rgb')

ancre.colours <- c("#0072B2", "#D55E00", "#999999")

files <- paste0("patient_", c(2658, 825, 7259, 7265, 13333, 13334, 13336, "cp-1362"), ".cens.csv")

data.files.rtt <- paste0("stats/", gsub(".csv", ".data.csv", files))
stats.files.rtt <- paste0("stats/", gsub(".csv", ".stats.csv", files))
data.rtt <- lapply(data.files.rtt, read.csv)
stats.rtt <- do.call(rbind, lapply(stats.files.rtt, read.csv))

good.rtt <- with(stats.rtt, Model.Fit == 1)
data.rtt <- data.rtt[good.rtt]

data.rtt <- lapply(data.rtt, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); subset(x, x$Censored == 1)})

g <- ggplot() + theme_bw() + my.theme2

i <- 1
for (x in data.rtt) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill="black", alpha=1/length(data.rtt), color="#00000000")
	
	pdf(paste0("plots/", gsub(".csv", ".dens.pdf", files[i])))
	print(ggplot() + theme_bw() + my.theme2 +  geom_density(aes(x=Scaled.Difference), data=x, fill='red', color="#00000000") + geom_vline(xintercept=0, linetype='dashed', colour='grey') + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 3)))
	dev.off()

	i <- i + 1
}

pdf("ancre.density.rtt.pdf")
g + geom_vline(xintercept=0, linetype='dashed', size=1, colour='black') + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 6))
dev.off()

data.files.ogr <- paste0("stats/", gsub(".csv", "-with_ref.data.csv", files))
stats.files.ogr <- paste0("stats/", gsub(".csv", "-with_ref.stats.csv", files))
data.ogr <- lapply(data.files.ogr, read.csv)
stats.ogr <- do.call(rbind, lapply(stats.files.ogr, read.csv))

good.ogr <- with(stats.ogr, Model.Fit == 1)
data.ogr <- data.ogr[good.ogr]

data.ogr <- lapply(data.ogr, function(x) {x$Scaled.Difference <- x$Date.Difference / (max(x$Collection.Date) - min(x$Collection.Date)); subset(x, x$Censored == 1)})

g <- ggplot() + theme_bw() + my.theme2 + geom_vline(xintercept=0, linetype='dashed', colour='grey')
i <- 1
for (x in data.ogr) {
	g <- g + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', alpha=1/length(data.ogr), color="#00000000")
	
	pdf(paste0("plots/", gsub(".csv", "-with_ref.dens.pdf", files[i])))
	print(ggplot() + theme_bw() + my.theme2 + geom_density(aes(x=Scaled.Difference), data=x, fill='purple', color="#00000000") + geom_vline(xintercept=0, linetype='dashed', colour='grey') + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 3)))
	dev.off()
	
	i <- i + 1
}

pdf("ancre.density.ogr.pdf")
g + scale_x_continuous(name="Scaled Date Difference", limits=c(-1.6, 1.6)) + scale_y_continuous(name="Density", limits=c(0, 6))
dev.off()

data.rtt <- lapply(data.files.rtt, read.csv)
data.rtt <- data.rtt[good.ogr & good.rtt]
data.ogr <- lapply(data.files.ogr, read.csv)
data.ogr <- data.ogr[good.ogr & good.rtt]

data.all <- do.call(rbind, lapply(1:length(data.rtt), function(i) cbind(Patient=stats.rtt[good.ogr & good.rtt, "Patient"][i], merge(data.rtt[[i]], data.ogr[[i]], by=c("ID", "Type", "Censored", "Collection.Date"), suffixes=c(".rtt", ".ogr")))))
data.all <- subset(data.all, Censored == 1)

data.all$Patient <- as.factor(data.all$Patient)
pat_levels <- levels(data.all$Patient)

pdf("ancre.comp.pdf")
p <- ggplot(data.all) + geom_abline(lty=2) + geom_point(aes(x=Estimated.Date.rtt / 365.25, y=Estimated.Date.ogr / 365.25, colour=Patient), size=3) + annotate("text", x=6, y=-1, label=sprintf("Concordance: %.2f", with(data.all, concord(Estimated.Date.rtt, Estimated.Date.ogr))), vjust=0, hjust=1, size=8) + scale_colour_manual(name="", values=ancre.colours, label=gsub("^.+_", "", gsub(".cens", "", pat_levels))) + scale_x_continuous(name="Estimated Date (RTT)", limits=c(-1, 6)) + scale_y_continuous(name="Estimated Date (OGR)", limits=c(-1, 6)) + theme_bw() + my.theme + guides(colour=guide_legend(override.aes=list(size=6)))
saveRDS(p, "ancre.comp.rds")
p
dev.off()

cat("Concordance:\n")
sup <- lapply(unique(data.all$Patient), function(x) cat(x, "; ", with(subset(data.all, Patient == x), concord(Estimated.Date.rtt, Estimated.Date.ogr)), "\n", sep=""))
cat("Total: ", with(data.all, concord(Estimated.Date.rtt, Estimated.Date.ogr)), "\n", sep="")