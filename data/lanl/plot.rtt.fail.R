source("../common/rtt.R")
source("../common/raxml.R")
source("../common/fasttree.R")
library(ape)

extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)

tr <- read.tree("tree/patient_820.tre")

tr <- drop.tip(tr, "REFERENCE")


dates <- extract_dates(tr$tip.label)
types <- extract_types(tr$tip.label)
evol <- node.depth.edgelength(tr)[1:length(tr$tip.label)]

x <- c(min(dates), max(dates))
y <- c(min(evol), max(evol))

a<-0.75
pdf("rtt.pdf", 20*a, 10*a)
par(mfrow=c(1,2))



t.dates <- dates[types=="PLASMA"]
p.dates <- dates[types=="PBMC"]

t.evo <- evol[types=="PLASMA"]
p.evo <- evol[types=="PBMC"]

lm <- glm(t.evo ~ t.dates)

plot(t.dates, 
     t.evo, 
     xlim=x, 
     ylim=y,
     xlab="Time Since Reference (days)", 
     ylab="Expected Number of Subs.", 
     pch=20,  cex=1.2, tck=.01) #,  axes=F)

points(p.dates, p.evo, col="red", pch=5,  cex=1.2)
abline(lm)

mtext("Outgroup Rooting", side=3, adj=0, line=1.1, cex=1.5, font=2); 
legend(1550, 0.030, c("Calibration dates (Plasma)", "Censored dates"), col = c("black", "red"),
        lty = c(-1, -1), pch = c(20, 5),
       merge = TRUE, bg = par("bg"), cex=1.2)

dates[types=="PBMC"] <- NA
tr <- rtt(tr, dates) 


dates <- extract_dates(tr$tip.label)
types <- extract_types(tr$tip.label)
evol <- node.depth.edgelength(tr)[1:length(tr$tip.label)]

x <- c(min(dates), max(dates))
y <- c(min(evol), max(evol))


t.dates <- dates[types=="PLASMA"]
p.dates <- dates[types=="PBMC"]

t.evo <- evol[types=="PLASMA"]
p.evo <- evol[types=="PBMC"]

lm <- glm(t.evo ~ t.dates)

plot(t.dates, 
     t.evo, 
     xlim=x, 
     ylim=y,
     xlab="Time Since Reference (days)", 
     ylab="Expected Number of Subs.", 
     pch=20,  cex=1.2, tck=.01) #,  axes=F)

mtext("RTT Regression Rooting", side=3, adj=0, line=1.1, cex=1.5, font=2); 
points(p.dates, p.evo, col="red", pch=5,  cex=1.2)
abline(lm)

legend(00, 0.045, c("Calibration dates (Plasma)", "Censored dates"), col = c("black", "red"),
        lty = c(-1, -1), pch = c(20, 5),
       merge = TRUE, bg = par("bg"), cex=1.2)

