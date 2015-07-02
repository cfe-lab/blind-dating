library(ape)
library(parallel)

#source('ape.patches.R')
source('include/rtt.R')
source('include/test.R')
source('include/raxml.R')
source('include/queue.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))

n.simulated <- 50

trees.read <- function(base.path) {
	ml.tree.read <- function(i, path){
		read.tree(sprintf("%s/HIV_ml_%d_out.nwk", path, i))
	}
	mclapply(1:n.simulated, ml.tree.read, base.path)
}

hiv.rna.read <- function(){
	trs <- dir('./7_ml_tree')[grep('_rna.nwk', dir('./7_ml_tree'))]
	ml.tree.read <- function(tr) {
		read.tree(paste('./7_ml_tree', tr, sep='/'))
	}
	lapply(trs, ml.tree.read)
}

extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))
extract_tag <- function(x) (gsub("(.+)_([0-9\\.]+)$", "\\1", x, perl=T))

extract_plasm_tag <- function(x) (gsub("(.+)_PLASMA_([0-9\\.]+)$", "\\1", x, perl=T))
types <- function(x) gsub("(.+)_((PLASMA)|(PBMC))_([0-9\\.]+)'?$", "\\2", x, perl=T)
extract_dates_pp <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)'?$", "\\2", x, perl=T))

#trees <- mclapply(1:n.simulated, ml.tree, mc.cores=1)
#trees <- hiv.rna.read()
trees <- trees.read("trees/")
n.simulated <- length(trees)
run_name <- "HIV RNA"
species <- "Simulated"
n.runs <- 10

#for(remove in c(1, 5, 10, 25)){
pdf(sprintf('hist.pdf', run_name, remove), width=11.5, height=8.5)
par(mfrow=c(1, 2), pty="s")

## Show an example of one of these simulated trees
tree <- trees[[25]]

for(j in 1:length(tree$tip.label)) {
	name <- tree$tip.label[j]
	date <- extract_dates(name)
	tag <- extract_tag(name)
	tree$tip.label[j] <- sprintf("%s_PLASMA_%s", tag, date)

}
tree <- latentize.tree(tree, runif(5) * 0.8 + 0.1, function(name){
	sprintf("%s_PBMC_%s", extract_plasm_tag(name), extract_dates(name))
});

tip.dates <- extract_dates_pp(tree$tip.label)
tip.types <- types(tree$tip.label)

# Mark the PBMC and PLASMA cells
tip.pbmc <- tip.types == "PBMC"
tip.plasma <- tip.types == "PLASMA"

plasma.dates <- tip.dates
plasma.dates[tip.pbmc] <- NA

tree <- rtt(tree, plasma.dates)

plasma.dates <- tip.dates[tip.plasma]
pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates

distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
plasma.dists <- distances[tip.plasma]
pbmc.dists <- distances[tip.pbmc]

model <- glm(plasma.dists ~ plasma.dates)

a<-model$coefficients[[1]]
b<-model$coefficients[[2]]

(pbmc.dists/b - a/b)

plot(
	plasma.dates, plasma.dists, 
	xlab="Time", 
	ylab="Expected Number of Subs.")
mtext("B) Latent Simulated Data", side=3, adj=0, line=1.1, cex=1.5, font=2); 
points(pbmc.s.dates, pbmc.dists, col="red")
abline(model)
##

plot(c(1,1), xlim=c(-0.1,0.2), ylim=c(0, 30), xlab="Normalized Error", ylab="Density")
add <- F
errs <- c()
dens <- queue(TRUE)

for(remove in c(50)) {
	for(i in 1:n.simulated){

		for(k in 1:n.runs) {
			tree <- trees[[i]]
		# plot(tree, cex=0.25)
			for(j in 1:length(tree$tip.label)) {
				name <- tree$tip.label[j]

				date <- extract_dates(name)
				tag <- extract_tag(name)
				tree$tip.label[j] <- sprintf("%s_PLASMA_%s", tag, date)

			}
			tree <- latentize.tree(tree, runif(5) * 0.8 + 0.1, function(name){
				sprintf("%s_PBMC_%s", extract_plasm_tag(name), extract_dates(name))
			});
			# plot(tree, cex=0.25);

			#####
			#
			#####
			tip.dates <- extract_dates_pp(tree$tip.label)
			tip.types <- types(tree$tip.label)

			# Mark the PBMC and PLASMA cells
			tip.pbmc <- tip.types == "PBMC"
			tip.plasma <- tip.types == "PLASMA"

			print(tip.pbmc)
			if(length(which(tip.pbmc==F)) <= 2){
				next
			}

			min_date <- min(tip.dates)
			max_date <- max(tip.dates)

			tip.dates <- (tip.dates - min_date)/(max_date - min_date)


			plasma.dates <- tip.dates
			plasma.dates[tip.pbmc] <- NA
			
			tree <- rtt(tree, plasma.dates)

			plasma.dates <- tip.dates[tip.plasma]
			pbmc.s.dates <- tip.dates[tip.pbmc] #sampled dates

			distances <- node.depth.edgelength(tree)[1:length(tip.dates)]
			plasma.dists <- distances[tip.plasma]
			pbmc.dists <- distances[tip.pbmc]

			model <- glm(plasma.dists ~ plasma.dates)

			a<-model$coefficients[[1]]
			b<-model$coefficients[[2]]

			err <- pbmc.s.dates - (pbmc.dists/b - a/b)
			d <- density(err, bw=0.015)
			polygon(d, col=rgb(0.5, 0.0, 0.9, 1/(2*n.simulated)), xlim=c(-0.1,0.2),  border=rgb(1,1,1,0))
			errs <- c(err, errs)
		}
	}
}

abline(h = 0, col = "black", lwd = 1)
abline(v = mean(errs), col = "black", lwd = 2) 
abline(v = median(errs), col = "red", lwd = 2) 
warnings()
dev.off()


