library(ape)
library(optparse)

args.all <- commandArgs(trailingOnly = F)

if (any(grep("--file=", args.all))) {
	source.dir <- dirname(sub("--file=", "", args.all[grep("--file=", args.all)]))
} else {
	file.arg <- F

	for (i in 1:length(args.all)) {
		if (file.arg) {
			source.dir <- dirname(args.all[i])
		
			break
		}
		
		file.arg <- args.all[i] == '-f'
	}
}

source(file.path(source.dir, 'raxml.R'), chdir=T)

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--raxml", type='character')
op <- add_option(op, "--threads", type='numeric', default=2)
op <- add_option(op, "--seed", type='numeric', default=1989)
op <- add_option(op, "--tmpdir", type='character', default=NULL)
op <- add_option(op, "--settings", type='character', default=NA)
args <- parse_args(op)

settings.file <- args$settings
if (!is.na(settings.file)) {
	settings <- readLines(settings.file)
	settings.filter <- unlist(lapply(op@options, function(x) settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]))
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}

fasta.file <- args$fasta
tree.file <- args$tree
exec <- args$raxml
threads <- args$threads
seed <- args$seed
tmp.dir <- args$tmpdir

tmp <- !is.null(tmp.dir)

tree <- raxml(fasta.file, N=100, threads=threads, executable=exec, parsimony.seed=seed, bootstrap.seed=seed, clear=!tmp, tmp.dir=if(tmp) tmp.dir else tempdir())
write.tree(tree, tree.file)
