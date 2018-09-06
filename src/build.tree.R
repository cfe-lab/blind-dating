library(ape)
library(optparse)

bd.src <- Sys.getenv("BDSRC", ".")
source(file.path(bd.src, "raxml.R"), chdir=T)

get.val <- function(x, default) if (is.null(x))	default	else x

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--tree", type='character')
op <- add_option(op, "--info", type='character')
op <- add_option(op, "--raxml", type='character')
op <- add_option(op, "--threads", type='numeric')
op <- add_option(op, "--seed", type='numeric')
op <- add_option(op, "--tmpdir", type='character')
op <- add_option(op, "--model", type='character')
op <- add_option(op, "--settings", type='character', default=NA)
args <- parse_args(op)

settings.file <- args$settings
if (!is.na(settings.file)) {
	settings <- readLines(settings.file)
	settings.filter <- unlist(
		lapply(op@options, function(x)
			settings[grepl(paste0("^", x@long_flag, "(=|$)"), settings)]
		)
	)
	args.settings <- parse_args(op, args=settings.filter)
	args <- c(args, args.settings)
}

fasta.file <- args$fasta
tree.file <- args$tree
exec <- get.val(args$raxml, "raxml")
threads <- get.val(args$threads, 2)
seed <- get.val(args$seed, 1989)
model <- get.val(args$model, "")
tmp.dir <- args$tmpdir

tmp <- !is.null(tmp.dir)

tree <- raxml(
	fasta.file,
	N=100,
	threads=threads,
	executable=exec,
	parsimony.seed=seed,
	bootstrap.seed=seed,
	model=model clear=!tmp,
	tmp.dir=if(tmp) tmp.dir else tempdir()
)
write.tree(tree, tree.file)
