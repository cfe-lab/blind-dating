#!/usr/bin/Rscript

library(seqinr)
library(optparse)

echo <- function(...) if (verbose) cat(..., '\n', sep="")
mylapply <- function(...) if (cpu.cores == 1) lapply(...) else mclapply(..., mc.cores=cpu.cores)

do.remove.missing <- function(f, cutoff) {
	f[sapply(f, function(x) sum(! x %in% c('a', 'c', 'g', 't', 'A', 'C', 'G', 'T')) / length(x) > cutoff)]
}

do.remove.stop.codons <- function(f) {
	f.a <- mylapply(f, function(x) translate(gsub('-', 'n', f), ambiguous=T))
	f[sapply(f.a, function(x) all(x != '*'))]
}

do.remove.duplicates <- function(f) {
	f.dup <- mylapply(f, function(x) which(sapply(f, function(y) all(x == y))))
	f.dup 
}

op <- OptionParser()
op <- add_option(op, c("-f", "--fasta"), type='character')
op <- add_option(op, c("-c", "--cleaned-fasta"), type='character')
op <- add_option(op, c("-s", "--keep-stop-codons"), type='logical', action='store_false', default=T) # stores reverse
op <- add_option(op, c("-d", "--keep-duplicates"), type='logical', action='store_false', default=T) # stores reverse
op <- add_option(op, c("-n", "--missing-cutoff"), type='numeric', default=0.05)
op <- add_option(op, c("-t", "--cpu-cores", type='numeric', default=getOption('mc.cores', 1))
op <- add_option(op, c("-v", "--verbose"), type='logical', action='store_true', default=F)
args <- parse_args(op)

fasta.name <- op[["fasta"]]
cleaned.fasta.name <- op[["cleaned-fasta"]]
remove.stop.codons <- op[["keep-stop-codons"]] # reversed
remove.duplicates <- op[["keep-duplictes"]] # reversed
missing.cutoff <- op[["missing-ctuoff"]]
cpu.cores <- op[["cpu-cores"]]
verbose <- op[["verbose"]]

echo("clean_fasta.R")
echo("Params:")
echo("fasta: ", fasta.name)
echo("cleaned-fasta: ", cleaned.fasta.name)
echo("keep-stop-codons: ", !remove.stop.codons)
echo("keep-duplicates: ", !remove.duplicates)
echo("missing-cutoff: ", missing.cutoff)
echo("cpu-cores: ", cpu.cores)
echo("verbose: ", verbose)
echo()

if (cpu.cores != 1)
	library(parallel)

f <- read.fasta(fasta.name)

f.cleaned <- do.remove.missing(f, missing.cutoff)

if (remove.stop.codons)
	f.cleaned <- do.remove.stop.codons(f.cleaned)

if (remove.duplicates)
	f.cleaned <- do.remove.duplicates(f.cleaned)
	
write.fasta(f.cleaned, names(f.cleaned), cleaned.fasta.name)