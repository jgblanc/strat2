# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# get all contrasts
sdi_eur <- strsplit(strsplit(args[1], "/")[[1]][5], "_")[[1]][1]
out_prefix <- args[2]

# Go through all pairwise comparisons

# SDI-EUR
nA <- 662
nB <- 27
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, sdi_eur, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, sdi_eur, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)


