# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# get all contrasts
eng_ni <- strsplit(strsplit(args[1], "/")[[1]][5], "_")[[1]][1]
eng_roi <- strsplit(strsplit(args[2], "/")[[1]][5], "_")[[1]][1]
eng_sct <- strsplit(strsplit(args[3], "/")[[1]][5], "_")[[1]][1]
eng_wal <- strsplit(strsplit(args[4], "/")[[1]][5], "_")[[1]][1]
ni_roi <- strsplit(strsplit(args[5], "/")[[1]][5], "_")[[1]][1]
ni_sct <- strsplit(strsplit(args[6], "/")[[1]][5], "_")[[1]][1]
ni_wal <- strsplit(strsplit(args[7], "/")[[1]][5], "_")[[1]][1]
roi_sct <- strsplit(strsplit(args[8], "/")[[1]][5], "_")[[1]][1]
roi_wal <- strsplit(strsplit(args[9], "/")[[1]][5], "_")[[1]][1]
sct_wal <- strsplit(strsplit(args[10], "/")[[1]][5], "_")[[1]][1]
out_prefix <- args[11]

# Go through all pairwise comparisons

# ENG-NI
nA <- 142215
nB <- 1120
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eng_ni, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eng_ni, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# ENG-ROI
nA <- 142215
nB <- 1802
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eng_roi, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eng_roi, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# ENG-SCT
nA <- 142215
nB <- 14591
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eng_sct, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eng_sct, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# ENG-WAL
nA <- 142215
nB <- 8125
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eng_wal, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eng_wal, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NI-ROI
nA <- 1120
nB <- 1802
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, ni_roi, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, ni_roi, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NI-SCT
nA <- 1120
nB <- 14591
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, ni_sct, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, ni_sct, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NI-WAL
nA <- 1120
nB <- 8125
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, ni_wal, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, ni_wal, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# ROI-SCT
nA <- 1802
nB <- 14591
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, roi_sct, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, roi_sct, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# ROI-WAL
nA <- 1802
nB <- 8125
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, roi_wal, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, roi_wal, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# SCT-WAL
nA <- 14591
nB <- 8125
N <- nA + nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, sct_wal, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, sct_wal, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)



















