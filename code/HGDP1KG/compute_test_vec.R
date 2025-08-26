# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# get all contrasts
eas_nfe <- strsplit(strsplit(args[1], "/")[[1]][5], "_")[[1]][1]
eas_sas <- strsplit(strsplit(args[2], "/")[[1]][5], "_")[[1]][1]
eas_afr <- strsplit(strsplit(args[3], "/")[[1]][5], "_")[[1]][1]
eas_amr <- strsplit(strsplit(args[4], "/")[[1]][5], "_")[[1]][1]
nfe_sas <- strsplit(strsplit(args[5], "/")[[1]][5], "_")[[1]][1]
nfe_afr <- strsplit(strsplit(args[6], "/")[[1]][5], "_")[[1]][1]
nfe_amr <- strsplit(strsplit(args[7], "/")[[1]][5], "_")[[1]][1]
sas_afr <- strsplit(strsplit(args[8], "/")[[1]][5], "_")[[1]][1]
sas_amr <- strsplit(strsplit(args[9], "/")[[1]][5], "_")[[1]][1]
afr_amr <- strsplit(strsplit(args[10], "/")[[1]][5], "_")[[1]][1]
out_prefix <- args[11]

# Go through all pairwise comparisons

# EAS-NFE
nA <- 825
nB <- 689
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eas_nfe, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eas_nfe, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-SAS
nA <- 825
nB <- 790
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eas_sas, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eas_sas, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AFR
nA <- 825
nB <- 1003
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eas_afr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eas_afr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AMR
nA <- 825
nB <- 552
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, eas_amr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, eas_amr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-SAS
nA <- 689
nB <- 790
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, nfe_sas, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, nfe_sas, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AFR
nA <- 689
nB <- 1003
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, nfe_afr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, nfe_afr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AMR
nA <- 689
nB <- 552
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, nfe_amr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, nfe_amr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AFR
nA <- 790
nB <- 1003
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, sas_afr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, sas_afr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AMR
nA <- 790
nB <- 552
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, sas_amr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, sas_amr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)

# AFR-AMR
nA <- 1003
nB <- 552
N <- nA+nB
tvec <- c(rep(N/(2*nA), nA) , rep(-1 * (N/(2*nB))))
varTvec <- var(tvec)
fwrite(as.data.frame(tvec),paste0(out_prefix, afr_amr, ".txt") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(varTvec),paste0(out_prefix, afr_amr, "_Var.txt") ,row.names=F,quote=F,sep="\t", col.names = T)



















