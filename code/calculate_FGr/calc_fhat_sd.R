## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

FGr_file = args[1]
SNP_file = args[2]
out_Fhat = args[3]
out_sd = args[4]

# Read in Fmat
dfMat <- as.matrix(fread(FGr_file))
fhat_raw <- apply(dfMat, 1, sum)
print(length(fhat_raw))

# Scale
dfSNP <- fread(SNP_file)
L <- sum(dfSNP$nSNP)
print(L)
fhat <- fhat_raw / (L -1)

# Get SD
sdF <- sd(fhat)
varF <- var(fhat)

# Construct output
dfOut <- data.frame(sigmaf = sdF, sigmaF2 = varF, L = L)
fwrite(dfOut, out_sd, quote = F, row.names = F, sep = "\t")

# FGR output
fwrite(as.data.frame(fhat), out_Fhat, quote = F, row.names = F, sep = "\t")


