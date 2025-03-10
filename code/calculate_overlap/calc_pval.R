## Compute H pval using sliding blocks

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_pval.R <reps> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

inFile = args[1]
outFile = args[2]

# Read in reps of H
df <- fread(inFile)

# Set params
trueH <- as.numeric(df[1,1])
L <- as.numeric(df[1,2])

# Get variance of reps
reps <- df[2:nrow(df), 1]
varH <- var(reps$allH)
meanH <- mean(reps$allH)

# P-value from sims
pvalSim <- sum(reps$allH > trueH) / length(reps$allH)
pvalNorm <- pnorm(trueH, mean = mean(reps$allH), sd=sqrt(varH), lower.tail = F)

# Construct output
dfOut <- data.frame(H = trueH, L = L, meanH = meanH, varH = varH, pvalNorm = pvalNorm, pvalSim = pvalSim)
fwrite(dfOut, outFile, quote = F, row.names = F, sep = "\t")



















