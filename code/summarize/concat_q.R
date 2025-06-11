# This script concatenates all the overlap statistics from a given dataset

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_OverlapStats.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 10)

for (i in 2:length(args)) {

  print(i)

  # Get results file name
  filename = args[i]

  # Extract dataset
  dataset <- strsplit(filename, "/")[[1]][3]

  # Extract which GWAS
  gwas <- strsplit(filename, "/")[[1]][4]

  # Extract covars
  covar <- strsplit(filename, "/")[[1]][5]

  # Extract GWAS model
  gtype <- strsplit(filename, "/")[[1]][6]
  print(gtype)

  # Extract threshold
  threshold <- strsplit(filename, "/")[[1]][7]

  # Extract contrasts
  constrasts <- strsplit(filename, "/")[[1]][8]

  # Extract phenotype
  tmp <- strsplit(filename, "/")[[1]][9]
  phenotype <- strsplit(tmp, ".results")[[1]][1]

  # Read in results
  df <- fread(filename)
  names_from_file <- colnames(df)
  df$dataset <- dataset
  df$gwas <- gwas
  df$contrasts <- constrasts
  df$covar <- covar
  df$gtype <- gtype
  df$threshold <- threshold
  df$phenotype <- phenotype
  colnames(dfOut) <- c(names_from_file, "dataset", "gwas", "contrasts", "covar", "gtype", "threshold", "phenotype")
  print(head(df))
  dfOut <- rbind(dfOut, df)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




