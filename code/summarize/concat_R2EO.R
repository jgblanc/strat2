# This script concatenates all the overlap statistics from a given dataset

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_OverlapStats.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 19)

for (i in 2:length(args)) {

  print(i)

  # Get results file name
  filename = args[i]

  # Extract dataset
  dataset <- strsplit(filename, "/")[[1]][3]

  # Extract which GWAS
  gwas <- strsplit(filename, "/")[[1]][4]

  # Extract constrasts
  tmp <- strsplit(filename, "/")[[1]][6]
  constrasts <- strsplit(tmp, ".txt")[[1]][1]

  # Read in results
  df <- fread(filename)
  names_from_file <- colnames(df)
  print(names_from_file)
  df$dataset <- dataset
  df$gwas <- gwas
  df$contrasts <- constrasts
  colnames(dfOut) <- c(names_from_file, "dataset", "gwas", "contrasts")
  dfOut <- rbind(dfOut, df)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




