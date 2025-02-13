## This script calculates the empirical variance of the GWAS panel genotypes

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_GWAS_variance.R <dosage counts> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

inFile = args[1]
outFile = args[2]


# Initialize an empty list to store results
results_list <- list()
index <- 1

con <- file(inFile, open = "r")

# Read and discard the first line to skip it
readLines(con, n = 1, warn = FALSE)

# Read and process each line
while(TRUE) {

  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break  # Exit loop if end of file

  # Get fields
  fields <- strsplit(line, "\\s+")[[1]]
  dosages <- as.numeric(fields[7:length(fields)])
  ID <- fields[2]

  # Convert to dosage of ALT allele
  dosages <- 2 - dosages

  # Calculate variance
  variance <- var(dosages)

  # Save results to list
  results_list[[index]] <- data.table(ID = ID, Var = variance)
  index <- 1+ index
}

# Close the connection
close(con)

# Combine the list into a data.table and save the output
dfOut <- rbindlist(results_list)
fwrite(dfOut, outFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
