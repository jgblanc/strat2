suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

args = commandArgs(TRUE)

if(length(args) < 3){
  stop("Usage: Rscript residualize_one_by_one.R <snp_file> <pc_file> <out_file>")
}

snp_file = args[1]
pc_file = args[2]
out_file = args[3]

# Read in PCs
dfPCs <- fread(pc_file)
colnames(dfPCs)[1] <- "FID"

# Extract PCs as covariates
covars <- as.matrix(dfPCs %>% select(starts_with("PC")))
print("got covars")

# Initialize an empty list to store results
results_list <- list()
index <- 1

con <- file(snp_file, open = "r")

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



