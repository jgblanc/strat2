suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

args = commandArgs(TRUE)

if(length(args) < 4){
  stop("Usage: Rscript residualize_one_by_one.R <snp_file> <pc_file> <out_file>")
}

snp_file = args[1]
pc_prefix = args[2]
chr_num = as.numeric(args[3])
out_file = args[4]

# Read in PCs for correct chromosomes
if (chr_num %in% c(1,3,5,7,9,11,13,15,17,19,21)) {
  print("Odd Chr")
  pca_file_path <- paste0(pc_prefix, "even_PCA.eigenvec")
  dfPCs <- fread(pca_file_path)
  colnames(dfPCs)[1] <- "FID"
  ind_ids <- dfPCs$FID

} else if (chr_num %in% c(2,4,6,8,10,12,14,16,18,20,22)) {
  print("Even chr")
  pca_file_path <- paste0(pc_prefix, "odd_PCA.eigenvec")
  dfPCs <- fread(pca_file_path)
  colnames(dfPCs)[1] <- "FID"
  ind_ids <- dfPCs$FID

}

# Extract PCs as covariates
covars <- as.matrix(dfPCs %>% select(starts_with("PC")))
#covars <- cbind(1, covars)
print("got covars")

# Initialize an empty list to store results
results_list <- list()
index <- 1

con <- file(snp_file, open = "r")

# Read and discard the first line to skip it
readLines(con, n = 1, warn = FALSE)

# Read and process each line
while(index <= 100) {

  print(index)
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break  # Exit loop if end of file

  # Get fields
  fields <- strsplit(line, "\\s+")[[1]]
  dosages <- matrix(as.numeric(fields[7:length(fields)]), ncol = 1)
  ID <- fields[2]

  # Convert to dosage of ALT allele
  #resids <- resid(lm.fit(x=covars, y=dosages))
  resids <- resid(lm(dosages ~ covars))

  # Save results to list
  results_list[[index]] <- data.table(ID = ID, t(resids))
  index <- 1+ index
}

# Close the connection
close(con)

# Combine the list into a data.table and save the output
dfOut <- rbindlist(results_list)
colnames(dfOut) <- c("ID", ind_ids)
fwrite(dfOut, out_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



