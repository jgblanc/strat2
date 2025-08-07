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

# Read FID and IID to use later for merging
ids <- dfPCs[,1:2]

# Get all SNP column names
all_cols <- colnames(fread(snp_file, nrows = 0))
snp_cols <- all_cols[7:length(all_cols)]

# First write: FID, IID only
fwrite(ids, out_file, sep = "\t", quote = FALSE)

# Residualize and append each SNP
for (snp in snp_cols) {
  print(paste("Processing:", snp))

  # Read FID, IID, and the SNP column only
  y <- fread(snp_file, select = snp)[[1]]
  print(y)

  # Residualize
  model <- lm(y ~ covars)
  resids <- resid(model)

  # Write to temp file
  temp_file <- tempfile()
  fwrite(data.table(resids), temp_file, col.names = snp)

  # Paste column onto the existing output file
  system(sprintf("paste %s %s > %s", out_file, temp_file, paste0(out_file, ".tmp")))
  system(sprintf("mv %s %s", paste0(out_file, ".tmp"), out_file))

  # Clean up
  unlink(temp_file)
}




