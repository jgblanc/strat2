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

# Read FID and IID to use later for merging
dfSNPs <- fread(snp_file)
print(dim(dfSNPs))

# Get all SNP column names
all_cols <- colnames(fread(snp_file, nrows = 0))
snp_cols <- all_cols[7:length(all_cols)]
L <- length(snp_cols)
M <- nrow(dfPCs)

# First write: FID, IID only
fwrite(ids, out_file, sep = "\t", quote = FALSE)

# Residualize and append each SNP
for (i in 1:L) {
  print(paste("Processing:", i))

  # Read FID, IID, and the SNP column only
  y <- scan(snp_file, what = numeric(), skip = 1, sep = "\t")[seq(from = i+7, to = L * M, by = L+7)]
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




