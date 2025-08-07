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


# Read FID and IID to use later for merging
ids <- fread(snp_file, select = c("FID", "IID"))

# Merge with PCs once
df_base <- inner_join(ids, dfPCs, by = c("FID", "IID"))
print(head(df_base))

# Get all SNP column names
all_cols <- colnames(fread(snp_file, nrows = 0))
snp_cols <- all_cols[7:length(all_cols)]
head(snp_cols)

snp_cols <- snp_cols[1:10]

# Create output file with header
fwrite(data.table(FID=character(), IID=character(), SNP=character(), Residual=numeric()),
       out_file, sep="\t")

# Residualize and append each SNP
for (snp in snp_cols) {
  print(paste("Processing:", snp))

  # Read FID, IID, and the SNP column only
  df_snp <- fread(snp_file, select = c("FID", "IID", snp))

  # Merge with PCs
  df <- inner_join(df_snp, dfPCs, by = c("FID", "IID"))

  # Extract genotype vector
  y <- df[[snp]]

  # Extract PCs as covariates
  covars <- as.matrix(df %>% select(starts_with("PC")))

  # Residualize
  model <- lm(y ~ covars)
  resids <- resid(model)

  # Save output
  result <- data.table(FID = df$FID, IID = df$IID, SNP = snp, Residual = resids)
  fwrite(result, out_file, sep="\t", append=TRUE)
}




