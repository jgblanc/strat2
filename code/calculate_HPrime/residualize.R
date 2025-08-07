## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

snp_file = args[1]
pc_file = args[2]
out_file = args[3]


# Read in dosages
dfSNPs <- fread(id_file)
print(head(dfSNPs))

# Read in PCs
dfPCs <- nrow(dfIDs)
print(head(dfPCs))

# Combine
dfALL <- inner_join(dfSNPs,dfPCs)

# Extract PCs
pc_matrix <- as.matrix(dfALL[ , grep("^PC", names(dfALL)) ])
print(dim(pc_matrix))

# Extract genotypes
geno_matrix <- as.matrix(dfALL[ , grepl("^:", names(dfALL)) ])
print(dim(geno_matrix))

# Function to residualize a vector on covariates
residualize <- function(y, covars) {
  print("LM!")
  lm.fit <- lm(y ~ covars)
  return(resid(lm.fit))
}

# Apply across all SNPs
geno_resid <- apply(geno_matrix, 2, residualize, covars = pc_matrix)

# Add back FID and IID
geno_resid <- cbind(FID = dfALL$FID, IID = dfALL$IID, geno_resid)

# Construct output
write.table(geno_resid, out_file, sep = "\t", quote = FALSE, row.names = FALSE)



