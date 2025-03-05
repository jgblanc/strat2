## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
var_prefix = args[3]
out_prefix = args[4]
snp_file = args[5]
id_file = args[6]
out_file = args[7]


# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r and all GWAS variances

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
var_file_name <- paste0(var_prefix,1, ".txt")
df <- fread(r_file_name)
df$CHR <- 1
dfVar <- fread(var_file_name)
df <- inner_join(df, dfVar)

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Read in Var file
  var_file_name <- paste0(var_prefix,i, ".txt")
  dfVar <- fread(var_file_name)
  tmp <- inner_join(tmp, dfVar)

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))
print(head(df))

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID")

# Combine SNP and R files
dfALL <- inner_join(dfSNP, df, by = "ID") %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)

# Divide SNPs into 1000 equal blocks
block_size <- floor(nrow(dfALL) / num_blocks)
print(block_size)

# Standardize r values and divide by GWAS variance
dfALL$r <- scale(dfALL$r)
dfALL$r <- dfR$r / sqrt(dfR$Var)

# Make a collector for all values of H
allH <- rep(NA, 21)
dfR <- dfALL

for (i in 1:length(allH)) {

  # Subset Rs and save
  dfR_tmp <- dfR %>% select("ID", "ALT", "r")
  tmp_r_name <- paste0(out_prefix, ".rvec")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Subset SNP IDs
  dfSNP_tmp <- dfR_tmp %>% select("ID")
  tmp_snp_name <- paste0(out_prefix, ".snp")
  fwrite(dfSNP_tmp, tmp_snp_name, quote = F, row.names = F, sep = "\t")

  # Set up plink command
  tmp_outfile <- out_prefix
  plink_prefix_chr <- paste0(plink_prefix)
  plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr, " --keep ", id_file, " --extract ", tmp_snp_name ," --threads 8 ",
                      " --score ", tmp_r_name, " center header-read cols=dosagesum,scoresums --out ", tmp_outfile)
  system(plink_cmd)

  # Read in plink output
  df<- fread(paste0(tmp_outfile, ".sscore"))
  rawFGr <- as.matrix(df[,3])

  # Calculate FGr
  print(paste0("The raw var is ", var(rawFGr)))
  FGr <- FGr_raw * (1/(sqrt(L-1)))
  print(paste0("The scaled var is ", var(FGr)))

  # Calculate H
  H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
  print(paste0("H is ", H))
  allH[i] <- H

  # Shift the dfR dataframe
  dfR <- dfR %>% mutate(r = c(tail(r, i * blockSize), head(r, -i * blockSize)))

}

# Calculate p-value
realH <- allH[1]
varH <- var(allH[2:length(allH)])
se = sqrt(varH)
meanH <- var(allH[2:length(allH)])
pvalNorm <- pnorm(H ,mean =meanH, sd = se, lower.tail = FALSE)
pvalSim <- sum(realH >= allH[2:length(allH)]) / (length(allH) - 1)

dfOut <- as.data.frame(c(realH, L, meanH, varH, pvalNorm, pvalSim))
colnames(dfOut) <- c("H", "L", "meanH", "varH", "pvalNorm", "pvalSim")
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")


