## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
out_prefix = args[3]
pca_prefix = args[4]
snp_file = args[5]
id_file = args[6]
out_file = args[7]


# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
df <- fread(r_file_name)
df$CHR <- 1

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID", "block")


# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
L <- nrow(dfALL)
print(paste0("There are ", L, " snps total"))


# Compute Gr

# Extrac SNPs
dfSNP_tmp <- dfALL %>% select("ID")
snp_name <- paste0(out_prefix, ".snp")
fwrite(dfSNP_tmp, snp_name, quote = F, row.names = F, sep = "\t")

# Write r file
tmp_r_name <- paste0(out_prefix, ".rvec")
dfR_tmp <- dfALL %>% select("ID", "ALT", "r")
fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

# Run plink command
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --extract ", snp_name ," --threads 8 --score ", tmp_r_name, " header-read variance-standardize cols=dosagesum,scoresums --out ", out_prefix)
system(plink_cmd)

# Read in GR
df <- fread(paste0(out_prefix, ".sscore"))
fhat_raw <- as.matrix(df$SCORE1_SUM)

# calculate r^\top t
r_center <- dfALL$r - mean(dfALL$r)
rTr <- t(as.matrix(r_center)) %*% as.matrix(r_center)

# calculate \hat{f}
fhat<- fhat_raw / c(rTr)

# calculate sigma2f
fhat_center <- fhat - mean(fhat)
numerator <- as.numeric(t(fhat_center) %*% fhat_center)
sigma2f <- numerator / (M-1)
print(paste0("Sigma2f is ", sigma2f))
print(paste0("the var is  is ", var(fhat_center)))


# Regress out PCs from fhat
pca_file_path_common <- paste0(pc_prefix, "even_PCA.eigenvec")
pca_file_path_rare   <- paste0(pc_prefix, "rare_even_PCA.eigenvec")

dfPCs_common <- fread(pca_file_path_common)
setnames(dfPCs_common, 1:2, c("FID","IID"))

dfPCs_rare <- fread(pca_file_path_rare)
setnames(dfPCs_rare, 1:2, c("FID","IID"))

dfPCs <- inner_join(dfPCs_common, dfPCs_rare, by = c("FID","IID"))
covars_df <- dfPCs %>% select(starts_with("PC"))

y <- as.numeric(fhat_center)
fhat_resid <- resid(lm(y ~ ., data = covars_df))

# calculate sigma2f prime
fhat_resid_center <- fhat_resid - mean(fhat_resid)
numerator <- as.numeric(t(fhat_resid_center) %*% fhat_resid_center)
sigma2f_prime <- numerator / (M-1)
print(paste0("Sigma2f is ", sigma2f_prime))
print(paste0("the var is  is ", var(fhat_resid_center)))

# Residualize GR and then mulitple by rTr
y <- as.numeric(fhat_raw)
gr_resid <- resid(lm(y ~ ., data = covars_df))
fhat_resid<- gr_resid / c(rTr)

# calculate sigma2f prime
fhat_resid_center <- fhat_resid - mean(fhat_resid)
numerator <- as.numeric(t(fhat_resid_center) %*% fhat_resid_center)
sigma2f_prime_post <- numerator / (M-1)
print(paste0("Sigma2f is ", sigma2f_prime_post))
print(paste0("the var is  is ", var(fhat_resid_center)))




# Save SNP file
dfOut <- data.frame(sigma2f, sigma2f_prime, sigma2f_prime_post)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



