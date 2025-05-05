## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)
library(matrixStats)
library(future.apply)  # parallel bootstrap
plan(multisession)
options(future.globals.maxSize = 2 * 1024^3)

common_pca_file = args[1]
rare_pca_file = args[2]
fgr_file = args[3]
snp_file = args[4]
out_file = args[5]

# Read inputs
cPCs <- fread(common_pca_file)
cPC_nums <- scale(as.matrix(cPCs[, 3:ncol(cPCs)]))
rPCs <- fread(rare_pca_file)
rPC_nums <- scale(as.matrix(rPCs[, 3:ncol(rPCs)]))
dfSNP <- fread(snp_file)

# Get even chr
chrFGr <- seq(2, 22, 2)
dfSNP_filter <- dfSNP %>% filter(CHR %in% chrFGr)
blockIndex <- which(dfSNP$CHR %in% chrFGr)
L <- sum(dfSNP_filter$nSNP)
numBlocks <- length(blockIndex)
print(numBlocks)

# Read in and compute FGr
dfFGr <- as.matrix(fread(fgr_file))[, blockIndex]
FGr_raw <- rowSums(dfFGr)
FGr <- FGr_raw / sqrt(L - 1)
M <- length(FGr)

# Calculate H
H <- (1 / (M * (L - 1))) * sum(FGr^2)
print(H)

# Compute jackknife estimates
mi_vec <- dfSNP_filter$nSNP
FGri_mat <- FGr_raw - dfFGr
scale_factors <- 1 / sqrt(L - mi_vec - 1)
FGri_scaled <- sweep(FGri_mat, 2, scale_factors, `*`)
Hi_vec <- colSums(FGri_scaled^2) / (L - mi_vec - 1) / M
allHs <- ((L - mi_vec) / mi_vec) * (H - Hi_vec)^2
varH <- mean(allHs)

# Jackknife variance of each FGr element
FGr_mat <- matrix(rep(FGr, numBlocks), nrow = M)
FGri_mat2 <- FGri_scaled
jckFGr <- ((FGr_mat - FGri_mat2)^2) * matrix((L - mi_vec) / mi_vec, M, numBlocks, byrow = TRUE)
meanJCK <- rowMeans(jckFGr)
numerator <- mean(meanJCK)
varFGr <- var(FGr)
error <- numerator / varFGr
signal <- 1 - error


# Construct output
PC_nums <- cbind(cPC_nums, rPC_nums)
print(dim(PC_nums))
dfOut <- matrix(NA, nrow = ncol(PC_nums), ncol = 11)
colnames(dfOut) <- c("H","varH", "Signal","PC", "B2", "R2", "Ratio", "lc", "uc", "se", "estimate")
FGr_scale <- scale(FGr)

# Loop through PCs
for (i in 1:ncol(PC_nums)) {

  B2 <- cov(FGr_scale, PC_nums[,i])^2
  mod <- lm(FGr_scale ~ PC_nums[,1:i])
  R2 <- summary(mod)$r.squared
  Ratio <- R2 / signal

  #ci_result <- bootstrap_ratio_ci(dfFGr, PC_nums[,1:i], n_boot = 1000, conf = 0.95)

  dfOut[i,] <- c(H, varH, signal, i, B2, R2, Ratio, NA, NA, NA, NA)
  cat("Finished PC", i, "\n")
}

# Write output
fwrite(as.data.table(dfOut), out_file, sep = "\t", quote = FALSE)






