## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)
library(matrixStats)
library(future.apply)  # parallel bootstrap
plan(multisession)
plan(multicore, workers = 8)
options(future.globals.maxSize = 2 * 1024^3)

pca_file = args[1]
fgr_file = args[2]
snp_file = args[3]
out_file = args[4]

# Read inputs
PCs <- fread(pca_file)
PC_nums <- scale(as.matrix(PCs[, 3:ncol(PCs)]))
dfSNP <- fread(snp_file)

# Get even chr
chrFGr <- seq(2, 22, 2)
dfSNP_filter <- dfSNP %>% filter(CHR %in% chrFGr)
blockIndex <- which(dfSNP$CHR %in% chrFGr)
L <- sum(dfSNP_filter$nSNP)
numBlocks <- length(blockIndex)

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


# Bootstrap helper
compute_Ratio <- function(Fmat, PC) {

  FGr_raw <- rowSums(Fmat)
  FGr <- FGr_raw / sqrt(L - 1)
  mod <- lm(scale(FGr) ~ PC)
  R2 <- summary(mod)$r.squared

  # Jackknife error for this Fmat
  FGri_mat <- FGr_raw - Fmat
  FGri_scaled <- sweep(FGri_mat, 2, scale_factors, `*`)
  jckFGr <- ((FGr - FGri_scaled)^2) * matrix((L - mi_vec) / mi_vec, nrow(Fmat), numBlocks, byrow = TRUE)
  numerator <- mean(rowMeans(jckFGr))
  signal <- 1 - (numerator / var(FGr))

  # Ratio
  ratio <- R2/signal

  return(list(R2 = R2, Ratio = ratio))
}

# Construct output
dfOut <- matrix(NA, nrow = ncol(PC_nums), ncol = 14)
colnames(dfOut) <- c("PC","H", "varH","signal", "B", "lcB", "ucB", "B2", "R2", "lcR2", "ucR2", "Ratio", "lcRatio", "ucRatio")

# Loop through PCs
for (i in seq_len(ncol(PC_nums))) {

  # Get Single PC stats
  mod <- lm(scale(FGr)~ PC_nums[,i])
  B <- mod$coefficients[2]
  lcB <- confint(mod)[-1, ][1]
  ucB <- confint(mod)[-1, ][2]
  B2 <- B^2

  # Get Combined R^2 stats
  tmp <- compute_Ratio(Fmat = dfFGr, PC = PC_nums[,1:i])
  R2 <- tmp$R2

  # Run bootstrappin
  #boot_ratios <- future_replicate(100, {
  #  idx <- sample(M, replace = TRUE)
  #  compute_Ratio(dfFGr[idx,], tmpPC <- PC_nums[idx,1:i])
  #})
  #boot_ratios_matrix <- matrix(unlist(boot_ratios), ncol = 2, byrow = TRUE)

  # Get CI's for R2
  #boot_R2_values <- boot_ratios_matrix[,1]
  #ci <- quantile(boot_R2_values, probs = c(0.025, 0.975))
  #lcR2 <- ci[1]
  #ucR2 <- ci[2]

  # Get Ratio
  Ratio <- tmp$Ratio

  # Get CI's for Ratio
  #boot_ratio_values <- boot_ratios_matrix[,2]
  #ci <- quantile(boot_ratio_values, probs = c(0.025, 0.975))
  #lcRatio <- ci[1]
  #ucRatio <- ci[2]

  dfOut[i,] <- c(i, H, varH, signal, B, lcB, ucB, B2, R2, NA, NA, Ratio, NA, NA)
  print(dfOut)
  cat("Finished PC", i, "\n")
}


# Write output
fwrite(as.data.table(dfOut), out_file, sep = "\t", quote = FALSE)






