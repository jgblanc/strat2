## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_R2_EO.R")}


library(data.table)
library(tidyverse)
library(matrixStats)
library(future.apply)  # parallel bootstrap
plan(multisession)
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

# Bootstrap helper
compute_Ratio <- function(Fmat, PC) {

  PC <- scale(PC)
  FGr_raw <- rowSums(Fmat)
  FGr <- FGr_raw / sqrt(L - 1)
  Fvec <- scale(FGr)
  mod <- lm(Fvec ~ PC)
  R2 <- summary(mod)$r.squared

  # Jackknife error for this Fmat
  FGri_mat <- FGr_raw - Fmat
  FGri_scaled <- sweep(FGri_mat, 2, scale_factors, `*`)
  jckFGr <- ((FGr - FGri_scaled)^2) * matrix((L - mi_vec) / mi_vec, nrow(Fmat), numBlocks, byrow = TRUE)
  numerator <- mean(rowMeans(jckFGr))
  signal <- 1 - numerator / var(FGr)

  R2 / signal
}

bootstrap_ratio_ci <- function(Fmat, PC, n_boot = 1000, conf = 0.95) {
  PC <- as.matrix(PC)
  n <- nrow(Fmat)
  boot_ratios <- future_replicate(n_boot, {
    idx <- sample(n, replace = TRUE)
    compute_Ratio(Fmat[idx, , drop = FALSE], PC[idx, , drop = FALSE])
  })

  se <- sd(boot_ratios)
  estimate <- compute_Ratio(Fmat, PC)
  z <- qnorm(1 - (1 - conf)/2)
  ci <- estimate + c(-1, 1) * z * se

  list(
    estimate = estimate,
    se = se,
    ci = ci
  )
}

# Construct output
dfOut <- matrix(NA, nrow = ncol(PC_nums), ncol = 11)
colnames(dfOut) <- c("H","varH", "Signal","PC", "B2", "R2", "Ratio", "lc", "uc", "se", "estimate")
FGr_scale <- scale(FGr)

# Loop through PCs
R2_cum <- 0
for (i in seq_len(ncol(PC_nums))) {

  B2 <- cov(FGr_scale, PC_nums[,i])^2
  R2_cum <- R2_cum + B2
  Ratio <- R2_cum / signal

  ci_result <- bootstrap_ratio_ci(dfFGr, PC_nums[,1:i], n_boot = 1000, conf = 0.95)

  dfOut[i,] <- c(H, varH, signal, i, B2, R2_cum, Ratio, ci_result$ci[1], ci_result$ci[2], ci_result$se, ci_result$estimate)
  cat("Finished PC", i, "\n")
}

# Write output
fwrite(as.data.table(dfOut), out_file, sep = "\t", quote = FALSE)






