## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

pca_file = args[1]
fgr_file = args[2]
snp_file = args[3]
out_file = args[4]


# Read PCA
PCs <- fread(pca_file)
PC_IDs <- PCs[,1:2]
PC_nums <- PCs[,3:ncol(PCs)]

# Scale the PCs to have variance 1
PC_nums <- scale(PC_nums)

# Figure out number of FGr Chrs
chrFGr <- seq(2, 22, 2)

# Read in SNP num files
dfSNP <- fread(snp_file)

# Extract only block nums with correct chr
blockIndex <- which(dfSNP$CHR %in% chrFGr)
dfSNP_filter <- dfSNP %>% filter(CHR %in% chrFGr)
L <- sum(dfSNP_filter$nSNP)
numBlocks <- length(blockIndex)

# Read in and compute FGr
dfFGr <- as.matrix(fread(fgr_file))
dfF_select <- dfFGr[,blockIndex]
FGr_raw <- apply(dfF_select, 1, sum)
FGr <- FGr_raw * (1/(sqrt(L-1)))
M <- length(FGr)

# Calculate H
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)

# Compute SE for H
allHs <- rep(NA, numBlocks)
allHi <- rep(NA, numBlocks)
for (i in 1:numBlocks) {

  mi <- as.numeric(dfSNP_filter[i,2])
  FGri <- (FGr_raw - dfF_select[,i]) * (1/sqrt(L-mi-1))
  Hi <- sum(FGri^2) * (1/M) * (1/(L-mi-1))
  allHs[i] <- ((L - mi)/mi) * (H - Hi)^2
  allHi[i] <- Hi

}
varH <- mean(allHs)

# Compute Jacknife of each FGR
jckFGr <- matrix(NA, nrow = M, ncol = numBlocks)
for (i in 1:numBlocks) {

  mi <- as.numeric(dfSNP_filter[i,2])
  FGri <- (FGr_raw - dfF_select[,i]) * (1/sqrt(L-mi-1))
  jckFGr[,i] <- (FGr - FGri)^2  * ((L - mi)/mi)
}

# Compute Numerator for error
meanJCK <- rowMeans(jckFGr)
numerator <- mean(meanJCK)
print(paste0("The numerator is ",numerator))

# Compute Denominator
varFGr <- var(FGr)

# Find Error
error <- numerator / varFGr

# Final signal
signal <- 1 - error


# Find R2 for each PC

# Compute R2 function
compute_Ratio <- function(FGr_raw, PC) {

  FGr <- FGr_raw * (1/(sqrt(L-1)))

  # Compute Jacknife of each FGR
  jckFGr <- matrix(NA, nrow = M, ncol = numBlocks)
  for (i in 1:numBlocks) {

    mi <- as.numeric(dfSNP_filter[i,2])
    FGri <- (FGr_raw - dfF_select[,i]) * (1/sqrt(L-mi-1))
    jckFGr[,i] <- (FGr - FGri)^2  * ((L - mi)/mi)
  }

  # Compute Numerator for error
  meanJCK <- rowMeans(jckFGr)
  numerator <- mean(meanJCK)

  # Compute Denominator
  varFGr <- var(FGr)

  # Find Error
  error <- numerator / varFGr

  # Final signal
  signal <- 1 - error

  # Get final Fvec
  Fvec <- scale(FGr)
  mod <- lm(Fvec ~ PC)
  R2 <- summary(mod)$r.squared
  Ratio <- R2/signal

  return(Ratio)
}

# Bootstrap function
bootstrap_ratio_ci <- function(FGr_raw, PC, n_boot = 1000, conf = 0.95) {

  PC <- as.matrix(PC)
  n <- nrow(FGr_raw)

  boot_ratios <- replicate(n_boot, {
    idx <- sample(n, replace = TRUE)
    compute_Ratio(FGr_raw[idx, , drop = FALSE], PC[idx, , drop = FALSE])
  })

  alpha <- 1 - conf
  ci <- quantile(boot_ratios, probs = c(alpha / 2, 1 - alpha / 2))
  return(list(
    estimate = compute_Ratio(Fvec, PC, signal),
    se = sd(boot_ratios),
    ci = ci
  ))
}

# Construct output
dfOut <- matrix(NA, nrow=ncol(PC_nums), ncol = 10)
colnames(dfOut) <- c("H","varH", "Signal","PC", "B2", "R2", "Ratio", "lc", "uc", "se")
FGr_scale <- scale(FGr)


# Loop through PCs
for (i in 1:ncol(PC_nums)) {

  print(i)
  # Calc cov
  B <- cov(FGr_scale, PC_nums[,i])
  B2 <- B^2

  # Collect output
  dfOut[i, 1] <- H
  dfOut[i, 2] <- varH
  dfOut[i,3] <- signal
  dfOut[i, 4] <- i
  dfOut[i,5] <- B2
  R2 <- sum(dfOut[1:i, 5])
  dfOut[i,6] <- R2
  Ratio <- R2/signal
  dfOut[i,7] <- Ratio

  # Get CI
  result <- bootstrap_ratio_ci(FGr_raw, PC_nums[,1:i], n_boot = 1000, conf = 0.95)
  dfOut[i,10] <- result$se
  dfOut[i,8] <- as.numeric(result$ci[1])
  dfOut[i,9] <- as.numeric(result$ci[2])

}

# Save output file
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



