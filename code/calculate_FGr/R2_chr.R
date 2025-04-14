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

# Get the PCA Chrs
tmp <- strsplit(out_file, "_")[[1]][4]
chrPCA <- as.numeric(strsplit(tmp, ".txt")[[1]][1])

# Figure out number of FGr Chrs
chrFGr <-seq(1, 22)[ seq(1, 22) >  chrPCA]

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


# Compute Jacknife of each FGR
jckFGr <- matrix(NA, nrow = M, ncol = numBlocks)
for (i in 1:numBlocks) {

  print(paste0("This is rep number ",i))
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

# Construct output
dfOut <- matrix(NA, nrow=ncol(PC_nums), ncol = 5)
colnames(dfOut) <- c("PC", "B2", "R2", "Signal", "Ratio")
FGr_scale <- scale(FGr)

# Loop through PCs
for (i in 1:ncol(PC_nums)) {

  # Calc cov
  B <- cov(FGr_scale, PC_nums[,i])
  B2 <- B^2

  # Collect output
  dfOut[i, 1] <- i
  dfOut[i,2] <- B2
  R2 <- sum(dfOut[1:i, 2])
  dfOut[i,3] <- R2
  dfOut[i,4] <- signal
  dfOut[i,5] <- R2/signal

}

# Save output file
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



