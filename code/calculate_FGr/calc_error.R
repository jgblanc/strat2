## Compute FGr using one PC SNPs

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

r_prefix = args[1]
Gr_mat_file = args[2]
snp_num_file = args[3]
out_file = args[4]
pc_snp_file = args[5]


####################################
########## Functions ###############
####################################


calc_fhat <- function(dfMat, r ) {

  fhat_raw <- apply(dfMat, 1, sum)
  rTr <- as.numeric(t(as.matrix(r)) %*% as.matrix(r))
  fhat <- fhat_raw / c(rTr)

  return(fhat)
}


####################################
########## Main  ###################
####################################


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
dfSNP <- fread(pc_snp_file) %>% select("ID", "block")
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Read in Fmat
dfMat <- as.matrix(fread(Gr_mat_file))
M <- nrow(dfMat)

# Read in block info
dfSNP_Num <- fread(snp_num_file)
numBlocks <- as.numeric(ncol(dfMat))


# Compute "real" fhat
fhat <- calc_fhat(dfMat, dfALL$r)

# Leave one out f;s
locoFGr <- matrix(NA, nrow = M, ncol = numBlocks)
for (i in 1:numBlocks) {

  # Block num
  blockNum <- as.numeric(dfSNP_Num[i,1])

  print(paste0("This is rep number ",i))
  dfR_not_i <- dfALL %>% filter(block != blockNum)
  fhat_i <- calc_fhat(dfMat[,-i],dfR_not_i$r)
  locoFGr[,i] <- fhat_i
}



# Compute numerator
mean_loco <- rowMeans(locoFGr)
jckFGr <- matrix(NA, nrow = M, ncol = numBlocks)
for (i in 1:numBlocks) {

  print(paste0("This is rep number ",i))
  jckFGr[,i] <- (locoFGr[,i] - mean_loco)^2

}
tmp <- rowSums(jckFGr) * ((numBlocks - 1)/numBlocks)
numerator <- mean(tmp)
print(paste0("The numerator is ",numerator))


# Compute Denominator
varFGr <- var(fhat)
print(paste0("The denominator is ",varFGr))

# Find Error
error <- numerator / varFGr

# Final signal
signal <- 1 - error


# Construct output
dfOut <- data.frame(error = error, signal = signal)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")






