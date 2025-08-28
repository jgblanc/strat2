## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

r_prefix = args[1]
FGr_prefix = args[2]
snpnum_file = args[3]
tvec_file = args[4]
out_file = args[5]
snp_file = args[6]

####################################
########## Functions ###############
####################################

calc_sigma2_f <- function(fhat, M) {

  numerator <- t(fhat) %*% fhat
  out <- numerator / (M-1)

  return(out)
}


calc_sigma2_r <-function(rmat, L) {

  rTr <- t(as.matrix(r)) %*% as.matrix(r)
  out <- rTr / (L - 1)

  return(out)
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
dfSNP <- fread(snp_file) %>% select("ID", "block")
dfSNP <- inner_join(dfPvar, dfSNP) %>% select("ID", "block")
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)


# Read in Fmat
dfMat <- as.matrix(fread(FGr_file))
M <- nrow(dfMat)

### Calculate H Real
fhat <- apply(fMat, 1, sum)
sigma2F <- calc_sigma2_f(dfMat, M)
sigma2r <- calc_sigma2_r(dfALL$r, L)
H <- sigma2F * sigma2r
print(paste0("H is ", H))

# Read in block info
dfSNPs <- fread(snpnum_file)
num_blocks <- ncol(dfMat)

#### Calculate SE via block jackknife
allHs <- rep(NA, numBlocks)
for (i in 1:numBlocks) {

  # Block num
  blockNum <- dfSNPs[i,1]

  # SNP num
  mi <- dfSNPs[i,2]
  print(paste0("The SNP number is ", mi))

  # Calc H
  fhat_i <- fhat - dfFGr_mat[,i]
  sigma2F_i <- calc_sigma2_f(fhat_i, M)
  dfR_i <- dfALL %>% filter(block == blockNum)
  print(paste0("The number of rows in dfR_i is ", nrow(dfR_i)))
  sigma2r_i <- calc_sigma2_r(dfR_i$r,mi)
  Hi <- sigma2F_i * sigma2r_i
  allHs[i] <- ((L - mi)/mi) * (H - Hi)^2


}

# Calculate SE
varH <- mean(allHs)

# P-value from jacknife
pvalNorm <- pnorm(H, mean = 1/L, sd=sqrt(varH), lower.tail = FALSE)

# Get N
dfTvec <- fread(tvec_file)
N <- nrow(dfTvec)

# Construct output
dfOut <- data.frame(H = H, L = L, varH = varH, pvalNorm = pvalNorm, sigma2r=sigma2r, sigma2F=sigma2F, N=N)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")




