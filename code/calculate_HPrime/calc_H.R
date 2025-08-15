## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

resid_prefix = args[1]
r_prefix = args[2]
snp_file = args[3]
id_file = args[4]
out_file = args[5]


# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r

## First Chr
r_file_name <- paste0(r_prefix, 22, ".rvec")
df <- fread(r_file_name)
df$CHR <- 22

## All other Chrs
#for (i in 2:22) {

  # Read in R file
#  r_file_name <- paste0(r_prefix, i, ".rvec")
#  tmp <- fread(r_file_name)
#  tmp$CHR <- i

  # Combine
#  df <- rbind(df, tmp)
#}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID", "block")
dfSNP <- dfSNP %>%
  separate(id, into = c("CHR", "POS"), sep = ":", remove = FALSE)
print(paste0("Number of PC SNPs ", nrow(dfSNP)))

# Combine SNP and R files
dfALL <- inner_join(df, dfSNP) %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)
print(L)

# Standardize r values variance
dfALL$r <- dfALL$r / sd(dfALL$r)
print(paste0("The variance of r is ", var(dfALL$r)))
print(paste0("The mean of r is ", mean(dfALL$r)))

# Set up data frame to collect SNP numbers
numBlocks <- length(unique(dfALL$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 2, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP")
print(paste0("The total number of blocks is ", numBlocks))
dfFGr_mat <- matrix(NA, nrow = M, ncol = numBlocks)




for (i in 22:22) {

  # Get r values on Chr i
  dfR_chr <- dfALL %>% filter(CHR == i)
  chr_blocks <- unique(dfR_chr$block)
  num_blocks_chr <- length(chr_blocks)
  print(paste0("The number of blocks on chr ", i, " is ", num_blocks_chr))

  # Get SNP ids on Chr i
  traw_file <- paste0(resid_prefix, i, ".traw")
  cmd_file<- paste0( "cut -f1 ",traw_file)
  dfSNPchr <- fread(cmd_file,header = TRUE)

  for (b in chr_blocks) {

    print(paste0("We are on block", b))
    # Select r values on block b
    dfR_block <- dfR_chr %>% filter(block == b)
    r_ids <- dfR_block$ID

    # Get the right row values
    row_ids <- which(dfSNPchr$SNP %in% r_ids)
    print(paste0("There are ", length(row_ids), " SNPs on the block" ))
    print(head(row_ids))

    # read only correct row IDs
    row_nums_str <- paste(row_ids + 1, collapse = ",")  # +1 for header row
    cmd_block <- paste0("awk 'NR==1 || NR==", gsub(",", " || NR==", row_nums_str), "' ", traw_file)
    dfBlock <- fread(cmd_block)

    # do matrix multiplication
    matBlock <- t(as.matrix(dfBlock))
    print(paste0("The dim of matBlock is ", dim(matBlock)))

    rawFGr <- matBlock %*% dfR_block$r
    dfFGr_mat[,b] <- rawFGr

    # Save number of SNPs
    nsnp_in_block <- nrow(dfR_block)
    dfSNPs[b,1] <- blockNum
    dfSNPs[b,2] <- nrow(dfR_block)

  }
}

# Calculate FGr
FGr_raw <- apply(dfFGr_mat, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Scale by 1/sqrt(L-1)
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

# Compute SE for H
allHs <- rep(NA, numBlocks)
allHi <- rep(NA, numBlocks)
for (i in 1:numBlocks) {

  print(paste0("This is rep number ",i))
  mi <- as.numeric(dfSNPs[i,2])
  FGri <- (FGr_raw - dfFGr_mat[,i]) * (1/sqrt(L-mi-1))
  Hi <- sum(FGri^2) * (1/M) * (1/(L-mi-1))
  allHs[i] <- ((L - mi)/mi) * (H - Hi)^2
  allHi[i] <- Hi

}
varH <- mean(allHs)
meanH <- mean(allHi)

# P-value from sims
pvalNorm <- pnorm(H, mean = 1/L, sd=sqrt(varH), lower.tail = FALSE)

# Construct output
dfOut <- data.frame(H = H, L = L, meanH = meanH, varH = varH, pvalNorm = pvalNorm)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



