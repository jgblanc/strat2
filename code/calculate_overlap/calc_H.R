## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
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
print(paste0("There are ", M, " individuals in the GWAS panel"))
dfFinal <- matrix(NA, nrow = M, ncol = 1)


# Read in all values of r and all GWAS variances

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
var_file_name <- paste0(var_prefix,1, ".txt")
df <- fread(r_file_name)
df$CHR <- 1
dfVar <- fread(varFile_name)
df <- inner_join(df, dfVar)

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Read in Var file
  var_file_name <- paste0(prefix,i, ".txt")
  dfVar <- fread(var_file_name)
  tmp <- inner_join(tmp, dfVar)

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(SNP_file)

# Combine SNP and R files
df <- inner_join(dfSNP, df)
print(paste0("There are ", nrow(df), " SNPs in all the R files combined with the pruned SNPs"))

# Select 1 SNP per 1703 blocks
dfALL <- df %>% group_by(block) %>% sample_n(1)

# Standardize r values
dfALL$r <- scale(dfALL$r)

# Check dimensions of dfALL
print(paste0("The dimensions of dfAll are ", dim(dfAll)))

# Individually score each of the 1703 SNPs
for (j in 1:22) {

  # Read in R File
  dfR <- dfALL %>% filter(CHR == j)

  # Number blocks in chromosome
  nBlock_chr <- length(unique(dfR$block))
  print(paste0("There are ", nBlock_chr, " blocks on the Chr"))
  dfFGr_mat <- matrix(NA, nrow = M, ncol = nBlock_chr)

  # Divide R by sd of GWAS variance
  dfR$r <- dfR$r / sqrt(dfR$Var)

  ## Loop through blocks
  for (i in 1:nBlock_chr) {

    # Block num
    blockNum <- unique(dfR$block)[i]
    print(blockNum)

    # Subset Rs and save
    dfR_tmp <- dfR %>% filter(block == blockNum) %>% select("ID", "ALT", "r")
    tmp_r_name <- paste0(out_prefix, blockNum, ".rvec")
    fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

    # Set up plink command
    tmp_outfile <- paste0(out_prefix, blockNum)
    plink_prefix_chr <- paste0(plink_prefix, j, "_v3")
    plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr, " --keep ", id_file, " --threads 8 ",
                        " --score ", tmp_r_name, " center header-read cols=dosagesum,scoresums --out ", tmp_outfile)
    system(plink_cmd)


    ## Read in plink output
    df<- fread(paste0(tmp_outfile, ".sscore"))
    rawFGr <- as.matrix(df[,3])
    dfFGr_mat[,i] <- rawFGr

    ## Remove tmp files
    rm_cmd <- paste0("rm ", tmp_outfile, ".*")
    system(rm_cmd)

  }
  dfFinal <- cbind(dfFinal, dfFGr_mat)
}

dfFinal <- dfFinal[,2:ncol(dfFinal)]
print(paste0("The dimensions of dfFinal is ", dim(dfFinal)))

# Calculate FGr
FGr_raw <- apply(dfFinal, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Scale by 1/sqrt(L-1)
L <- nrow(dfALL)
print(paste0("The total number of SNPs is ", L))
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
M <- nrow(dfFGr)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

# Compute SE for H
nblocks <- ncol(dfFGr)
allHs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  FGri <- dfFGr[,i]
  Hi <- (sum(FGri^2)) * (1/M) * (1 / (L -1))
  allHs[i] <- (1 / (L - 1)) * (H - Hi)^2

}
varH <- mean(allHs)
se <- sqrt(varH)

pval <- pnorm(H ,mean =(1/L), sd = se, lower.tail = FALSE)
print(pval)


## Save FGr
dfOut <- as.data.frame(cbind(H[1,], pval, var(FGr), varH))
colnames(dfOut) <- c("H", "pval", "VarFGr", "varH")
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")
