## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
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

# Read in all values of r and all GWAS variances

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
var_file_name <- paste0(var_prefix,1, ".txt")
df <- fread(r_file_name)
df$CHR <- 1
dfVar <- fread(var_file_name)
df <- inner_join(df, dfVar)

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Read in Var file
  var_file_name <- paste0(var_prefix,i, ".txt")
  dfVar <- fread(var_file_name)
  tmp <- inner_join(tmp, dfVar)

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))
print(head(df))

# Read in SNP file
dfSNP <- fread(snp_file)

# Combine SNP and R files
dfALL <- inner_join(dfSNP, df, by = "ID") %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)

# Standardize r values
dfALL$r <- scale(dfALL$r)

print(head(dfALL))
# Calculate FGr
dfFGr_mat <- matrix(NA, nrow = M, ncol = 22)
for (j in 1:22) {

  # Read in R File
  dfR <- dfALL %>% filter(CHR == j)

  # Divide R by sd of GWAS variance
  dfR$r <- dfR$r / sqrt(dfR$Var)

  # Subset Rs and save
  dfR_tmp <- dfR %>% select("ID", "ALT", "r")
  tmp_r_name <- paste0(out_prefix, ".rvec")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Subset SNP IDs
  dfSNP_tmp <- dfR_tmp %>% select("ID")
  tmp_snp_name <- paste0(out_prefix, ".snp")
  fwrite(dfSNP_tmp, tmp_snp_name, quote = F, row.names = F, sep = "\t")

  # Set up plink command
  tmp_outfile <- out_prefix
  plink_prefix_chr <- paste0(plink_prefix, j, "_v3")
  plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr, " --keep ", id_file, " --extract ", tmp_snp_name ," --threads 8 ",
                      " --score ", tmp_r_name, " center header-read cols=dosagesum,scoresums --out ", tmp_outfile)
  system(plink_cmd)

  # Read in plink output
  df<- fread(paste0(tmp_outfile, ".sscore"))
  rawFGr <- as.matrix(df[,3])
  dfFGr_mat[,j] <- rawFGr

  # Remove tmp files
  rm_cmd <- paste0("rm ", tmp_outfile, ".*")
  system(rm_cmd)
}

# Calculate FGr
FGr_raw <- apply(dfFGr_mat, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))





















# Set up data frame to collect SNP numbers
numBlocks <- length(unique(dfALL$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 2, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP")
print(paste0("The total number of blocks is ", numBlocks))


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
    print(nrow(dfR_tmp))
    tmp_r_name <- paste0(out_prefix, blockNum, ".rvec")
    fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

    # Save number of SNPs
    nsnp_in_block <- nrow(dfR_tmp)
    dfSNPs[blockNum,1] <- blockNum
    dfSNPs[blockNum,2] <- nrow(dfR_tmp)
    print(head(dfSNPs))

    # Subset SNP IDs
    dfSNP_tmp <- dfR_tmp %>% select("ID")
    tmp_snp_name <- paste0(out_prefix, blockNum, ".snp")
    fwrite(dfSNP_tmp, tmp_snp_name, quote = F, row.names = F, sep = "\t")

    # Set up plink command
    tmp_outfile <- paste0(out_prefix, blockNum)
    plink_prefix_chr <- paste0(plink_prefix, j, "_v3")
    plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr, " --keep ", id_file, " --extract ", tmp_snp_name ," --threads 8 ",
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

dfOut <- as.data.frame(dfFinal)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")

fwrite(dfSNPs, out_file_L, quote = F, row.names = F, sep = "\t")

