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
out_file_L = args[8]

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
dfSNP <- dfSNP[,1]
colnames(dfSNP) <- "ID"
print(head(dfSNP))

# Combine SNP and R files
dfALL <- inner_join(dfSNP, df) %>% drop_na()
print(paste0("There are ", nrow(df), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)

# Standardize r values
dfALL$r <- scale(dfALL$r)

# Check dimensions of dfALL
print(paste0("The dimensions of dfAll are ", dim(dfALL)))

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
    print(head(dfR_tmp))
    tmp_r_name <- paste0(out_prefix, blockNum, ".rvec")
    fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

    # Subset single SNP ID
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

dfL <- as.data.frame(L)
fwrite(dfL, out_file_L, quote = F, row.names = F, sep = "\t")

