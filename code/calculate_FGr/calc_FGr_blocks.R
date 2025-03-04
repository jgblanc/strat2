## Compute FGr for a single Chromosome with all SNPs

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr_single_chr.R <prefix to plink files> <freq file> <r> <FGR>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

plink_prefix = args[1]
out_prefix = args[2]
FGr_outfile = args[3]
SNP_outfile = args[4]
var_prefix = args[5]
r_file = args[6]
id_file = args[7]

# Read in R
dfR_ALl <- fread(r_file)
colnames(dfR_ALl)[1] <- "CHR"

# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)
print(M)
dfFinal <- matrix(NA, nrow = M, ncol = 1)
SNPcounter <- rep(0, 581)

# Read in all GWAS variance files
dfVar <- fread(paste0(var_prefix, "1.txt"))
for (j in 2:22) {
  tmp <- fread(paste0(var_prefix, j, ".txt"))
  dfVar <- rbind(dfVar, tmp)
}
print(paste0("There are ", nrow(dfVar), " variants in the GWAS var file"))

# Combine R and Var file
dfR_ALL <- inner_join(dfR_ALL, dfVar)
print(head(dfR_ALL))

total_blocks <- unique(dfR_ALL$block)

## Loop through chromosomes
for (j in 1:22) {

  ## Read in R File
  dfR <- dfR_ALL %>% filter(CHR == j)

  ## Number blocks in chromosome
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
    SNPcounter[blockNum] <- nrow(dfR_tmp)
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
    print(var(rawFGr * (1/sqrt(nrow(dfR_tmp)))))
    print(var(rawFGr * sqrt(100000/1286958) * (1/sqrt(nrow(dfR_tmp)))))
    dfFGr_mat[,i] <- rawFGr

    ## Remove tmp files
    rm_cmd <- paste0("rm ", tmp_outfile, ".*")
    system(rm_cmd)

  }
  dfFinal <- cbind(dfFinal, dfFGr_mat)
}

dfFinal <- dfFinal[,2:ncol(dfFinal)]
print(paste0("The dimensions of dfFinal is ", dim(dfFinal)))
print(SNPcounter)


# Save FGr of each block
dfOut <- as.data.frame(cbind(df[,1], dfFinal))
fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")

# Save SNP counts
dfSNP <- as.data.frame(cbind(seq(1,581), SNPcounter))
fwrite(dfSNP, SNP_outfile, quote = F, row.names = F, sep = "\t")

# Calculate FGr
#FGr_raw <- apply(dfFinal, 1, sum)
#print(paste0("The raw var is ", var(FGr_raw)))

## Scale by 1/sqrt(L-1)
#L <- sum(SNPcounter)
#print(L)
#FGr <- FGr_raw * (1/(sqrt(L-1)))
#print(paste0("The scaled var is ", var(FGr)))

# Calculate H
#M <- nrow(df)
#H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
#print(paste0("H is ", H))
#print(paste0("1/L is ", 1/L))

## Save FGr
#dfOut <- as.data.frame(cbind(df[,1], FGr))
#fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")
