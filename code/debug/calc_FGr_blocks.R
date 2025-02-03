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


## Collect data
dfFinal <- matrix(NA, nrow = 9999, ncol = 1)
SNPcounter <- rep(0, 22)

## Loop through chromosomes
for (j in 1:22) {

  ## Read in R File
  r_file <- paste0(out_prefix, "r", j, "_standardize_blocks.rvec")
  dfR <- fread(r_file)
  SNPcounter[j] <- nrow(dfR)

  ## Number blocks in chromosome
  nBlock_chr <- length(unique(dfR$block))
  print(paste0("There are ", nBlock_chr, " blocks on the Chr"))
  dfFGr_mat <- matrix(NA, nrow = 9999, ncol = nBlock_chr)

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
    plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr,
                        " --score ", tmp_r_name, " variance-standardize header-read cols=dosagesum,scoresums --out ", tmp_outfile)
    system(plink_cmd)

    ## Read in plink output
    df<- fread(paste0(tmp_outfile, ".sscore"))
    rawFGr <- as.matrix(df[,3])
    dfFGr_mat[,i] <- rawFGr

  }
  dfFinal <- cbind(dfFinal, dfFGr_mat)
}

dfFinal <- dfFinal[,2:ncol(dfFinal)]
print(paste0("The dimensions of dfFinal is ", dim(dfFinal)))
print(SNPcounter)

# Calculate FGr
FGr_raw <- apply(dfFinal, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

## Scale by 1/sqrt(L-1)
L <- sum(SNPcounter)
print(L)
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
M <- nrow(df)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

## Save FGr
dfOut <- as.data.frame(cbind(df[,1], FGr))
fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")
