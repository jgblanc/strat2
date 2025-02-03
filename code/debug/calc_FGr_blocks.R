## Compute FGr for a single Chromosome with all SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_FGr_single_chr.R <prefix to plink files> <freq file> <r> <FGR>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

plink_prefix = args[1]
out_prefix = args[2]
r_outfile = args[3]
FGr_outfile = args[4]

## Read in R File
dfR <- fread(freq_file)

## Set up loop through blocks
nBlocks <- length(unique(dfR$block))
print(paste0("There are ", nBlocks, " blocks on the Chr"))
dfFGr_mat <- matrix(NA, nrow = 9999, ncol = nBlocks)

## Loop through blocks
for (i in 1:nBlocks) {

  # Block num
  blockNum <- unique(dfR$block)[i]
  print(blockNum)

  # Subset Rs and save
  dfR_tmp <- dfR %>% filter(block == blockNum)
  tmp_r_name <- paste0(out_prefix, blockNum, ".rvec")
  fwrite(dfR_tmp, r_outfile, quote = F, row.names = F, sep = "\t")

  # Set up plink command
  tmp_outfile <- paste0(out_prefix, blockNum)
  plink_cmd <- paste0("plink2 --pfile ", plink_prefix,
                      " --score ", tmp_r_name, " variance-standardize header-read cols=dosagesum,scoresums --out ", tmp_outfile)
  system(plink_cmd)

  ## Read in plink output
  df<- fread(paste0(tmp_outfile, ".sscore"))
  print(head(df))
  rawFGr <- as.matrix(df[,3])
  dfFGr_mat[,i] <- rawFGr
}


# Calculate FGr
FGr_raw <- apply(dfFGr_mat, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

## Scale by 1/sqrt(L-1)
L <- nrow(dfR)
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
