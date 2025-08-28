## Compute Add together all LD blocks

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

FGr_file = args[1]
out_Fhat = args[2]
r_prefix = args[4]

# Read in Fmat
dfMat <- as.matrix(fread(FGr_file))
fhat_raw <- apply(dfMat, 1, sum)
print(length(fhat_raw))

# Read in r values

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


## calculate r^\top t
rTr <- t(df$r) %*% r
print(rTr)

## calculate \hat{f}
fhat<- fhat_raw / rTr

# FGR output
fwrite(as.data.frame(fhat), out_Fhat, quote = F, row.names = F, sep = "\t")
