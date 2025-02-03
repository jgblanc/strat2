## Compute FGr for a single Chromosome with all SNPs

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr_single_chr.R <prefix to plink files> <freq file> <r> <FGR>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

FGr_file = args[1]
SNP_file = args[2]
H_outfile = args[3]

# Read in FGR across blocks
dfFGr <- fread(FGr_file)

# Calculate FGr
FGr_raw <- apply(dfFGr, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Read in SNP file
dfSNP <- fread(SNP_file)

# Scale by 1/sqrt(L-1)
L <- sum(dfSNP$SNPcounter)
print(L)
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
M <- nrow(df)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

## Save FGr
#dfOut <- as.data.frame(cbind(df[,1], FGr))
#fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")
