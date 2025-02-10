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
dfFGr <- as.data.frame(dfFGr[,2:ncol(dfFGr)])

# Read in SNP file
dfSNP <- fread(SNP_file)

# Check variance of each block
vars <- rep(0, 574)
for (i in 1:ncol(dfFGr)) {
  L <- dfSNP$SNPcounter[i]
  vars[i] <- var(dfFGr[,i] * (1/(sqrt(L-1))))
}
mean(vars)

# Calculate FGr
FGr_raw <- apply(dfFGr, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Scale by 1/sqrt(L-1)
L <- sum(dfSNP$SNPcounter)
print(L)
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
M <- nrow(dfFGr)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

# Compute SE for D
nblocks <- ncol(dfFGr)
print(nblocks)
allHs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  mi <- as.numeric(dfSNP[i, 2])
  FGri <- dfFGr[,i] * (1/sqrt(mi-1))
  Hi <- (sum(FGri^2)) * (1/M) * (1 / (mi -1))
  allHs[i] <- (mi / (L - mi)) * (H - Hi)^2

}
varH <- mean(allHs)
se <- sqrt(varH)

pval <- pnorm(H ,mean =(1/L), sd = se, lower.tail = FALSE)
print(pval)


## Save FGr
#dfOut <- as.data.frame(cbind(df[,1], FGr))
#fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")
