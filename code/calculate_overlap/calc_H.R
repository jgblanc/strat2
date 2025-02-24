## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_H.R <FGr file> <out File>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

inFile = args[1]
outFile = args[2]

# Read in FGr
dfFinal <- fread(inFile)
M <- nrow(dfFinal)

# Calculate FGr
FGr_raw <- apply(dfFinal, 1, sum)
print(paste0("The raw var is ", var(FGr_raw)))

# Scale by 1/sqrt(L-1)
L <- 1703
FGr <- FGr_raw * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

# Compute SE for H
nblocks <- ncol(dfFinal)
allHs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  FGri <- apply(dfFinal[,-i], 1, sum)
  FGri <- FGri * (1/sqrt(L-2))
  Hi <- (sum(FGri^2)) * (1/M) * (1/(L-1))
  allHs[i] <- (nblocks - 1) * (H - Hi)^2

}
print(allHs)

varH <- mean(allHs)
se <- sqrt(varH)

pval <- pnorm(H ,mean =(1/L), sd = se, lower.tail = FALSE)
print(pval)


## Save FGr
dfOut <- as.data.frame(cbind(H[1,], pval, var(FGr), varH))
colnames(dfOut) <- c("H", "pval", "VarFGr", "varH")
print(dfOut)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")
