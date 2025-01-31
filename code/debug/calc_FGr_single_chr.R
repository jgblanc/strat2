## Compute FGr for a single Chromosome with all SNPs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_FGr_single_chr.R <prefix to plink files> <freq file> <r> <FGR>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

plink_prefix = args[1]
freq_file = args[2]
r_outfile = args[3]
FGr_outfile = args[4]

## Read in freq file and add random r
dfFreq <- fread(freq_file)
dfFreq <- dfFreq %>% select("ID", "ALT")
dfFreq$r <- scale(runif(nrow(dfFreq)))

## Save R file
fwrite(dfFreq, r_outfile, quote = F, row.names = F, sep = "\t")

## Set up and plink command
plink_cmd <- paste0("plink2 --pfile ", plink_prefix,
                    " --score ", r_outfile, " variance-standardize header-read cols=dosagesum,scoresums --out ", plink_prefix)
system(plink_cmd)

## Read in plink output
df<- fread(paste0(plink_prefix, ".sscore"))
print(head(df))
rawFGr <- df[,3]

## Raw variance
print(paste0("The raw var is ", var(rawFGr)))

## Scale by 1/sqrt(L-1)
L <- nrow(dfFreq)
FGr <- rawFGr * (1/(sqrt(L-1)))
print(paste0("The scaled var is ", var(FGr)))

# Calculate H
M <- nrow(df)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))

## Save FGr
dfOut <- as.data.frame(cbind(df[,1], FGr))
fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")
