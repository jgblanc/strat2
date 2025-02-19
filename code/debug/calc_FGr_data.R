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
r_file = args[4]


## Read in R File
dfR <- fread(r_file)

## Combine with SNPs in GWAS panel
dfFreq <- fread(paste0(plink_prefix), ".afreq")
dfR <- inner_join(dfR, dfFreq)

## Filter to SNPs with nonzero r values
dfR <- dfR %>% filter(!is.na(r) & r != 0)
dfR$r <- scale(dfR$r)

## Print out L
L <- nrow(dfR)
print(paste0("L is ", L))

# Subset Rs and save
dfR_tmp <- dfR %>% filter(block == blockNum) %>% select("ID", "ALT", "r")
tmp_r_name <- paste0(out_prefix, ".rvec")
fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

# Set up plink command
plink_cmd <- paste0("plink2 --pfile ", plink_prefix,
                        " --score ", tmp_r_name, " variance-standardize header-read cols=dosagesum,scoresums --out ", out_prefix)
    system(plink_cmd)

# Read in plink output
df<- fread(paste0(tmp_outfile, ".sscore"))
rawFGr <- as.matrix(df[,3])
M <- nrow(df)

# Print stats
paste0("The raw variance is ", var(rawFGr))
paste0("The scales variance is ", var(rawFGr * (1/sqrt(L-1))))

# Save FGr of each block
dfOut <- as.data.frame(cbind(df[,1], rawFGr * (1/sqrt(L-1))))
fwrite(dfOut, FGr_outfile, quote = F, row.names = F, sep = "\t")


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
