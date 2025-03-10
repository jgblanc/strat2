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
r_prefix = args[4]
nsnp = as.numeric(args[5])
pc_snps = args[6]

# Read in IDs
dfSAM <- fread(paste0(plink_prefix, "1_v3.psam"))
M <- nrow(dfSAM)
dfFinal <- matrix(NA, nrow = M, ncol = 22)

## Read in R File
dfR <- fread(paste0(r_prefix, "1.rvec"))
for (i in 2:22) {
  tmp <- fread(paste0(r_prefix,i,".rvec"))
  dfR <- rbind(dfR, tmp)
}
print(paste0("There are ", nrow(dfR), " rows in the whole R file"))

## Combine with genotypes
dfFreq <- fread(paste0(plink_prefix, "1_v3.afreq"))
for (i in 2:22) {
  tmp <- fread(paste0(plink_prefix,i,"_v3.afreq"))
  dfFreq <- rbind(dfFreq, tmp)
}
dfR <- inner_join(dfR, dfFreq)
print(paste0("There are ", nrow(dfR), " rows in the whole R file"))
colnames(dfR)[1] <- "CHR"

## Combine with PC snps
dfSNPs <- fread(pc_snps)
dfR <- inner_join(dfR, dfSNPs)
print(paste0("There are ", nrow(dfR), " rows in the r combined with pc snp file"))

## Filter to SNPs with nonzero r values
dfR <- dfR %>% filter(!is.na(r) & r != 0) %>% sample_n(nsnp)
#dfR <- dfR %>% mutate(r = sample(r))
#dfR$r <- rnorm(nrow(dfR))
dfR$r <- scale(dfR$r)

## Print out L
L <- nrow(dfR)
print(paste0("L is ", L))

# Subset Rs and save
for (j in 1:22) {

  # Get SNPs on chr j
  dfR_tmp <- dfR %>% filter(CHR == j)  %>% select("ID", "ALT", "r")
  tmp_r_name <- paste0(out_prefix, j, ".rvec")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  if (nrow(dfR_tmp)== 0) {
     dfFinal[,j] <- as.matrix(rep(0, M))
     next
  }

  # Set up plink command
  plink_prefix_chr <- paste0(plink_prefix, j, "_v3")
  tmp_outfile <- paste0(out_prefix, j)
  plink_cmd <- paste0("plink2 --pfile ", plink_prefix_chr,
                      " --score ", tmp_r_name, " variance-standardize header-read cols=dosagesum,scoresums --out ", tmp_outfile)
  system(plink_cmd)

  # Read in plink output
  df<- fread(paste0(tmp_outfile, ".sscore"))
  rawFGr <- as.matrix(df[,3])
  dfFinal[,j] <- rawFGr

  # Remove tmp files
  rm_cmd <- paste0("rm ", tmp_outfile, ".*")
  system(rm_cmd)

}

print(head(dfFinal))

# Read in plink output
FGr_sum <- apply(dfFinal, 1, sum)
FGr <- FGr_sum * (1/sqrt(L-1))

# Print stats
varRaw <- var(FGr_sum)
paste0("The raw variance is ", varRaw)
varScale <- var(FGr)
paste0("The scales variance is ", varScale)
H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
print(paste0("H is ", H))
print(paste0("1/L is ", 1/L))
print(paste0("1/M is ", 1/M))
print(paste0("1/N is ", 1/1828))

# Save FGr of each block
dfOut <- as.data.frame(cbind(varRaw, varScale, H, (1/L), (1/M), (1/1/1828), nsnp))
colnames(dfOut) <- c("varRaw", "varScale", "H", "1/L", "1/M", "1/N", "L")
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
