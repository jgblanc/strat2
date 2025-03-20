## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
out_prefix = args[3]
snp_file = args[4]
id_file = args[5]
out_file = args[6]
nsnp = as.numeric(args[7])


# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r

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

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID")

# Combine SNP and R files
dfALL <- inner_join(dfSNP, df, by = "ID") %>% drop_na() %>% sample_n(nsnp)
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)

# Divide SNPs into 1000 equal blocks
block_size <- floor(nrow(dfALL) / 1000)
print(block_size)

# Standardize r values and divide by GWAS variance
dfALL$r <- dfALL$r / sd(dfALL$r)
print(paste0("The variance of r is ", var(dfALL$r)))
print(paste0("The mean of r is ", mean(dfALL$r)))

# Make a collector for all values of H
allH <- rep(NA,1001)
dfR <- dfALL

# Subset SNP IDs
dfSNP_tmp <- dfR %>% select("ID")
snp_name <- paste0(out_prefix, ".snp")
fwrite(dfSNP_tmp, snp_name, quote = F, row.names = F, sep = "\t")

# Get freq file
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --extract ", snp_name ," --threads 8 ",
                      " --freq --out ", out_prefix)
system(plink_cmd)

# Set up plink command
freq_file <- paste0(out_prefix, ".afreq")
tmp_r_name <- paste0(out_prefix, ".rvec")
plink_cmd <- paste0("plink2 --pfile ", plink_prefix, " --keep ", id_file, " --extract ", snp_name ," --threads 8 --read-freq ", freq_file,
                      " --score ", tmp_r_name, " header-read variance-standardize cols=dosagesum,scoresums --out ", out_prefix)


for (i in 1:length(allH)) {

  # Subset Rs and save
  dfR_tmp <- dfR %>% select("ID", "ALT", "r")
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Run plink
  system(plink_cmd)

  # Read in plink output
  df<- fread(paste0(out_prefix, ".sscore"))
  print(head(df))
  rawFGr <- as.matrix(df[,3])

  # Calculate FGr
  print(paste0("The raw var is ", var(rawFGr)))
  FGr <- rawFGr * (1/(sqrt(L-1)))
  print(paste0("The scaled var is ", var(FGr)))

  # Calculate H
  H <- (1/(M * (L-1))) * (t(FGr) %*% FGr)
  print(paste0("H is ", H))
  allH[i] <- H

  # Shift the dfR dataframe
  dfR <- dfALL %>% mutate(r = c(tail(r, i * block_size), head(r, -i * block_size)))

}

# Remove tmp files
rm_cmd <- paste0("rm ", out_prefix, ".*")
system(rm_cmd)

# Compute output
trueH <- allH[1]

# Get variance of reps
reps <- allH[2:length(allH)]
varH <- var(reps)
meanH <- mean(reps)

# P-value from sims
pvalSim <- sum(reps > trueH) / length(reps)
pvalNorm <- pnorm(trueH, mean = mean(reps), sd=sqrt(varH), lower.tail = F)

# Construct output
dfOut <- data.frame(H = trueH, L = L, meanH = meanH, varH = varH, pvalNorm = pvalNorm, pvalSim = pvalSim)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")


