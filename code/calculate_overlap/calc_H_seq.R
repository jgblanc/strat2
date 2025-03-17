## Compute H using one SNP per independent block

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_H.R <prefix to plink files> <r prefix> <var prefix> <out prefix> <snp> <ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

plink_prefix = args[1]
r_prefix = args[2]
var_prefix = args[3]
out_prefix = args[4]
snp_file = args[5]
id_file = args[6]
out_file = args[7]


# Read in IDs
dfIDs <- fread(id_file)

# Set data collectors
M <- nrow(dfIDs)

# Read in all values of r and all GWAS variances

## First Chr
r_file_name <- paste0(r_prefix, 1, ".rvec")
var_file_name <- paste0(var_prefix,1, ".txt")
df <- fread(r_file_name)
df$CHR <- 1
dfVar <- fread(var_file_name)
df <- inner_join(df, dfVar)

## All other Chrs
for (i in 2:22) {

  # Read in R file
  r_file_name <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file_name)
  tmp$CHR <- i

  # Read in Var file
  var_file_name <- paste0(var_prefix,i, ".txt")
  dfVar <- fread(var_file_name)
  tmp <- inner_join(tmp, dfVar)

  # Combine
  df <- rbind(df, tmp)
}
print(paste0("There are ", nrow(df), " SNPs in all the R files"))

# Read in SNP file
dfSNP <- fread(snp_file) %>% select("ID")

# Combine SNP and R files
dfALL <- inner_join(dfSNP, df, by = "ID") %>% drop_na()
print(paste0("There are ", nrow(dfALL), " SNPs in all the R files combined with the pruned SNPs"))
L <- nrow(dfALL)

# Divide SNPs into 1000 equal blocks
block_size <- floor(nrow(dfALL) / 1000)
print(block_size)

# Standardize r values and divide by GWAS variance
dfALL$r <- dfALL$r / sd(dfALL$r)
print(paste0("The variance of r is ", var(dfALL$r)))
print(paste0("The mean of r is ", mean(dfALL$r)))
#dfALL$r <- dfALL$r / sqrt(dfALL$Var)

# Make a collector for all values of H
allH <- rep(NA, 2)
dfR <- dfALL
dfR$r <- dfR$r / sqrt(dfR$Var)

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
                      " --score ", tmp_r_name, " header-read center cols=dosagesum,scoresums --out ", out_prefix)


for (i in 1:length(allH)) {

  print(paste0("This is rep ", i))

  # Subset Rs and save
  dfR_tmp <- dfR %>% select("ID", "ALT", "r")
  print(paste0("The variance of r is ", var(dfR_tmp$r)))
  print(paste0("The mean of r is ", mean(dfR_tmp$r)))
  fwrite(dfR_tmp, tmp_r_name, quote = F, row.names = F, sep = "\t")

  # Set up plink command
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
  dfR$r <- dfR$r / sqrt(dfR$Var)

}

# Remove tmp files
#rm_cmd <- paste0("rm ", out_prefix, ".*")
#system(rm_cmd)

#dfOut <- as.data.frame(cbind(allH, rep(L, length(allH))))

# Calculate p-value
#realH <- allH[1]
#varH <- var(allH[2:length(allH)])
#se <- sqrt(varH)
#meanH <- mean(allH[2:length(allH)])
#pvalNorm <- pnorm(realH ,mean =meanH, sd = se, lower.tail = FALSE)
#pvalSim <- sum(realH >= allH[2:length(allH)]) / (length(allH) - 1)
#dfOut <- data.frame(H = realH, L = L, meanH = meanH, varH = varH,
#                    pvalNorm = pvalNorm, pvalSim = pvalSim)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")


