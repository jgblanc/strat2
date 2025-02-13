## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript standardized_r_IDs.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

r_prefix = args[1]
outfile = args[2]
snps = args[3]

num_chr=22

# Read in pruned SNPs
psnps <- fread(snps, header = FALSE)
colnames(psnps) <- "ID"
print(head(psnps))

# Read in all Rs
r_file <- paste0(r_prefix, "1.rvec")
r <- fread(r_file)
for (i in 2:num_chr) {
  r_file <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file)
  r <- rbind(r, tmp)
}
print(head(r))
print(tail(r))


# Join dataframes and standrdize everything
df <- inner_join(r, psnps)
df$r[is.na(df$r)] <- 0
df$r <- scale(df$r)

# Print relevant output
print(paste0("The var of r is ", var(df$r)))
print(paste0("L is ", nrow(df)))

# Save output
fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")

