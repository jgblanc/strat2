## Compute FGr for a single Chromosome with all SNPs

args=commandArgs(TRUE)

if(length(args)<1){stop("Rscript calc_FGr_single_chr.R <prefix to plink files> <freq file> <r> <FGR>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

prefix = args[1]

## Read in all block info

file_name <- paste0(prefix, 1, "_blocks.txt")
df <- fread(file_name)

for (i in 2:22) {

  # Read in file
  file_name <- paste0(prefix, i, "_blocks.txt")
  tmp <- fread(file_name)
  df <- rbind(df, tmp)

}

## Total number of SNPs
print(paste0("The total number of SNPs is ", nrow(df)))

## Add r
df$r <- scale(runif(nrow(df)))

## Save each chromosome separately

for (i in 1:22) {

  tmp <- df %>% separate(ID, c("CHR", "POS"), remove = FALSE) %>% filter(CHR == i)
  tmp <- tmp %>% select("ID", "ALT", "block", "r")
  file_name <- paste0(prefix,"r", i, "_standardize_blocks.rvec")
  fwrite(tmp, file_name, quote = F, row.names = F, sep = "\t")

}

