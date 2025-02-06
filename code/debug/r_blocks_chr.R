## Add block info random rs

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript r_blocks_chr.R  <freq file> <ld blocks> <r>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))

freq_file = args[1]
ld_file = args[2]
r_outfile = args[3]


## Read in freq file
df <- fread(freq_file)
df <- df %>% separate(ID, c("CHR", "POS"), remove = FALSE)
print(paste0("There are ", nrow(df), " blocks before assigning info"))
print(head(df))

# Read in LD info
ld <- fread(ld_file)

# Convert types
df$POS <- as.numeric(df$POS)
df$CHR <- as.numeric(df$CHR)

# Assign SNPs to blocks
assign_SNP_to_block <- function(CHR, BP, block = ld) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Add block info - takes a while
small2 <- small %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[3]), as.numeric(params[4])))) %>%
  drop_na()
print(paste0("Now df blocks has", nrow(df), " rows"))
print(head(df))

# Format R file
dfOut <- df  %>% select("ID", "ALT", "block")

## Save R file
fwrite(dfOut, r_outfile, quote = F, row.names = F, sep = "\t")


