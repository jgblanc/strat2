args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript block_nps.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

inFile = args[1]
outFile = args[2]

# Read in all pruned SNPs
df <- fread(inFile)
colnames(df) <- "ID"
print(head(df))

# Separate in to 1000 blocks
num_blocks <- 1000
block_size <- ceiling(nrow(df) / num_blocks)
df <- df %>% mutate(block = rep(1:num_blocks, each = block_size, length.out = n()))

# Save file
fwrite(df, outFile ,row.names=F,quote=F,sep="\t", col.names = F)



