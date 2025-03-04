args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript block_nps.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

inFile = args[1]
outFile = args[2]

# Read in all pruned SNPs
df <- fread(inFile)
colnames(df) <- "ID"

# Separate in to 1000 blocks
df <- df %>% mutate(block = ntile(ID, 1000))

# Save file
fwrite(df, outFile ,row.names=F,quote=F,sep="\t", col.names = F)



