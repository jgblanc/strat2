# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


## Load in IDs
dfIDs <- fread(args[1])
dfIDs <- dfIDs %>% select(FID, IID, pop)


# Go through all pairwise comparisons

# sdi-eur
df <- dfIDs %>% filter(pop %in% c("sdi", "eur")) %>% mutate(tvec = case_when(pop == "sdi" ~ 1, pop == "eur" ~ -1))
df$tvec<- scale(df$tvec)
print(paste0("The mean of Tvec is ", mean(df$tvec)))
print(paste0("The variance of Tvec is ", var(df$tvec)))
fwrite(df, args[2] ,row.names=F,quote=F,sep="\t", col.names = T)

















