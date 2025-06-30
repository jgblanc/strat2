# This script computes the allele frequency difference between sardinians and all other european pops

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript compute_afr-wbs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# Read in allele frequencies
dfSDI <- fread(args[1])
dfEUR <- fread(args[2])


# Go through all pairwise comparisons

# SDI - EUR
dfOut <- dfSDI %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfSDI$ALT_FREQS - dfEUR$ALT_FREQS
fwrite(dfOut,args[3] ,row.names=F,quote=F,sep="\t", col.names = T)


