# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript compute_afr-wbs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# Read in allele frequencies
dfEAS <- fread(args[1])
dfNFE <- fread( args[2])
dfSAS <- fread(args[3])
dfAFR <- fread(args[4])
dfAMR <- fread(args[5])

# Go through all pairwise comparisons

# EAS - NFE
dfOut <- dfEAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfEAS$ALT_FREQS - dfNFE$ALT_FREQS
fwrite(dfOut,args[6] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS - SAS
dfOut <- dfEAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfEAS$ALT_FREQS - dfSAS$ALT_FREQS
fwrite(dfOut,args[7] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS - AFR
dfOut <- dfEAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfEAS$ALT_FREQS - dfAFR$ALT_FREQS
fwrite(dfOut,args[8] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS - AMR
dfOut <- dfEAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfEAS$ALT_FREQS - dfAMR$ALT_FREQS
fwrite(dfOut,args[9] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE - SAS
dfOut <- dfNFE %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfNFE$ALT_FREQS - dfSAS$ALT_FREQS
fwrite(dfOut,args[10] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE - AFR
dfOut <- dfNFE %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfNFE$ALT_FREQS - dfAFR$ALT_FREQS
fwrite(dfOut,args[11] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE - AMR
dfOut <- dfNFE %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfNFE$ALT_FREQS - dfAMR$ALT_FREQS
fwrite(dfOut,args[12] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS - AFR
dfOut <- dfSAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfSAS$ALT_FREQS - dfAFR$ALT_FREQS
fwrite(dfOut,args[13] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS - AMR
dfOut <- dfSAS %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfSAS$ALT_FREQS - dfAMR$ALT_FREQS
fwrite(dfOut,args[14] ,row.names=F,quote=F,sep="\t", col.names = T)

# AFR - AMR
dfOut <- dfAFR %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfAFR$ALT_FREQS - dfAMR$ALT_FREQS
fwrite(dfOut,args[15] ,row.names=F,quote=F,sep="\t", col.names = T)


