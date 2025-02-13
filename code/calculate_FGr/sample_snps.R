# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

outSample = args[1]
nsnp = strsplit(args[2], "L-")[[1]][2]

if (nsnp == "all") {
   df <- fread(args[25], header=FALSE)
   print(args[25])
   for (i in 26:length(args)) {
       tmp <- fread(args[i], header=FALSE)
       df <- rbind(df, tmp)
       }
   fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = F)
} else if (nsnp == "pruneall") {
   df <- fread(args[3], header=FALSE)
   for (i in 4:24) {
       tmp <- fread(args[i], header=FALSE)
       df <- rbind(df, tmp)
       }
   fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = F)
} else {
  nsnp = as.numeric(nsnp)
  df <- fread(args[3], header=FALSE)
  for (i in 4:24) {
      tmp <- fread(args[i], header=FALSE)
      df <- rbind(df, tmp)
      }
  df <- df %>% sample_n(nsnp)
  fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = F)
}



