# This script computes the allele frequency difference between two contintental "populations"

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript compute_test_vec.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


## Load in IDs
dfIDs <- fread(args[11])
dfIDs <- dfIDs %>% select(FID, IID, pop)


# Go through all pairwise comparisons

# EAS-NFE
df <- dfIDs %>% filter(pop %in% c("eas", "nfe"))
nA <- nrow(df %>% filter(pop %in% c("eas")))
nB <-  nrow(df %>% filter(pop %in% c("nfe")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[1] ,row.names=F,quote=F,sep="\t", col.names = T)


# EAS-SAS
df <- dfIDs %>% filter(pop %in% c("eas", "sas"))
nA <- nrow(df %>% filter(pop %in% c("eas")))
nB <-  nrow(df %>% filter(pop %in% c("sas")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[2] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AFR
df <- dfIDs %>% filter(pop %in% c("eas", "afr"))
nA <- nrow(df %>% filter(pop %in% c("eas")))
nB <-  nrow(df %>% filter(pop %in% c("afr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[3] ,row.names=F,quote=F,sep="\t", col.names = T)

# EAS-AMR
df <- dfIDs %>% filter(pop %in% c("eas", "amr"))
nA <- nrow(df %>% filter(pop %in% c("eas")))
nB <-  nrow(df %>% filter(pop %in% c("amr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[4] ,row.names=F,quote=F,sep="\t", col.names = T)


# NFE-SAS
df <- dfIDs %>% filter(pop %in% c("nfe", "sas"))
nA <- nrow(df %>% filter(pop %in% c("nfe")))
nB <-  nrow(df %>% filter(pop %in% c("sas")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[5] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AFR
df <- dfIDs %>% filter(pop %in% c("nfe", "afr"))
nA <- nrow(df %>% filter(pop %in% c("nfe")))
nB <-  nrow(df %>% filter(pop %in% c("afr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[6] ,row.names=F,quote=F,sep="\t", col.names = T)

# NFE-AMR
df <- dfIDs %>% filter(pop %in% c("nfe", "amr"))
nA <- nrow(df %>% filter(pop %in% c("nfe")))
nB <-  nrow(df %>% filter(pop %in% c("amr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[7] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AFR
df <- dfIDs %>% filter(pop %in% c("sas", "afr"))
nA <- nrow(df %>% filter(pop %in% c("sas")))
nB <-  nrow(df %>% filter(pop %in% c("afr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[8] ,row.names=F,quote=F,sep="\t", col.names = T)

# SAS-AMR
df <- dfIDs %>% filter(pop %in% c("sas", "amr"))
nA <- nrow(df %>% filter(pop %in% c("sas")))
nB <-  nrow(df %>% filter(pop %in% c("amr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[9] ,row.names=F,quote=F,sep="\t", col.names = T)

# AFR-AMR
df <- dfIDs %>%  filter(pop %in% c("afr", "amr"))
nA <- nrow(df %>% filter(pop %in% c("afr")))
nB <-  nrow(df %>% filter(pop %in% c("amr")))
tvec <- c(rep(1, nA) , rep(-1, nB))
tvec <- scale(tvec)
print(paste0("The mean of Tvec is ", mean(tvec)))
print(paste0("The variance of Tvec is ", var(tvec)))
dfOut <- cbind(df[,1:2], tvec)
fwrite(dfOut, args[10] ,row.names=F,quote=F,sep="\t", col.names = T)



















