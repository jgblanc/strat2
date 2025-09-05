suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

args = commandArgs(TRUE)

if(length(args) < 4){
  stop("Usage: Rscript residualize_one_by_one.R")
}

fgr_file = args[1]
pc_prefix = args[2]
out_file = args[3]
snp_file =args[4]


#------------------------------
# Read in GR Matrix
#-----------------------------

dfMat <- fread(fgr_file, drop = 1)
print(dim(dfMat))
dfSNP <- fread(snp_file)
chrs <- as.character(dfSNP$CHR)
colnames(dfMat) <- chrs


#----------------------------
# Residualize all GRs
#---------------------------
dfMat_resids <- matrix(NA, ncol = ncol(dfMat), nrow = nrow(dfMat))
print(dim(dfMat_resids))

for (i in seq_along(dfMat)) {

  print(i)
  chr_num <- as.numeric(colnames(dfMat)[i])
  print(chr_num)

  if (chr_num %% 2 == 1) {  # odd chromosome
    pca_file_path_common <- paste0(pc_prefix, "even_PCA.eigenvec")
    pca_file_path_rare   <- paste0(pc_prefix, "rare_even_PCA.eigenvec")
  } else {
    pca_file_path_common <- paste0(pc_prefix, "odd_PCA.eigenvec")
    pca_file_path_rare   <- paste0(pc_prefix, "rare_odd_PCA.eigenvec")
  }

  dfPCs_common <- fread(pca_file_path_common)
  setnames(dfPCs_common, 1:2, c("FID","IID"))

  dfPCs_rare <- fread(pca_file_path_rare)
  setnames(dfPCs_rare, 1:2, c("FID","IID"))

  dfPCs <- inner_join(dfPCs_common, dfPCs_rare, by = c("FID","IID"))
  covars_df <- dfPCs %>% select(starts_with("PC"))

  y <- as.numeric(dfMat[[i]])
  resids <- resid(lm(y ~ ., data = covars_df))
  dfMat_resids[,i] <- resids
}


#----------------------------
# Save output
#---------------------------
tmp <- fread(fgr_file)
dfOut <- as.data.frame(dfMat_resids)
dfOut <- cbind(IID = tmp[[1]], dfOut)
fwrite(dfOut, out_file, quote = F, row.names = F, sep = "\t")



