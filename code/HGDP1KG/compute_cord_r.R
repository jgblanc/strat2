## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript project_Tvec_chr.R <test panel prefix> <gwas panel prefix> <test vec file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

test_prefix = args[1]
tvec_file_lat = args[2]
tvec_file_long=args[3]
out_prefix = args[4]
overlap_snps = args[5]
outfile_lat = args[6]
outfile_long = args[7]

####################
## Functions #######
####################

# Compute t(X) %*% T
compute_b <- function(path_to_test, testvec_file, test_type, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  pheno_file <- testvec_file
  cmd_XT <- paste("sh code/HGDP1KG/compute_XT.sh", path_to_test, pheno_file, test_type, outfile_XT, overlap_snps, sep = " ")
  system(cmd_XT)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "xt_temp." , test_type ,".glm.linear"))
  count_plink <- fread(paste0(outpath, "xt_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  #nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  #counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  #mean_gc <- counts / nOBS
  #length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  #betas_plink_norm <- beta_plink$BETA * length_mc_genos

  #  Re-write .linear file with correct betas
  #beta_plink$BETA <- betas_plink_norm
  beta_reformat <- beta_plink %>% dplyr::select("#CHROM","ID", "REF", "ALT",  "BETA")
  beta_reformat[is.na(beta_reformat)] <- 0

  return(beta_reformat)
}


#####################
##     Main       ###
#####################


# Gather parameters
testID <- fread(paste0(test_prefix, ".psam"))
head(testID)
colnames(testID) <- c("IID",  "Sex")
n <- nrow(testID)

# Compute b for latitude
r_lat = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_lat, test_type = "lat", outpath = out_prefix)
r_lat$BETA  <- r_lat$BETA * (1/n)
colnames(r_lat)[5] <- "r"
print(head(r_lat))
fwrite(r_lat, outfile_lat, row.names = F, col.names = T, quote = F, sep = "\t")

# Compute b for latitude
r_long = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_long, test_type = "long", outpath = out_prefix)
r_long$BETA  <- r_long$BETA * (1/n)
colnames(r_long)[5] <- "r"
print(head(r_long))
fwrite(r_long, outfile_long, row.names = F, col.names = T, quote = F, sep = "\t")

