suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

args = commandArgs(TRUE)

if(length(args) < 4){
  stop("Usage: Rscript residualize_one_by_one.R <snp_file> <pc_prefix> <chr_num> <out_file> <out_psam> <out_pvar>")
}

snp_file = args[1]
pc_prefix = args[2]
chr_num = as.numeric(args[3])
out_file = args[4]


#-------------------------------
# Select correct PCA file
#-------------------------------
if (chr_num %in% c(1,3,5,7,9,11,13,15,17,19,21)) {
  print("Odd Chr")
  pca_file_path_common <- paste0(pc_prefix, "even_PCA.eigenvec")
  pca_file_path_rare <- paste0(pc_prefix, "rare_even_PCA.eigenvec")
} else if (chr_num %in% c(2,4,6,8,10,12,14,16,18,20,22)) {
  print("Even Chr")
  pca_file_path_common <- paste0(pc_prefix, "odd_PCA.eigenvec")
  pca_file_path_rare <- paste0(pc_prefix, "rare_odd_PCA.eigenvec")
} else {
  stop("Chromosome number not recognized.")
}

dfPCs_common <- fread(pca_file_path_common)
colnames(dfPCs_common)[1] <- "FID"

dfPCs_rare <- fread(pca_file_path_rare)
colnames(dfPCs_rare)[1] <- "FID"

dfPCs <- inner_join(dfPCs_common, dfPCs_rare, by = c("IID", "FID"))
ind_ids <- dfPCs$FID


#-------------------------------
# Prepare covariates (with intercept)
#-------------------------------
covars <- as.matrix(dfPCs %>% select(starts_with("PC")))
covars_with_intercept <- cbind(1, covars)
print(paste0("The number of cols in covars_with_intercept ", ncol(covars_with_intercept)))


# Precompute QR decomposition for speed
qr_covars <- qr(covars_with_intercept)

print("Got covariates")

#-------------------------------
# Write header to output file
#-------------------------------
header <- c("SNP", "A1", "A2", ind_ids)
fwrite(as.list(header), out_file, col.names = FALSE, sep = "\t", quote = FALSE)


#-------------------------------
# Process SNP file line-by-line
#-------------------------------
con <- file(snp_file, open = "r")
readLines(con, n = 1, warn = FALSE)  # skip header

line_count <- 0
while(TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break  # EOF

  fields <- strsplit(line, "\\s+")[[1]]
  dosages <- as.numeric(fields[7:length(fields)])
  ID <- fields[2]
  REF <- fields[5]
  ALT <- fields[6]

  # Compute residuals using precomputed QR
  #resids <- resid(lm(dosages ~ covars))

  #fit <- lm.fit(x = covars_with_intercept, y = dosages)
  #resids <- fit$residuals

  fitted_vals <- qr.fitted(qr_covars, dosages)
  resids <- dosages - fitted_vals

  # Prepare and write row
  row_out <- data.table(ID, REF, ALT, t(resids))
  setnames(row_out, c("SNP", "A1", "A2", ind_ids))
  fwrite(row_out, out_file, append = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

  line_count <- line_count + 1
  #if (line_count %% 1000 == 0) {
  #  print(paste("Processed", line_count, "SNPs"))
  #}
  print(line_count)

}

close(con)
print(paste("Finished processing", line_count, "SNPs"))

