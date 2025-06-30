## Get Lat/Long Test Vector
## This script uses HGDP metadata and a psam to format test vectors for latitide and longitude

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_LatLong_Tvec.R <psam> <outfile Lat> <outfile Long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

md_file = args[1]
outfile_lat = args[2]
outfile_long = args[3]

# Read in metadata
df <- fread(md_file)

# Set up outfiles
dfLat <- df %>% select("FID","IID", "lat")
dfLong <- df %>% select("FID","IID", "long")

# Save output
fwrite(dfLat,outfile_lat, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfLong,outfile_long, row.names = F, col.names = T, quote = F, sep = "\t")
