## Get Lat/Long Test Vector
## This script uses HGDP metadata and a psam to format test vectors for latitide and longitude

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript get_LatLong_Tvec.R <psam> <outfile Lat> <outfile Long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

infile_lat= args[1]
outfile_lat = args[2]
infile_long = args[3]
outfile_long = args[4]

# Read in lat
dfLat <- fread(infile_lat)
varTvec_lat <- var(dfLat$lat)

# Read in lat
dfLong <- fread(infile_long)
varTvec_long <- var(dfLong$long)

# Save output
fwrite(data.frame(varTvec = varTvec_lat),outfile_lat, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(data.frame(varTvec = varTvec_long),outfile_long, row.names = F, col.names = T, quote = F, sep = "\t")
