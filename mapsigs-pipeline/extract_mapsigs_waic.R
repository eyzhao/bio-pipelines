' extract_mapsigs_waic.R

Usage: extract_mapsigs_waic.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to MapSigs RDS output file
    -o --output OUTPUT      Path to file to dump the WAIC metric
' -> doc

library(docopt)
args <- docopt(doc)

readRDS(args[['input']])
