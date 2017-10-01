' compute_germline_zygosity.R -i GERMLINE -m METADATA -o OUTPUT
' -> doc

library('docopt')

args <- docopt(doc)

library('readr')

germline <- read_tsv(args[['GERMLINE']])
metadata <- read_tsv(args[['METADATA']])

metadata_germline <- 

merge(germline)
