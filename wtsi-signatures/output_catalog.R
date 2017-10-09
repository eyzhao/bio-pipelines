'Usage: output_catalog.R -c CATALOG -o OUTPUT -s STUDY
' -> doc

library('docopt')
library('R.matlab')

args <- docopt(doc)
catalog <- readRDS(args[['CATALOG']])

catalog <- catalog[order(rownames(catalog)), ]

types <- sapply(strsplit(rownames(catalog), ' '), function(z) {
    paste(substr(z[1], 1, 1), substr(z[1], 2, 2), sep='>')
})

subtypes <- sapply(strsplit(rownames(catalog), ' '), function(z) {
    paste0(substr(z[2], 1, 1), substr(z[1], 1, 1), substr(z[2], 3, 3))
})

cancerType <- as.character(args[['STUDY']])

originalGenomes <- as.matrix(catalog)
storage.mode(originalGenomes) <- "double" 

sampleNames <- colnames(catalog)

writeMat(con=args[['OUTPUT']],
         cancerType = cancerType,
         originalGenomes = originalGenomes,
         sampleNames = sampleNames,
         subtypes = subtypes,
         types = types)
