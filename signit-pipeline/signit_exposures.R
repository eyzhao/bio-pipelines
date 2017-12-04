' signit_exposures.R

Usage: signit_exposures.R -c CATALOG -o OUTPUT [ -s SIGNIT ]

Options:
    -c --catalog CATALOG        Path to mutation catalog. TSV with two columns: mutation_type (str) and count (int).

    -o --output OUTPUT          Path to .Rds file containing the serialized exposures object

    -s --signit SIGNIT          Path to SignIT R library files. If not provided, will assume that SignIT is installed
                                and load the package using library(signit) instead.
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    library(devtools)
    library(tidyverse)
    library(rjags)
    library(nnls)
    library(dbscan)
    library(Rtsne)

    load_all(args[['signit']])
}

mutation_catalog <- read_tsv(args[['catalog']])

colnames(mutation_catalog) <- c('mutation_type', 'count')

exposures <- get_exposures(mutation_catalog)

print('SignIT Analysis Complete.')

saveRDS(exposures, file = args[['output']])

print(paste0('Results saved as RDS serialized object in ', args[['output']]))
