' signit_pairwise_bleeds.R

Usage: signit_pairwise_bleeds.R -e EXPOSURES -o OUTPUT [ -s SIGNIT ]

Options:
    -e --exposures EXPOSURES    Path to SignIT exposures, serialized as RDS files

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
    library(rstan)
    library(nnls)

    load_all(args[['signit']])
}

exposures <- readRDS(args[['exposures']])
get_exposure_pairwise_correlations(exposures) %>%
    mutate(path = args[['exposures']]) %>%
    write_tsv(args[['output']])

print(paste0('Results saved as tsv in ', args[['output']]))
