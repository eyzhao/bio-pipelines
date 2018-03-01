' strand_bias_exposures.R

Usage: strand_bias_exposures.R -i INPUT -o OUTPUT [ --ref REF --ncores NCORES ]

Options:
    -i --input INPUT        Path to strand bias catalogs with columns mutation_type, strand, count
    -o --output OUTPUT      Path to output SignIT exposures object

    --ref REF               Path to reference signatures. If not provided, will use COSMIC 30 set
                            from signit::get_reference_signatures().
    --ncores NCORES         Number of cores to run MCMC across for SignIT exposures. Default: 1
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(signit)
library(foreach)

catalog <- read_tsv(args[['input']])

if (!is.null(args[['ref']])) {
    reference_signatures <- read_tsv(args[['ref']])
} else {
    reference_signatures <- get_reference_signatures()
}

if (is.null(args[['ncores']])) {
    n_cores = 1
} else {
    n_cores = as.numeric(args[['ncores']])
}

strand_chains <- foreach(strand_name = catalog$strand %>% unique) %do% {
    catalog %>%
        filter(strand == strand_name) %>%
        select(mutation_type, count) %>%
        get_exposures(reference_signatures, n_cores = min(n_cores, 4))
}

names(strand_chains) <- catalog$strand %>% unique

list(
    mutation_catalog = catalog,
    output = strand_chains
) %>%
    saveRDS(args[['output']])
