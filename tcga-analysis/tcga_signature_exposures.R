' tcga_signature_exposures.R - Computes mutation signature catalogs for many TCGA MAFs

Usage: tcga_signature_exposures.R -m MUTATIONS -o OUTPUT [ -c NCORES ]

Options:
    -m MUTATIONS    Path to file which contains the table of mutation catalogs (output by tcga_signature_catalogs.R)
    -o OUTPUT       Path to output file, into which mutation signature exposures will be deposited
    -c NCORES       Number of cores to use for parallel processing

It is assumed that the mutation catalog input contains two metadata columns cancer_type, sample.
The remaining 96 columns should denote mutation types.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(deconstructSigs)
library(doParallel)


### Register cores for multicore processing

if (!is.null(args[['c']])) {
    n_cores <- as.numeric(args[['c']])
    registerDoParallel(cores = n_cores)
    run_in_parallel = TRUE
}


### Read mutation catalogs and compute exposures

i = 0
catalogs <- read_tsv(args[['m']]) %>% arrange(sample)
n_samples <- dim(catalogs)[1]
exposures <- catalogs %>%
    plyr::ddply('sample', function(z) {
        cat(sprintf("\rCurrent sample: %s. Calculations %s percent complete  ", z$sample, round(i/n_samples * 100 * n_cores)))
        signatures <- z %>%
            select(-cancer_type) %>%
            column_to_rownames('sample') %>%
            whichSignatures(contexts.needed = TRUE,
                            tri.counts.method = 'exome2genome',
                            signatures.ref = signatures.cosmic)
        i <<- i + 1
        signatures$weights
    }, .parallel = run_in_parallel)
print('')

exposures <- exposures %>%
    left_join(catalogs %>% select(sample, cancer_type), by = 'sample') %>%
    select(cancer_type, sample, Signature.1:Signature.30) %>%
    arrange(cancer_type, sample)

exposures %>% write_tsv(args[['o']])
print(sprintf('Output exposures to %s', args[['o']]))
