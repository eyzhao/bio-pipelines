' signit_pop_wrapper.R

Usage: signit_pop_wrapper.R -i INPUT -o OUTPUT [ -p NPOP -r REF --signit SIGNIT ]

Options:
    -i --input INPUT        Input mutation table
    -o --output OUTPUT      Path to output RDS file
    -p --npop NPOP          Number of populations. If not specified, will default to automated model selection.
    -r --ref REF            Path to reference signatures table. A matrix with columns mutation_type, ...
                                where ... refers to columns named based on the signature name.
                                mutation_type column is formatted for example C[C>A]T for a C>A mutation
                                in CCT context.
    --signit SIGNIT         Path to SignIT package for loading, if it is not installed in R
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    library(devtools)
    load_all(args[['signit']])
}

library(nnls)
library(tidyverse)
library(rjags)
library(deconstructSigs)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rstan)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

maf <- read_tsv(args[['input']]) %>%
    select(
        chr,
        pos,
        ref,
        alt,
        total_depth,
        alt_depth,
        tumour_copy,
        tumour_content,
        normal_copy
    ) %>%
    filter(
        tumour_copy > 0
    )

if (is.null(args[['npop']])) {
    n_pop = NULL
} else {
    n_pop = as.integer(args[['npop']])
}

if (is.null(args[['ref']])) {
    ref_signatures <- get_reference_signatures()
} else {
    ref_signatures <- read_tsv(args[['ref']])
}

print(ref_signatures)

start_time <- Sys.time()
stan_object <- get_population_signatures(maf, reference_signatures = ref_signatures, method = 'vb', n_populations = n_pop)
end_time <- Sys.time()

output <- list(
    stan_obj = stan_object,
    run_time = as.double(end_time - start_time, units = 'secs')
)

saveRDS(output, args[['output']])
print(paste('Output to', args[['output']]))
