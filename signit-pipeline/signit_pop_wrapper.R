' signit_pop_wrapper.R

Usage: signit_pop_wrapper.R -i INPUT -o OUTPUT [ --signit SIGNIT ]

Options:
    -i --input INPUT        Input mutation table
    -o --output OUTPUT      Path to output RDS file
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
    ) %>%
    mutate(
        mutation_type = get_snv_mutation_type(chr, pos, ref, alt)
    ) %>%
    filter(
        ! grepl('N', mutation_type)
    )

start_time <- Sys.time()
stan_object <- get_population_signatures(maf, method = 'vb')
end_time <- Sys.time()

output <- list(
    stan_obj = stan_object,
    run_time = as.double(end_time - start_time, units = 'secs')
)

saveRDS(output, args[['output']])
print(paste('Output to', args[['output']]))
