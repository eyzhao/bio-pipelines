' Usage: sv_nnls.R -i INPUT -r REFPATH -o OUTPUT [ -b ITERATIONS ]

Options:
    -i INPUT        Path to SNV catalogs file (.tsv, with columns "mutation_type" and "value")
    -r REFPATH      Path to SNV signatures reference file (.tsv)
    -o OUTPUT       Path to output exposures table (.tsv)
    -b ITERATIONS   The number of iterations to run for NNLS Monte Carlo (default: 1000)
' -> doc

library('docopt')

args <- docopt(doc)

library(tidyverse)
library(nnls)

if (is.null(args[['c']])) {
    iterations = 1000
} else {
    iterations <- as.numeric(args[['b']])
} 

if (iterations < 1) {
    stop('Number of iterations must be at least 1')
}

nnls_montecarlo <- function(A, b, iterations = 1000) {
    sig_names <- colnames(A)
    naive_solution <- nnls(A, b)$x
    random_draws <- rmultinom(iterations, sum(b), b)
    result <- apply(random_draws, 2, function(b_rand) { nnls(A, b_rand)$x })
    exposure_list <- apply(result, 1, function(x) {
        return(data.frame(
            nnls_mean = mean(x),
            nnls_lCI = sort(x)[round(iterations * 0.025)],
            nnls_uCI = sort(x)[round(iterations * 0.975)]
        ))
    })
    exposure_df <- do.call('rbind', exposure_list) %>%
        `rownames<-`(sig_names) %>%
        rownames_to_column('signature') %>%
        mutate(nnls_solution = naive_solution)
    return(exposure_df)
}

snv <- read_tsv(args[['i']])
ref <- read_tsv(args[['r']]) %>% 
    select(mutation_type = `Somatic Mutation Type`, everything()) %>% 
    select(-`Substitution Type`, -`Trinucleotide`)
   
merged <- snv %>% inner_join(ref, by = 'mutation_type') %>%
          as_tibble()
A <- merged %>% select(-mutation_type, -value) %>% as.matrix
b <- merged[['value']]
exposures <- nnls_montecarlo(A, b, iterations)

write_tsv(exposures, args[['o']])

print(paste0('SNV NNLS Exposures written to ', args[['o']]))
