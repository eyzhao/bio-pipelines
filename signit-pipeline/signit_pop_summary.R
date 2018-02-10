' signit_pop_summary.R

Usage: signit_pop_summary.R -i INPUT -o OUTPUT [ --signit SIGNIT ]

Options:
    -i --input INPUT        Path to SignIT-Pop output. An RDS list object containing an item named
                            stan_obj, which holds the output from SignIT-Pop.

    -o --output OUTPUT      Path to output TSV file for writing SignIT Population Signature summary.

    --signit SIGNIT         Path to SignIT library, if not installed.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    devtools::load_all(args[['signit']])
}

signit_pop_data <- readRDS(args[['input']])
summarise_population_signatures(signit_pop_data$stan_obj) %>%
    mutate(
        n_mutations = signit_pop_data$stan_obj$mutation_table %>% nrow,
        n_populations = signit_pop_data$stan_obj$n_populations,
        bic = compute_population_signatures_bic(signit_pop_data$stan_obj),
        waic = compute_population_signatures_waic(signit_pop_data$stan_obj, parallel = TRUE, n_cores = 10)
    ) %>%
    write_tsv(args[['output']])
