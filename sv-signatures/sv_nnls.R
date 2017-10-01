' Usage: sv_nnls.R -s SVPATH -r REFPATH -o OUTPUT [ -b ITERATIONS ]

Options:
    -s SVPATH       Path to SV catalogs file (.RData, list containing catalogs)
    -r REFPATH      Path to SV reference table (.tsv, with leading columns type, subtype, class)
    -o OUTPUT       Path to output exposures table (.tsv)
    -b ITERATIONS   Number of iterations to run for Monte Carlo simulation. Defaults to 1000.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(nnls)

print(args)

if (is.null(args[['b']])) {
    iterations = 1000
} else {
    iterations = as.numeric(args['b'])
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

sv <- readRDS(args[['s']])
ref <- read_tsv(args[['r']]) %>%
    mutate(clustered = grepl('clustered', type),
           subtype = factor(subtype, levels = c('1 kb', '10 kb', '100 kb', '1 Mb', '10Mb', 'tra')),
           SVTYPE = gsub('clustered_', '', type)) %>%
    select(clustered, SVTYPE, subtype, everything()) %>%
    select(-class, -type)

exposures <- plyr::ldply(sv, function(z) {
  z <- z %>%
    as_tibble %>%
    filter(SVTYPE == 'TRA' | subtype != 0) %>%
    mutate(subtype = factor(subtype, levels = c(1000, 10000, 100000, 1000000, 10000000, 0)),
           SVTYPE = as.character(SVTYPE))
  levels(z$subtype) <- c('1 kb', '10 kb', '100 kb', '1 Mb', '10Mb', 'tra')

  merged <- ref %>% left_join(z, by = c('clustered', 'SVTYPE', 'subtype'))

  A <- merged %>% select(-clustered, -SVTYPE, -subtype, -value) %>% as.matrix()
  b <- merged[['value']]

  nnls_out <- nnls_montecarlo(A, b, iterations)
})

write_tsv(exposures, args[['o']])
print(paste0('SV NNLS Exposures written to ', args[['o']]))
