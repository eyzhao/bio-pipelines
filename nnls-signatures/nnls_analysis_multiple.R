' Usage: sv_nnls.R -i INPUT -r REFPATH -o OUTPUT [ -b ITERATIONS ]

Options:
    -i INPUT        Path to SNV catalogs file (.RData, list containing catalogs)
    -r REFPATH      Path to SNV signatures reference file
    -o OUTPUT       Path to output exposures table (.tsv)
    -b ITERATIONS   Number of iterations to run for Monte Carlo simulation. Defaults to 1000.
' -> doc

library('docopt')

args <- docopt(doc)

library(tidyverse)
library(nnls)

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

snv <- readRDS(args[['i']]) %>% t() %>% as.data.frame %>% rownames_to_column('sample') %>% as_tibble()
ref <- read_tsv(args[['r']]) %>% 
    mutate(class = paste(gsub('>', '', `Substitution Type`), paste(substr(Trinucleotide, 1, 1), substr(Trinucleotide, 3, 3), sep = '.'))) %>% 
    select(class, everything()) %>% 
    select(-`Substitution Type`, -Trinucleotide, -`Somatic Mutation Type`)
   
exposures <- plyr::ddply(snv, 'sample', function(z) {
    cat(paste0('\rCurrently processing: ', z[['sample']][1]))
    merged <- z[1, ] %>% 
              select(-sample) %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column('class') %>%
              rename(count = `1`) %>%
              inner_join(ref, by = 'class') %>%
              as_tibble()
    A <- merged %>% select(-class, -count) %>% as.matrix
    b <- merged[['count']]
    nnls_montecarlo(A, b)
})

write_tsv(exposures, args[['o']])

cat('\n')
print(paste0('SNV NNLS Exposures written to ', args[['o']]))
