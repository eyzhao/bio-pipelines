' split_snv_by_copystate.R

Usage: split_snv_by_copystate.R -i SNV -m METADATA -o OUTPUTDIR

Options:
    -i SNV          Input path - SNVs by copystate TSV to split
    -m METADATA     Metadata path
    -o OUTPUTDIR    Output directory to dump SNV files into
' -> doc

library(docopt)
args <- docopt(doc)

library(doParallel)
registerDoParallel(30)

library(tidyverse)

mcn_max_likelihood <- function(a, N, tcn, purity, ncn = 2) {
    mcn_domain <- 1:tcn
    likelihood <- sapply(mcn_domain, function(mcn) {
        theoretical_vaf <- purity * mcn / (ncn * (1-purity) + tcn * purity)
        dbinom(a, N, theoretical_vaf)
    })
    return(mcn_domain[which(likelihood == max(likelihood))][1])
}

mcn_max_likelihood_vector <- function(a, N, tcn, purity, ncn = 2) {
    if (length(ncn) == 1) { ncn = rep(ncn, length(a)) }
    sapply(1:length(a), function(i) {
        return(mcn_max_likelihood(a[i], N[i], tcn[i], purity[i], ncn[i]))
    }) %>% as.numeric()
}

metadata <- read_tsv(args[['m']]) %>%
     mutate(biofx_ploidy = if_else(grepl('diploid', biofx_ploidy), 'diploid', biofx_ploidy)) %>%
     filter(!is.na(biofx_tc) & !is.na(biofx_ploidy))

snv <- read_tsv(args[['i']], progress = TRUE) %>%
     filter(chr %in% (1:22 %>% as.character())) %>% # Retain only the autosomes
     inner_join(metadata %>% select(id, sample_prefix, biofx_tc, biofx_ploidy), by = c('pog_id' = 'id', 'sample_prefix')) %>%
     mutate(index = 1:length(pog_id),
            alt_depth = as.numeric(alt_depth),
            total_depth = as.numeric(total_depth),
            cnv_state = as.numeric(cnv_state),
            biofx_tc = as.numeric(biofx_tc) / 100,
            biofx_ploidy = factor(biofx_ploidy, levels = c('monoploid', 'diploid', 'triploid', 'tetraploid', 'pentaploid')) %>% as.numeric()
            ) %>%
     mutate(mcn = mcn_max_likelihood_vector(alt_depth, total_depth, cnv_state, biofx_tc, biofx_ploidy))

plyr::ddply(snv, c('pog_id', 'comparison'), function(snv_for_sample) { 
    path <- sprintf('meta/mcgranahan_input/%s_%s.txt', snv_for_sample$pog_id[1], snv_for_sample$comparison[1])
    snv_for_sample %>% write_tsv(path)
    print(paste('Output written to', path))
}, .parallel = T)
