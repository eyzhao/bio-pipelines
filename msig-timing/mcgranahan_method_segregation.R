' mcgranahan_method_segregation.R

Segregates mutations into clonal vs. non-clonal based on the analysis performed in
McGranahan et al. (2015): http://stm.sciencemag.org/content/7/283/283ra54.full

Usage: mcgranahan_method_segregation.R -i INPUT -m METADATA -s SAMPLEPREFIX -o OUTPUT [ -c CORES ]

Options:
    -i INPUT            Input file - TSV containing tibble of variant data with cnv_state column
    -m METADATA         Metadata file - Path to the flatfile
    -s SAMPLEPREFIX     Sample Prefix (i.e. biop1)
    -o OUTPUT           Output file
    -c CORES            Number of cores to run parallel processing with

Abbreviated variables:
    - vaf: variant allele fraction
    - ccf: cancer cell fraction
    - tcn: tumour copy number
    - ncn: normal copy number
    - mcn: mutation copy number, the number of alleles (<= tcn) harboring the mutation
    - purity: tumour content (cancer cellularity within the sample)
    - a: number of reads carrying the mutation
    - N: total number of reads at locus
' -> doc

library(docopt)
library(tidyverse)
library(doParallel)

cdf <- function(pdf) {
    sapply(1:length(pdf), function(i) {
        sum(pdf[1:i])
    })
}

vaf <- function(ccf, mcn, purity, tcn, ncn=2) {
    (purity * ccf * mcn) / ((1-purity) * ncn + purity * tcn)
}

vaf_table <- function(purity, mcn, tcn, ncn=2) {
    ccf_range = seq(0.01, 2, by = 0.01)
    vaf_vector <- sapply(ccf_range, function(ccf) {
        return(vaf(ccf, mcn, purity, tcn, ncn))
    })

    vaf_vector %>% 
        as_tibble %>%
        rename(vaf = value) %>%
        mutate(ccf = ccf_range) %>%
        mutate(tcn = tcn)
}

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

posterior <- function(a, N, mcn, tcn, purity, ncn=2) {
    vaf_table(purity, mcn, tcn, ncn) %>%
        mutate(posterior = sapply(vaf, function(v) { dbinom(a, N, prob = v) })) %>%
        mutate(posterior = posterior / sum(posterior))
}

register_cores <- function(number_of_cores) {
    if (is.null(number_of_cores)) {
        run_in_parallel = F
    } else {
        registerDoParallel(as.numeric(number_of_cores))
        run_in_parallel = T
    }

    return(run_in_parallel)
}

if (! interactive()) {
    args <- docopt(doc)

    run_in_parallel <- register_cores(args[['c']]) 

    metadata <- read_tsv(args[['m']]) %>%
        mutate(biofx_ploidy = if_else(grepl('diploid', biofx_ploidy), 'diploid', biofx_ploidy)) %>%
        filter(!is.na(biofx_tc) & !is.na(biofx_ploidy)) %>%
        filter(sample_prefix == args[['s']])

    tumour_content <- metadata$biofx_tc[1]
    ploidy <- metadata$biofx_ploidy[1]

    print('Reading SNV copystate file')

    snv <- read_tsv(args[['i']], progress = TRUE) %>%
        filter(chr %in% (1:22 %>% as.character())) %>% # Retain only the autosomes
        mutate(index = 1:length(chr),
               alt_depth = as.numeric(alt_depth),
               total_depth = as.numeric(total_depth),
               cnv_state = as.numeric(cnv_state),
               biofx_tc = tumour_content / 100,
               biofx_ploidy = factor(ploidy, levels = c('monoploid', 'diploid', 'triploid', 'tetraploid', 'pentaploid')) %>% as.numeric()
               ) %>%
        mutate(mcn = mcn_max_likelihood_vector(alt_depth, total_depth, cnv_state, biofx_tc, biofx_ploidy))

    print('Computing likelihoods')

    max_likelihood_snv <- plyr::ddply(snv, 'index', function (row) {
        post = posterior(row[['alt_depth']], row[['total_depth']], row[['mcn']], row[['cnv_state']], row[['biofx_tc']], row[['biofx_ploidy']]) %>%
            mutate(cumulative = cdf(posterior)) %>%
            mutate(`p(CCF>95%)` = sum(posterior[ccf > 0.95])) %>%
            mutate(`p(CCF<95%)` = sum(posterior[ccf < 0.95]))
        lCI = post %>% filter(abs(0.025 - cumulative) == min(abs(0.025 - cumulative))) %>% .$ccf
        uCI = post %>% filter(abs(0.975 - cumulative) == min(abs(0.975 - cumulative))) %>% .$ccf

        post %>% 
            filter(posterior == max(posterior)) %>%
            filter(row_number() == 1) %>%
            mutate(lCI = lCI, uCI = uCI)
    }, .parallel = run_in_parallel) %>%
        as_tibble()

    snv_ml <- snv %>%
        inner_join(max_likelihood_snv, by = 'index') %>%
        select(-posterior, -cumulative) %>%
        mutate(clonal = uCI > 1)

    write_tsv(snv_ml, args[['o']])

    print(sprintf('Wrote output to %s', args[['o']]))
}
