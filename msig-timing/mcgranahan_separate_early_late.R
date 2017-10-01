' mcgranahan_separate_early_late.R

Given mutation data table with boolean column called "clonal" and mutation copy number column "mcn",
this script separates mutations into early and late classes and exports them to separate files.

Usage: mcgranahan_separate_early_late.R -i INPUT -e EARLY -l LATE [ -s ]

Options:
    -i INPUT            Input file - TSV containing tibble of variant data with cnv_state column
    -e EARLY            Output path for early mutations
    -l LATE             Output path for late mutations
    -s                  Stringent condition: clonal requires prob(CCF > 0.95) >= 0.75 and subclonal prob(CCF < 0.95) >= 0.75
' -> doc

library(docopt)

args <- docopt(doc)

library(tidyverse)

snv <- read_tsv(args[['i']]) %>%
    mutate(mcn = as.numeric(mcn),
           tcn = as.numeric(tcn),
           biofx_ploidy = as.numeric(biofx_ploidy),
           `p(CCF>95%)` = as.numeric(`p(CCF>95%)`)) %>%
    mutate(stringent_clonal = `p(CCF>95%)` >= 0.75,
           stringent_subclonal = `p(CCF<95%)` >= 0.75) %>%
    mutate(event_timing = if_else(tcn > biofx_ploidy, if_else(mcn == 1, 'after', 'before'), 'none'))

if (args[['s']]) {
    early <- snv %>%
        filter(clonal & stringent_clonal & event_timing != 'after')
    late <- snv %>%
        filter(((!clonal) & stringent_subclonal) | event_timing == 'after')
} else {
    early <- snv %>%
        filter(clonal & event_timing != 'after')
    late <- snv %>%
        filter((!clonal) | event_timing == 'after')
}

early %>% write_tsv(args[['e']])
late %>% write_tsv(args[['l']])
