' filter_paths_by_cases.R

Usage: filter_paths_by_cases.R -i INPUT -c CASELIST -o OUTPUT

Options:
    -i INPUT        Path to file containing paths - must have columns pog_id and sample_prefix
    -c CASELIST     Path to file containing cases, one per line, i.e. POG144_biop1
    -o OUTPUT       Path to filtered output
' -> doc

library(docopt)
args <- docopt(doc)

library(dplyr)
library(readr)
library(tibble)
library(tidyr)

caselist <- tibble(id = read_lines(args[['c']])) %>%
    separate(id, c('pog_id', 'sample_prefix'))

read_tsv(args[['i']]) %>%
    inner_join(caselist, by = c('pog_id', 'sample_prefix')) %>%
    write_tsv(args[['o']])
