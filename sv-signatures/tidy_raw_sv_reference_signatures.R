' tidy_raw_sv_reference_signatures.R

Usage: tidy_raw_sv_reference_signatures.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT            Path to raw SV signatures as provided in supplementary tables
                                of Nik Zainal et al. (2016).
                                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4910866/

    -o --output OUTPUT          Path to output SV signatures after parsing
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

remap_sv_types <- setNames(c('DEL', 'DUP', 'INV', 'TRA'), c('del', 'tds', 'inv', 'trans'))

reference_sv <- read_tsv(args[['input']]) %>% 
    separate(Size, c('sv_type', 'size'), sep = ':') %>%
    replace_na(list(size = '')) %>% 
    mutate(
        sv_type = sapply(sv_type, function(t) remap_sv_types[t])
    ) %>% 
    unite('mutation_type', Type, sv_type, size, sep = '|')

numeric_portion <- reference_sv %>%
    select(-mutation_type) %>% 
    apply(2, function(z) {as.numeric(gsub('%', '', z)) / 100.0}) %>%
    as_tibble

bind_cols(reference_sv %>% select(mutation_type), numeric_portion) %>%
    write_tsv(args[['output']])
