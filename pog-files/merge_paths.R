' merge_paths.R 

Usage: merge_paths.R -o OUTPUT -c CNV -s CNVSEGS -l LOHSEGS -v SNV -m METADATA [ -u ]

Options:
    -o OUTPUT       Output path for merged data
    -c CNV          Copy number variants
    -s CNVSEGS      Copy number variant segs
    -l LOHSEGS      LOH output segs file
    -v SNV          Somatic small variants VCF
    -m METADATA     Metadata table
    -u              Unique values. For any POG IDs and sample prefixes with more than one row, select only the last row
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(lazyeval)

print(args)

process_paths_type_1 <- function(data, col) {
    # For paths of type: /projects/POG/POG_data/POG002/wgs/biop1_t_A10512_blood1_n_A10511/A10512_A10511/strelka/20773/bwa/results/passed.somatic.snvs.vcf

    arguments <- as.list(match.call())
    col <- eval(arguments$col, data)

    data %>% 
    mutate(pog_id = str_extract(col, 'POG[0-9]{3}'), 
           comparison = str_replace(col, '.*?wgs\\/(.*?)\\/.*', '\\1')) %>% 
    mutate(sample_prefix = str_extract(comparison, '[a-z0-9]+'), 
           sample_prefix_comparison = str_replace(comparison, '.*?_([bioparchblood0-9]+)_.*', '\\1'))
}

ploidy_map <- list('triploid' = 3,
                   'diploid' = 2,
                   'DIPLOID' = 2,
                   'tetraploid' = 4,
                   'diploid/triploid' = 2,
                   'diploid/tetratloid' = 2,
                   'pentaploid' = 5)

metadata <- read_tsv(args[['m']]) %>%
    filter(!is.na(biofx_ploidy) & !is.na(biofx_tc)) %>%
    mutate(pog_id = str_extract(path, 'POG[0-9]{3}')) %>%
    mutate(biofx_ploidy = if_else(grepl('diploid', tolower(biofx_ploidy)), 2,
        if_else(grepl('triploid', biofx_ploidy), 3,
        if_else(grepl('tetraploid', biofx_ploidy), 4,
        if_else(grepl('pentaploid', biofx_ploidy), 5, 0))))) %>%
    select(pog_id, sample_prefix, biofx_tc, biofx_ploidy)

cna_raw <- read_tsv(args[['c']], col_names = 'cna_raw') %>%
    process_paths_type_1(cna_raw)

cna_segs <- read_tsv(args[['s']], col_names = 'cna_segs') %>%
    process_paths_type_1(cna_segs)

loh_segs <- read_tsv(args[['l']], col_names = 'loh_segs') %>%
    process_paths_type_1(loh_segs)

snv <- read_tsv(args[['v']], col_names = 'snv') %>%
    process_paths_type_1(snv)

merged_paths <- cna_raw %>%
    inner_join(cna_segs, by = c('pog_id', 'comparison', 'sample_prefix', 'sample_prefix_comparison')) %>%
    inner_join(loh_segs, by = c('pog_id', 'comparison', 'sample_prefix', 'sample_prefix_comparison')) %>%
    inner_join(snv, by = c('pog_id', 'comparison', 'sample_prefix', 'sample_prefix_comparison')) %>%
    inner_join(metadata, by = c('pog_id', 'sample_prefix')) %>%
    distinct()

if (args[['u']]) {
    merged_paths <- merged_paths %>%
        group_by(pog_id, sample_prefix) %>%
        filter(row_number() == n()) %>%
        ungroup()
}

write_tsv(merged_paths, args[['o']])
