' collapse_multiple_sample_snv_summary.R

Takes a table of SNVs and "collapses" the sample IDs, such that the same
SNV observed in multiple samples will have a comma-separated sample value
containing both samples.

Usage: collapse_multiple_sample_snv_summary.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to TSV file containing mutations. One per line per sample.
                                Columns: sample, chr, pos, ref, alt.
                            Example:
                            sample       chr  pos      ref   alt
                            sample_a     1   54712     T     C   
                            sample_b     1   54712     T     C   
                            sample_a     1  726778     A     T   
                            sample_c     1  726778     A     T   

    -o --output OUTPUT      Path to output file with collapsed sample names. Will look like
                            sample                chr  pos      ref   alt
                            sample_a,sample_b     1   54712     T     C   
                            sample_a,sample_c     1  726778     A     T   
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

tsv_path <- args[['input']]

snv <- read_tsv(
        tsv_path, 
        col_types = cols(
            sample = col_character(),
            chr = col_character(),
            pos = col_number(),
            ref = col_character(),
            alt = col_character()
        )) %>% 
    dplyr::select(sample, chr, pos, ref, alt) %>%
    group_by(chr, pos, ref, alt) %>%
    summarise(sample = sample %>% unique %>% sort %>% paste(collapse = ',')) %>%
    ungroup() %>%
    write_tsv(args[['output']])
