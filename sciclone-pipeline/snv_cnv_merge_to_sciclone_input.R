' snv_cnv_merge_to_sciclone_input.R

Usage: snv_cnv_merge_to_sciclone_input.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to input data - snv_cnv_merge output for signature timing.
    -o --output OUTPUT      Path to output SciClone-formatted input data.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

input <- read_tsv(args[['input']])

input %>%
    select(
        chr,
        pos,
        total_depth,
        var_reads = alt_depth
    ) %>%
    mutate(
        ref_reads = total_depth - var_reads,
        vaf = var_reads / total_depth * 100
    ) %>%
    select(
        chr, pos, ref_reads, var_reads, vaf
    ) %>%
    write_tsv(args[['output']])
