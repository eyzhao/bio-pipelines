#!/usr/bin/env Rscript

' case_filter.R - Filters a set of paths from STDIN based on case list

Usage: case_filter.R -c CASELIST -o OUTPUT [ -i -n COLNAME ]

Options:
    -c CASELIST     File pointing to a list of case IDs (one per line, i.e. POG101_biop1). Filenames must contain unique case IDs.
    -o OUTPUT       Path to output resulting paths to
    -i              Optional: Include a column of matching IDs
    -n COLNAME      Optional: Name of the column to contain the paths (default: path)
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(dplyr))
library(tibble)
library(tidyr)
library(stringr)
library(readr)
args <- docopt(doc)

print('Reading from STDIN')
f <- file("stdin")
open(f)

print('Reading case list')
case_list <- tibble(id = read_lines(args[['c']])) %>%
    separate(id, into = c('pog_id', 'sample_prefix'))

print('Matching...')
df <- tibble(path = readLines(f)) %>%
    mutate(pog_id = str_extract(path, 'POG\\d\\d\\d')) %>%
    mutate(sample_prefix = str_extract(path, '[bioparch]+\\d')) %>%
    right_join(case_list, by = c('pog_id', 'sample_prefix')) %>%
    unite(sample_id, pog_id, sample_prefix)

if (!is.null(args[['n']])) {
    colnames(df)[colnames(df) == 'path'] <- args[['n']]
}

if (args[['i']]) {
    write_tsv(df[, c(2,1)], args[['o']])
} else {
    write_tsv(df[, 2], args[['o']])
}

print(paste0('Paths written to ', args[['o']]))
