' mutation_set_analysis.R

Usage: mutation_set_analysis.R -x SETX -y SETY -o OUTPUT [ -c COLNAME --true --false ]

Options:
    -x --setx SETX          TSV mutation file with chr, pos, ref, and alt as columns

    -y --sety SETY          TSV mutation file with chr, pos, ref, and alt as columns

    -o --output OUTPUT      Path to output mutation table

    -c --colname COLNAME    Name of new boolean column to introduce indicating whether
                                mutations from SETX are also in SETY

    --true                  Filters the dataset so that rows of SETX are kept only if
                                key are in SETY

    --false                 Filters the dataset so that rows of SETX are kept only if
                                they are NOT in SETY
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

if (args[['true']] && args[['false']]) {
    stop('Both --true and --false cannot be provided at the same time')
}

read_mutation_table <- function(path) {
    read_tsv(
        path,
        col_types = cols(
            chr = col_character(),
            pos = col_number(),
            ref = col_character(),
            alt = col_character()
        )
    ) %>% 
    mutate(
        pos = as.integer(pos)
    )
}

if (is.null(args[['colname']])) {
    merged_colname = 'in_mutation_set'
} else {
    merged_colname = args[['colname']]
}

mutation_set_x <- read_mutation_table(args[['setx']])
mutation_set_y <- read_mutation_table(args[['sety']]) %>%
    select(
        chr, pos, ref, alt
    ) %>%
    mutate(in_set_y = TRUE)

joined_mutation_set <- mutation_set_x %>%
    left_join(mutation_set_y, by = c('chr', 'pos', 'ref', 'alt')) %>%
    mutate(in_set_y = ifelse(is.na(in_set_y), FALSE, in_set_y))

if (args[['true']]) {
    joined_mutation_set <- joined_mutation_set %>%
        filter(in_set_y)
} else if (args[['false']]) {
    joined_mutation_set <- joined_mutation_set %>%
        filter(! in_set_y)
}

colnames(joined_mutation_set)[colnames(joined_mutation_set) == 'in_set_y'] = merged_colname

write_tsv(joined_mutation_set, args[['output']])
