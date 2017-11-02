'
Converts a table of merged mutation catalogs in long "tidy" format into a wide matrix ready to be parsed by
the mutSignatures (WTSI-style) mutation signature calculator.

Usage: catalog_tidy_to_mutSignatures.R -i CATALOGPATH -o OUTPUT [ -c COLUMN ]

Options:
    -i --input CATALOGPATH          Input path to catalog in long "tidy" format.
    -o --output OUTPUT              Output path to mutSignatures input .tsv file
    -c --column COLUMN              Name of the sample column (default: id_comparison)
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['column']])) {
    column_name = 'id_comparison'
} else {
    column_name = args[['column']]
}

library(tidyverse)

input <- read_tsv(args[['input']], col_types = cols(count = col_number()))

input <- input[, c(column_name, 'mutation_type', 'count')]

input %>%
    spread_(key_col = column_name, value_col = 'count') %>%
    as.data.frame %>%
    `rownames<-`(NULL) %>%
    column_to_rownames('mutation_type') %>%
    write.table(args[['output']], quote=FALSE, sep='\t')
