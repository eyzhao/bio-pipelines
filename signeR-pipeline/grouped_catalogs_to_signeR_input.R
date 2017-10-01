'
Converts a table of merged mutation catalogs in long "tidy" format into a wide matrix ready to be parsed by
the signeR mutation signature calculator.

Usage: catalog_tidy_to_mat.R -i CATALOGPATH -o OUTPUT

Options:
    -i --input CATALOGPATH          Input path to catalog in long "tidy" format.
    -o --output OUTPUT              Output path to WTSI input .mat file
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

input <- read_tsv(args[['input']], col_types = cols(value = col_number()))

input %>%
    select(id_comparison, mutation_type, value) %>%
    mutate(mutation_type = gsub('(.)\\[(.)>(.)\\](.)', '\\2>\\3:\\1\\2\\4', mutation_type)) %>%
    spread(mutation_type, value) %>%
    `rownames<-`(NULL) %>%
    column_to_rownames('id_comparison') %>%
    write.table(args[['output']], quote=FALSE, sep='\t')
