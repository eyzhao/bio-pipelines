'
Converts a table of merged mutation catalogs in long "tidy" format into a wide matrix ready to be parsed by
the WTSI matlab code.

Usage: catalog_tidy_to_mat.R -i INPUT -o OUTPUT 

Options:
    -i --input INPUT                Path to input catalog in .mat format.
    -o --output OUTPUT              Output path to tidy catalog table in TSV format
' -> doc

library(docopt)
library(R.matlab)
args <- docopt(doc)

library(tidyverse)

input <- readMat(args[['input']])

original_genomes <- input$originalGenomes
sample_names <- sapply(input$sampleNames, as.character)
types <- sapply(input$types, as.character)
subtypes <- sapply(input$subtypes, as.character)
mutation_types <- gsub('(.>.) (.).(.)', '\\2\\[\\1\\]\\3', paste(types, subtypes))

original_genomes %>%
    as_tibble %>% 
    `colnames<-`(sample_names) %>% 
    mutate(mutation_type = mutation_types) %>% 
    gather(sample, count, -mutation_type) %>%
    select(sample, mutation_type, count) %>%
    write_tsv(args[['output']])

message(paste('Output written to', args[['output']]))
