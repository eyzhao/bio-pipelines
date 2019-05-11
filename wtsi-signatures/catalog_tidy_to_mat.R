'
Converts a table of merged mutation catalogs in long "tidy" format into a wide matrix ready to be parsed by
the WTSI matlab code.

Usage: catalog_tidy_to_mat.R -i CATALOGPATH -o OUTPUT -s STUDY [ -t TYPE -c IDCOL ]

Options:
    -i --input CATALOGPATH          Input path to catalog in long "tidy" format.
    -o --output OUTPUT              Output path to WTSI input .mat file
    -s --study STUDY                Name of the study, which will be loaded into cancerType variable
    -t --type TYPE                  Type of mutation signature. Can be either "snv" (default) or "sv".
                                        SNV catalogs mutation_type column is formatted like C[C>A]T
                                        denoting a C>A mutation in CCT context.
                                        SV catalogs mutation_type column is formatted like "clustered|DUP|10-100kb"
                                        or "non-clustered|TRA|".
    -c --idcol IDCOL                Name of the sample ID column (default is "sample")
' -> doc

library(docopt)
library(R.matlab)
args <- docopt(doc)

library(tidyverse)

input <- read_tsv(args[['input']])

if (!is.null(args[['idcol']])) {
    colnames(input)[ colnames(input) == args[['idcol']] ] = 'sample'
}

if (! 'sample' %in% colnames(input)) {
    stop('Column name sample is not found in input table. Please rename columns or provide --idcol parameter.')
}

matrix_df <- input %>%
    select(sample, mutation_type, count) %>%
    spread(sample, count)

originalGenomes <- matrix_df %>%
    select(-mutation_type) %>%
    data.matrix
colnames(originalGenomes) <- NULL
rownames(originalGenomes) <- NULL
storage.mode(originalGenomes) <- 'integer'

if (args[['type']] == 'sv') {
    types <- gsub('(.*?)\\|(.*?)\\|(.*)', '\\1', matrix_df$mutation_type)
    subtypes <- gsub('(.*?)\\|(.*?)\\|(.*)', '\\2\\|\\3', matrix_df$mutation_type)
} else if (args[['type']] == 'snv' || is.null(args[['type']])) {
    types <- gsub('.\\[(...)\\].', '\\1', matrix_df$mutation_type)
    subtypes <- gsub('(.)\\[(.)..\\](.)', '\\1\\2\\3', matrix_df$mutation_type)
}

cancerType <- as.character(args[['study']])
sampleNames <- colnames(matrix_df %>% select(-mutation_type))

writeMat(con = args[['output']],
         cancerType = cancerType,
         originalGenomes = originalGenomes,
         sampleNames = sampleNames,
         subtypes = subtypes,
         types = types)
