' get_mutation_burden.R

Usage: get_mutation_burden.R -p PATHS -o OUTPUT
' -> doc

library(docopt)
library(tidyverse)
library(stringr)

args <- docopt(doc)

get_mutation_burden <- function(paths) {
    sapply(paths, function(p) {
        vcf <- readVcf(p)
        length(rowRanges(vcf))
    })
}

library(VariantAnnotation)

tibble(paths = read_lines(args[['PATHS']])) %>%
    mutate(id = str_extract(paths, 'POG\\d\\d\\d'),
           sample_prefix = str_extract(paths, '[bioparch]+\\d+'),
           mutation_burden = get_mutation_burden(paths)) %>%
    write_tsv(args[['OUTPUT']])
