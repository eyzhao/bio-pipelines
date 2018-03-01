' strand_bias_summary.R

Usage: strand_bias_summary.R -i INPUT -o OUTPUT [ --signit SIGNIT ]

Options:
    -i --input INPUT        Path to strand-separated SignIT exposures
    -o --output OUTPUT      Path to output summary table
    --signit SIGNIT         Path to SignIT package (if not installed)
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    devtools::load_all(args[['signit']])
}

library(magrittr)
library(readr)

plyr::ldply(readRDS(args[['input']])$output, get_exposure_summary_table, .id = 'strand') %>%
    write_tsv(args[['output']])
