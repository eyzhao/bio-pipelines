#!/usr/bin/env Rscript

' run_purdom_analysis.R

Usage: run_purdom_analysis.R -i TIMINGINPUT -f FLATFILE -s SAMPLEPREFIX -o OUTPUT [ -b ITERATIONS ]

Options:
    -i TIMINGINPUT    Path to a file containing timing input table
    -f FLATFILE       Path to metadata (POG flatfile)
    -s SAMPLEPREFIX   Sample prefix (i.e. biop1)
    -o OUTPUT         Path to output (.Rds object)
    -b ITERATIONS     Number of iterations for bootstrapping cancerTiming confidence intervals (default: 1000)
' -> doc

library('docopt')
library(tidyverse)

event_timing <- function(inputTable, contamination, iterations = 1000) {
    library('cancerTiming')
    inputList <- list(inputTable)

    eventArgs = list()
    eventArgs['returnAssignments'] = TRUE
    eventArgs['bootstrapCI'] = 'nonparametric'
    eventArgs['B'] = iterations

    print(sprintf('Running mutation timing algorithm with %s iterations...', iterations))
    timingObject <- eventTimingOverList(inputList, contamination, eventArgs)
    return(timingObject)
}

args <- docopt(doc)

purity <- read_tsv(args[['f']]) %>% 
    filter(sample_prefix == args[['s']],
           !is.na(biofx_tc)) %>% 
    .$biofx_tc
contamination = 1 - (as.numeric(purity[1]) / 100.0)

timing_input = read_tsv(args[['i']], col_types = cols(chr = col_character()))

print(sprintf('Contamination: %s', contamination))
print('Timing Input:')
print(timing_input)

iterations = if_else(is.null(args[['b']]), 1000, as.numeric(args[['b']]))

timing_input %>% event_timing(contamination, iterations) %>% saveRDS(args[['o']])
