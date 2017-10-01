' Sample Aggregation

Given a bunch of identically-formatted tables, one per sample, this script
aggregates the tables in long format, introducing a new column to account
for the path of the file.

Usage: aggregate_samples.R ( -p PATHS | -G GLOB ) -o OUTPUT [ -t ]

Options:
    -p --paths PATHS    Comma-separated list of file paths (i.e. path/to/file1,path/to/file2,path/to/file3)
    -G --glob GLOB      Globstring matching the desired input files.
    -o --output OUTPUT  Path to output file containing the aggregated table

    -t --trim           Trim the leading paths and file extension 
                            (i.e. remove everything before the last slash)

' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(doParallel)

registerDoParallel(30)

if (!is.null(args[['glob']])) {
    path <- Sys.glob(args[['glob']])
} else if (!is.null(args[['paths']])) {
    path <- strsplit(args[['paths']], ',')[[1]]
} else {
    stop('Must provide either --glob or --paths')
}

output <- plyr::ddply(data.frame(path), 'path', function(row) {
    suppressMessages(read_tsv(row[['path']][1] %>% as.character))
}, .parallel = TRUE) 

if (args[['trim']]) {
    output <- output %>% mutate(path = gsub('.*?([^/]+$)', '\\1', path))
}

output %>% write_tsv(args[['output']])
