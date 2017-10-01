'aggregate_metadata.R - Analysis of POG flatfiles to identify cases with multiple biopsies

Usage: aggregate_metadata.R -p FLATFILEPATHS -o OUTPUT
' -> doc

library(docopt)
library(readr)
library(dplyr)

args <- docopt(doc)

flatfile_paths <- read_lines(args[['FLATFILEPATHS']])
metadata <- plyr::ddply(data.frame(path = flatfile_paths), 'path', function(p) {
    cat(paste0('\r Reading File: ', as.character(p$path)))
    suppressMessages(read_tsv(as.character(p$path)))
})

message('Aggregating and filtering data...')

output <- metadata %>%
    mutate(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', path)) %>%
    distinct()

message('done.')

write_tsv(output, args[['OUTPUT']])

message(paste0('Wrote output to ', args[['OUTPUT']]))
