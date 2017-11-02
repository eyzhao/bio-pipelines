' row_bind_tables.R

Usage: row_bind_tables.R ( -p PATHS | -g GLOB ) -o OUTPUT

Options:
    -p --paths PATHS            Comma separated list of paths to TSV files
    -g --glob GLOB              Globstring matching paths
    -o --output OUTPUT          Path to output
' -> doc

library(docopt)
library(plyr)
library(readr)
library(doParallel)
args <- docopt(doc)

if ( is.null(args[['glob']]) ) {
    paths = strsplit(args[['paths']], ',')[[1]]
} else if ( is.null(args[['paths']]) ) {
    paths = Sys.glob(args[['glob']])
}

message('Merging files')

registerDoParallel()

options(readr.show_progress = FALSE)

output <- ddply(data.frame(paths), 'paths', function(z) {
      path = as.character(z$paths)
      message(paste0('Reading file: ', path))
      suppressMessages(read_tsv(path))
}, .parallel = TRUE)

message('Done merging')

write_tsv(output, args[['output']])

message(paste0('Output written to ', args[['output']]))
