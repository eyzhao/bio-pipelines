' aggregate_cluster_table.tsv

Usage: aggregate_cluster_table.tsv ( -p PATHS | -i INPUT | -G GLOB ) -o OUTPUT [ --trim-paths --index-col-name ICN ]

Options:
    -p --paths PATHS            Comma separated list of paths to TSV files
    -i --input INPUT            Path to a file containing paths to tables, one per line
    -G --glob GLOB              Globstring matching paths
    -o --output OUTPUT          Path to output
    --trim-paths                Includes only the file name instead of full path
    --index-col-name ICN        Name of the index column containing the paths. Defaults to "paths"
' -> doc

library(docopt)
library(plyr)
library(readr)
library(doParallel)
args <- docopt(doc)

if ( ! is.null(args[['paths']]) ) {
    paths = strsplit(args[['paths']], ',')[[1]]
} else if (! is.null(args[['input']])) {
    paths = readLines(args[['input']])
} else if ( ! is.null(args[['glob']]) ) {
    paths = Sys.glob(args[['glob']])
} else {
    stop('Must provide one of --paths, --input, or --glob.')
}

paths <- unique(paths)

message('Merging files')

registerDoParallel()

options(readr.show_progress = FALSE)

output <- ddply(data.frame(paths), 'paths', function(z) {
      path = as.character(z$paths)
      message(paste0('Reading file: ', path))

      file <- suppressMessages(
        read_tsv(
          path,
          skip = 1,
          col_names = c(
            'chr',
            'pos',
            'ref_depth',
            'var_depth',
            'vaf',
            'tumour_copy',
            'clean_copy',
            'depth',
            'adequate_depth',
            'cluster'
          )
        )
      )
}, .parallel = TRUE)

if (args[['trim-paths']]) {
    output <- dplyr::mutate(output,
        paths = sapply(as.character(paths), function(p) { tail(strsplit(p, '/')[[1]], 1) })
    )
}

if (is.null(args[['index-col-name']])) {
    icn = 'paths'
} else {
    icn = as.character(args[['index-col-name']])
}
colnames(output)[colnames(output) == 'paths'] <- icn

message('Done merging')

write_tsv(output, args[['output']])

message(paste0('Output written to ', args[['output']]))
