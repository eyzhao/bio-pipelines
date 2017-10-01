'
Usage: signeR-BIC.R ( -p PATHS | -G GLOB ) -o OUTPUT

Options:
    -p --paths PATHS            Comma-separated list of paths to signeR outputs serialized in .Rds files
    -G --glob GLOB              Globstring matching signeR output files serialized in .Rds files (i.e. "./signer/*.Rds")
    -o --output OUTPUT          Path to table with "tidy" BIC data ready to be visualized as a boxplot
' -> doc

library(docopt)
args <- docopt(doc)

library(signeR)
library(tidyverse)

if ( !is.null(args[['paths']]) ) {
    paths <- strsplit(args[['paths']], ',')[[1]]
} else if ( !is.null(args[['glob']]) ) {
    paths <- Sys.glob(args[['glob']])
} else {
    stop('Must provide either --glob or --paths')
}

signature_count_vector <- gsub('.*?(\\d+_signatures).*', '\\1', paths)
objects <- lapply(paths, function(f) { readRDS(f) })
names(objects) <- signature_count_vector

output <- plyr::ldply(objects, function(sig_obj) {
    return(data.frame(
        n_sig = sig_obj$Nsign,
        bic = sig_obj$BIC
    ))
}) %>% arrange(n_sig) %>% as_tibble

print(output)

output %>% write_tsv(args[['output']])
