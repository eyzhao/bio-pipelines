' export_wtsi_signatures_exposures.R

Usage: export_wtsi_signatures_exposures.R -p WTSIPATH -s SIGNATURES -e EXPOSURES 

Options:
    -p --path WTSIPATH          Comma-separated list of path
    -s --signatures SIGNATURES  Output path for signature table
    -e --exposures EXPOSURES    Output path for exposures table
' -> doc

library(docopt)
args <- docopt(doc)

library('R.matlab')
library('tidyverse')

get_exposure_table <- function(matlab_object) {
    sample_names <- unlist(matlab_object$input[, , 1]$sampleNames)
    exposures <- as.data.frame(t(matlab_object$exposures))
    exposures <- cbind(sample_names, exposures) %>% gather(signature, exposure, -sample_names)
    return(exposures)
}

get_signature_table <- function(matlab_object) {
    type <- unlist(matlab_object$input[, , 1]$types)
    subtype <- unlist(matlab_object$input[, , 1]$subtypes)
    class <- paste(type, gsub('(.).(.)', '\\1\\.\\2', subtype))

    signatures <- cbind(data.frame(type, subtype, class), as.data.frame(matlab_object$processes)) %>% 
        gather(signature, proportion, -type, -subtype, -class)
    return(signatures)
}

path <- strsplit(args[['paths']], ',')[[1]] %>% .[. != '']

wtsi_outputs <- plyr::dlply(data.frame(path), 'path', function(z) {
    file_path <- z$path
    output <- readMat(file_path)
})

wtsi_signatures <- plyr::ldply(wtsi_outputs, function(z) {
    get_signature_table(z)
}) %>%
    group_by(path) %>%
    mutate(n_signatures = signature %>% unique %>% length) %>%
    ungroup()

wtsi_exposures <- plyr::ldply(wtsi_outputs, function(z) {
    get_exposure_table(z)
}) %>%
    group_by(path) %>%
    mutate(n_signatures = signature %>% unique %>% length) %>%
    ungroup()

write.table(wtsi_signatures, file=args[['signatures']], col.names=T, row.names=F, sep='\t', quote=F)
write.table(wtsi_exposures, file=args[['exposures']], col.names=T, row.names=F, sep='\t', quote=F)
