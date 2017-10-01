'multiple_biopsy_metadata.R - Analysis of POG flatfiles to identify cases with multiple biopsies

Usage: multiple_biopsy_metadata.R -p FLATFILEPATHS -o OUTPUT --primary-available --include POGID
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
    filter(!is.na(sample_collection_time)) %>%
    filter(
        grepl('-FFPE-', source) |
        grepl('-OCT-', source) | 
        grepl('-FF-', source) |
        grepl('-FA-', source) |
        grepl('-Archival-', source)
    ) %>%
    filter(grepl('WGS', protocol)) %>%
    group_by(id) %>%
    filter(length(unique(sample_collection_time)) > 1 | id[1] == args[['POGID']])

if (args[['--primary-available']]) {
    output <- output %>% 
        filter(
            (
                'Primary' %in% pathology_occurrence & 
                'Metastatic' %in% pathology_occurrence
            ) | 
            id[1] == args[['POGID']]
        )
}

message('done.')

write_tsv(output, args[['OUTPUT']])

message(paste0('Wrote output to ', args[['OUTPUT']]))
