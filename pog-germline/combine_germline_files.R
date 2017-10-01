' Usage: combine_germline_files.R -p PATHS -o OUTPUT [ -n GENENAME ]
' -> doc

library(docopt)
library(readr)
library(plyr)
args = docopt(doc)

files <- readLines(args[['PATHS']])

combined <- ddply(data.frame(files), 'files', function(p) {
    a <- read_tsv(as.character(p$files[1]), comment = '##')
    exp_cols <- (grepl('percentile', colnames(a)) | grepl('Bodymap', colnames(a)) | grepl('targetExp', colnames(a)) | grepl('FC_', colnames(a)) )
    return(a[, !exp_cols]) 
})

if (!is.null(args[['GENENAME']])) {
    combined <- combined[combined$Gene == args[['GENENAME']], ]
}

colnames(combined) <- gsub('#', '', colnames(combined))

write_tsv(combined, args[['OUTPUT']])
