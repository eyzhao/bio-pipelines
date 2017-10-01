' Usage: similarity_analysis.R -i INPUT -r REFERENCE -o OUTPUT [ -t TYPECOLNAME -s SUBTYPECOLNAME ]
' -> doc

library('docopt')
args <- docopt(doc)

library('readr')
library('tibble')
library('reshape2')

reference <- read_tsv(args[['REFERENCE']])

input <- read_tsv(args[['INPUT']])

if (is.null(args[['TYPECOLNAME']])) {
    type_colname <- 'type'
} else {
    type_colname <- args[['TYPECOLNAME']]
}

if (is.null(args[['SUBTYPECOLNAME']])) {
    subtype_colname <- 'subtype'
} else {
    subtype_colname <- args[['SUBTYPECOLNAME']]
}

merged <- as_tibble(merge(reference, input, by.x = c(type_colname, subtype_colname),
                          by.y = c('type', 'subtype')))

ref_col_indices <- which(grepl('Signature', names(merged)))
sample_col_indices <- which(grepl('V\\d', names(merged)))

m <- apply(merged[, ref_col_indices], 2, function(ref_col) {
    apply(merged[, sample_col_indices], 2, function(sample_col) {
        return(ref_col %*% sample_col / (sqrt(sum(ref_col^2)) * sqrt(sum(sample_col^2))))
    })
})

correlations <- melt(m)
colnames(correlations) <- c('sample_signature', 'reference_signatures', 'cosine_similiarity')

write_tsv(correlations, path = args[['OUTPUT']])
