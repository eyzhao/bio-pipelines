'Usage: get_molecular_subtypes.R -p PATHSFILE -m METADATAFILE -o OUTPUT

Options:
    -p PATHSFILE        File with list of molecular subtype file paths, one per line
    -m METADATAFILE     Path to file containing metadata for POG cases
    -o OUTPUT           Path to output file
' -> doc

library('docopt')
library('plyr')

args <- docopt(doc)

metadata <- read.table(args[['-m']], header=T, stringsAsFactors=F, sep='\t')
paths <- read.table(args[['-p']], header=F, stringsAsFactors=F)[, 1]

library_id <- gsub('.*\\/([A-Z]\\d\\d\\d\\d\\d)\\/.*', '\\1', paths)

meta_merge <- merge(data.frame(library_id, paths), metadata,
                       by.x='library_id', by.y='library_name')

subtypes <- apply(meta_merge, 1, function(meta_row) {
    print(meta_row['pog_id'])
    t <- read.table(meta_row['paths'], header=T, sep='\t')
    cancer <- t[t$type == 'cancer', ]
    values <- ddply(cancer, 'diseases', function(z) {median(z$correlation)})
    med <- values[, 2]; names(med) <- values[, 1]
    return(med)
})

output = cbind(data.frame(pog_id = meta_merge$pog_id, sample_prefix = meta_merge$sample_prefix), t(subtypes))
output = output[order(output[, 1]), ]
write.table(output, file=args[['-o']], sep='\t', col.names=T, row.names=F, quote=F)
