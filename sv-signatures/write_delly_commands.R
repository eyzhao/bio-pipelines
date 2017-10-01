' Usage: write_delly_commands.R -p PATHSFILE -o OUTPUT
' -> doc

library('plyr')
library('docopt')
args <- docopt(doc)

paths <- readLines(args[['PATHSFILE']])

library('VariantAnnotation')

library_names <- sapply(paths, function(p) {
    print(p)
    vcf <- scanVcfHeader(p)
    library_name <- paste(collapse='_', sapply(samples(vcf), function(z) { strsplit(z, '_')[[1]][[1]] }))
    return(library_name)
})

pog_ids <- gsub('.*?(POG\\d\\d\\d).*', '\\1', paths)

dir_commands <- paste0('mkdir -p data/', pog_ids)
destination <- paste0('data/', pog_ids, '/', pog_ids, '_', library_names, '.delly.vcf')
link_commands <- paste0('ln -sf ', paths, ' ', destination)[library_names != '']

writeLines(c(dir_commands, link_commands), con = args[['OUTPUT']])
