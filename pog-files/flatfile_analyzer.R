' Usage: flatfile_analyzer.R -c CASELIST -o OUTPUT
' -> doc

library('docopt')

args = docopt(doc)

cases <- read.table(args[['CASELIST']])[, 1]
flatfile_paths <- sprintf('/projects/POG/POG_data/%s/flatfile/%s.tab', cases, cases)
flatfile_paths <- flatfile_paths[file.exists(flatfile_paths)]

output <- lapply(flatfile_paths, function(path) {
                 print(path)
    table <- read.table(path, header=T, sep='\t', stringsAsFactors=F, fill=T, strip.white=T, quote="")

    if ('sample_prefix' %in% names(table)) {
        prefix = table$sample_prefix
    } else {
        prefix = rep('', dim(table)[1])
    }

    df <- data.frame(
        pog_id = gsub('(POG\\d\\d\\d).*', '\\1', table$source),
        sample_id = table$source,
        library_name = table$library_name,
        sample_prefix = prefix,
        path_tc = table$path_tc,
        biofx_tc = table$biofx_tc,
        ploidy = table$biofx_ploidy,
        physician = table$physician
    )

    return(df)
})

combined <- do.call('rbind', output)

write.table(combined, file=args[['OUTPUT']], col.names=T, row.names=F, quote=F, sep='\t')
