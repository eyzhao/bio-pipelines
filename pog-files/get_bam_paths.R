' Usage: get_bam_paths.R -c CASELIST -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

flatfile_paths <- Sys.glob('/projects/POG/POG_data/POG???/flatfile/POG???.tab')
caselist <- readLines(args[['CASELIST']])
flatfile_paths <- flatfile_paths[gsub('.*(POG\\d\\d\\d).*', '\\1', flatfile_paths) %in% caselist]

flatfile_df <- do.call('rbind', lapply(flatfile_paths, function(path) {
    flatfile <- read.table(path, header=T, quote="", sep='\t', stringsAsFactors=F, fill=T)
    return(flatfile[, c('source', 'library_name', 'diseased_status', 'protocol', 'merged_bam', 'sample_prefix')])
}))

flatfile_df$pog_id <- sapply(strsplit(flatfile_df$source, '-'), function(z) {z[1]})

normal <- flatfile_df[flatfile_df$diseased_status == 'Normal'
    & grepl('WGS', flatfile_df$protocol)
    & ! grepl('FFPE', flatfile_df$source), c('pog_id', 'library_name', 'merged_bam', 'sample_prefix')]

tumour <- flatfile_df[flatfile_df$diseased_status == 'Diseased'
    & grepl('WGS', flatfile_df$protocol)
    & ! grepl('FFPE', flatfile_df$source), c('pog_id', 'library_name', 'merged_bam', 'sample_prefix')]

colnames(normal) <- c('pog_id', 'normal_library', 'normal_bam', 'normal_prefix')
colnames(tumour) <- c('pog_id', 'tumour_library', 'tumour_bam', 'tumour_prefix')
normal$normal_bam <- paste0(normal$normal_bam, sapply(normal$normal_bam, function(p) { dir(p, pattern='.bam$') }))
tumour$tumour_bam <- paste0(tumour$tumour_bam, sapply(tumour$tumour_bam, function(p) { dir(p, pattern='.bam$') }))

merged <- merge(normal, tumour, by='pog_id')

write.table(merged, file=args[['OUTPUT']], sep='\t', col.names=T, row.names=F, quote=F)
