' Usage: brca_config_table.R -v VCFGLOB -c CNVGLOB -m METAFILE -o OUTPUT
' -> doc

library('docopt')
library('plyr')
args <- docopt(doc)

vcf_paths <- Sys.glob(args[['VCFGLOB']])
cnv_paths <- Sys.glob(args[['CNVGLOB']])

vcf_df <- data.frame(id = gsub('.*(POG\\d\\d\\d\\_.*?)\\..*', '\\1', vcf_paths), vcf_paths)
cnv_df <- data.frame(id = gsub('.*(POG\\d\\d\\d\\_.*?)\\..*', '\\1', cnv_paths), cnv_paths)
merged <- merge(vcf_df, cnv_df, by='id')

meta <- read.table(args[['METAFILE']], sep='\t', header=T)
meta$id <- paste(meta$pog_id, meta$sample_prefix, sep='_')

final <- ddply(merge(merged, meta, by='id'), 'id', function(z) {
    data.frame(vcf_paths = z$vcf_paths[1],
               cnv_paths = z$cnv_paths[1],
               tumor_content = mean(as.numeric(z$biofx_tc))
               )
})

write.table(final, file = args[['OUTPUT']], sep='\t', col.names=T, row.names=F, quote=F)
