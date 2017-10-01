' Usage: copystate_tables.R -c CONFIGDATA -o OUTPUT

The CONFIGDATA file is a TSV table with headers id, vcf_paths, and cnv_paths (and other columns as desired).
OUTPUT is an RData file
' -> doc

library('docopt')

args <- docopt(doc)

library('devtools'); load_all('~/Projects/signatures/svn/mutation-timing/mustache/')

vcf_paths_df <- read.table(args[['CONFIGDATA']], header=T, stringsAsFactors=F, sep='\t', quote="")

copystates <- apply(vcf_paths_df, 1, function(path_line) {
    print(path_line['id'])
    vcf <- import_vcf(path_line['vcf_paths'], 'hg19')
    cnv <- import_cnv(path_line['cnv_paths'])
    copystate_table <- get_variant_copystate_table(vcf, cnv)
    return(copystate_table)
})

names(copystates) <- vcf_paths_df$id

saveRDS(copystates, file=args[['OUTPUT']])
