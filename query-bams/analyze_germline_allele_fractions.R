' Usage: analyze_germline_allele_fractions.R -i INPUT -m METADATA -o OUTPUT

Options:
    -i INPUT        TSV file output from query_bams_by_row used on germline output
    -m METADATA     Path to metadata extracted from flatfiles
    -o OUTPUT       Path to output table
' -> doc

library('docopt')
args <- docopt(doc)

library('dplyr')

germline_bam_data <- read.table(args[['i']], sep='\t', header=T)
print(head(germline_bam_data))

by_variant <- group_by(germline_bam_data, tumour_id, chr, pos)
variant_data <- summarise(by_variant,
    ref = Ref[1],
    alt = Alt[1],
    normal_depth = depth[which(bam_type == 'normal')],
    tumour_depth = depth[which(bam_type == 'tumour')],
    normal_var_depth = var_depth[which(bam_type == 'normal')],
    tumour_var_depth = var_depth[which(bam_type == 'tumour')],
    gene = Gene[1],
    variant = Variant[1],
    clinvar = ClinVar[1],
    type = Type[1],
    flagged = Flagged[1],
    gmaf = as.numeric(gsub('GMAF=', '', GMAF[1]))
)

# Add Tumour Content Column to the Data

metadata <- read.table(args[['m']], sep='\t', header=T)
metadata <- metadata[metadata$sample_prefix == 'biop1', ]
metadata$tumour_id <- paste(metadata$pog_id, metadata$library_name, sep='_')

merged <- merge(variant_data, metadata[, c('tumour_id', 'biofx_tc', 'ploidy')], by='tumour_id')

write.table(merged, file=args[['o']], col.names=T, row.names=F, sep='\t', quote=F)
