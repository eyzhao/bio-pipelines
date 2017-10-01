' Usage: filter_vranges.R -v VROBJ -l REGIONS -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(GenomicRanges)
library(VariantAnnotation)

vr <- readRDS(args[['VROBJ']])

regions_gr <- GRanges(read.table(args[['REGIONS']], header=F, sep='\t', col.names=c('chr', 'start', 'end', 'gene', 'ensg')))
seqlevels(regions_gr) <- paste0('chr', seqlevels(regions_gr)) 

overlaps <- findOverlaps(vr, regions_gr)
filtered <- vr[queryHits(overlaps), ]
values(filtered)$gene <- values(regions_gr)[subjectHits(overlaps), 'gene']

output <- as.data.frame(filtered, row.names=NULL)
saveRDS(output, file=args[['OUTPUT']])
