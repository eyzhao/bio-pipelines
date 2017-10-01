'Usage: mutation_catalog.R -v VCFDATA -o OUTPUT
' -> doc

library('docopt')
library('SomaticSignatures')
library('BSgenome')

mutation_catalog <- function(vr) {
    genome <- getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    context <- mutationContext(vr, genome)
    catalog <- motifMatrix(context, group='sampleNames', normalize=FALSE)
    colnames(catalog) <- gsub('.*?(POG\\d\\d\\d\\_biop\\d).*', '\\1', colnames(catalog))
    return(catalog)
}

args <- docopt(doc)
vr <- readRDS(args[['VCFDATA']])
seqlevels(vr, force=TRUE) <- seqlevels(vr)[seqlevels(vr) != 'chrM']
catalog <- mutation_catalog(vr)
saveRDS(catalog, file=args[['OUTPUT']])
