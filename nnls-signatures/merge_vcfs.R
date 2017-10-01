'Usage: merge_vcfs.R -v VCFPATHS -s STUDY -o OUTPUT
' -> doc

library('docopt')
library('VariantAnnotation')

import_vcf <- function(vcfPath, sampleName=NA, genomeVersion='hg19') {
    snvRanges <- readVcf(vcfPath, 'hg19')

    vrObj = rowRanges(snvRanges)

    # Remove rows where ALT has more than one character
    rowsToKeep = elementNROWS(vrObj$ALT) == 1
    vrObj <- vrObj[rowsToKeep, ]

    vrObj$REF = as.character(vrObj$REF)
    vrObj$ALT = unlist(CharacterList(vrObj$ALT))
    vrObj$totalDepth = geno(snvRanges)$DP[rowsToKeep, 2]

    depths = data.frame(A=geno(snvRanges)$AU[rowsToKeep, 2, 1],
                        C=geno(snvRanges)$CU[rowsToKeep, 2, 1],
                        G=geno(snvRanges)$GU[rowsToKeep, 2, 1],
                        T=geno(snvRanges)$TU[rowsToKeep, 2, 1]
                        )

    vrObj$refDepth = sapply(1:length(vrObj$REF), function(i) {depths[i, vrObj$REF[i]]} )
    vrObj$altDepth = sapply(1:length(vrObj$ALT), function(i) {depths[i, vrObj$ALT[i]]} )
    vrObj$sampleNames = rep(sampleName, length(vrObj$REF))

    if ('EFF' %in% names(info(snvRanges))) {
        vrObj$eff = info(snvRanges)[rowsToKeep, 'EFF']
    }

    return(makeVRangesFromGRanges(vrObj))
}

load_vcfs <- function(vcf_paths_file, study_name = 'POG') {
    paths <- as.character(read.table(vcf_paths_file, header=F)[, 1])

    vr_objects <- lapply(paths, function(p) {
        print(p)
        import_vcf(p, sampleName = p, 'hg19')
    })

    vr <- do.call('c', vr_objects)

    seqlevels(vr, force=TRUE) <- seqlevels(vr)[1:25]
    seqlevels(vr) <- paste0('chr', seqlevels(vr))
    seqlevels(vr)[25] <- 'chrM'

    values(vr)['study'] = study_name

    return(vr)
}

args <- docopt(doc)
vr <- load_vcfs(args[['VCFPATHS']], args[['STUDY']])
saveRDS(vr, file=args[['OUTPUT']])
