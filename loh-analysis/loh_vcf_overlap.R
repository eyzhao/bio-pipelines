#!/usr/bin/env Rscript

' loh_vcf_overlap.R - Identifies variant positions contained within LOH ranges

Usage: loh_vcf_overlap.R

' -> doc

library(docopt)
library(plyr)
library(hrdtools)
library(VariantAnnotation)
library(GenomicRanges)

get_files_table <- function() {
    loh_files = Sys.glob('data/POG*/*apolloh_out_segs.txt')
    germline_files = Sys.glob('data/POG*/germline_depth/*.bam.txt')
    pathogenic_files = Sys.glob('data/POG*/pathogenic_germline_depth/*.bam.txt')

    loh_ids = sapply(regmatches(loh_files, gregexpr('POG\\d+', loh_files)), function(z) {z[1]})
    germline_ids = sapply(regmatches(germline_files, gregexpr('POG\\d+', germline_files)), function(z) {z[1]})
    pathogenic_ids = sapply(regmatches(pathogenic_files, gregexpr('POG\\d+', pathogenic_files)), function(z) {z[1]})

    loh = data.frame(id=loh_ids, path=loh_files)
    germline = data.frame(id=germline_ids, path=germline_files)
    pathogenic = data.frame(id=pathogenic_ids, path=pathogenic_files)

    print('a')
    print(head(loh))
    print(head(germline))
    germline_merged = merge(loh, germline, by='id')
    pathogenic_merged = merge(loh, pathogenic, by='id')
    print('b')
    return(list(germline=germline_merged, pathogenic=pathogenic_merged))
}

import_snvs <- function(path, genomeVersion='hg19', hasHeader=FALSE) {
    table = read.table(path, header=hasHeader, stringsAsFactors=FALSE, sep='\t', fill=TRUE)
    table = table[table[, 1] != 'MT', ]
    
    chr = as.character(table[,1])

    ranges = IRanges(as.numeric(table[, 2]), as.numeric(table[, 2]))
    ref = as.character(table$Ref)
    alt = as.character(table$Alt)
    sampleNames = rep(path, length(chr))
    seqinfo = 'Genomes for HRD SNV signature analysis'

    vr <- VRanges(
        seqnames = as.character(chr),
        ranges = ranges,
        ref = ref,
        alt = alt,
        sampleNames = sampleNames,
        gene = as.character(table[, 6]),
        dbSNP = as.character(table[, 7]),
        zygosity = as.character(table[, 8]),
        impact = as.character(table[, 9]),
        cancertype = as.character(table[, 10])
    )

    return(vr)
}

has_overlaps <- function(lohPath, varPath) {
    loh <- import_ranges(lohPath)
    var <- import_snvs(varPath, hasHeader=TRUE)
    var_df <- read.table(varPath, header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
    overlaps = as.data.frame(findOverlaps(loh, var))

    if (dim(overlaps)[1] == 0) {
        return(FALSE)
    } else {
        lohRegions = loh[overlaps[, 1], ]
        if ('NLOH' %in% lohRegions$lohtype) {
            indices = which(lohRegions$lohtype == "NLOH")
            affectedGenes = var_df[overlaps[indices, 2], ]
            write.table(affectedGenes, paste0(varPath, '.cnlohoverlap'), quote=FALSE, sep="\t",
                        row.names=FALSE)
            print(paste0('Wrote to file ', varPath, '.cnlohoverlap'))
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
}

args <- docopt(doc)

files <- get_files_table()
germline = files$germline
pathogenic = files$pathogenic 

germline_overlaps <- apply(germline, 1, function(z) {has_overlaps(z[2], z[3])})
pathogenic_overlaps <- apply(pathogenic, 1, function(z) {has_overlaps(z[2], z[3])})

germline <- data.frame(id=germline[, 1], loh=germline[, 2], germline.overlap=germline_overlaps)
pathogenic <- data.frame(id=pathogenic[, 1], loh=pathogenic[, 2], pathogenic.overlap=pathogenic_overlaps)

output <- merge(germline, pathogenic, by='id')

countsTable <- ddply(output, "id", function(z) {
    id = z$id[1]
    g = z$germline.overlap           
    p = z$pathogenic.overlap
    return(data.frame(id=id, any.germline.overlap=any(g), any.pathogenic.overlap=any(p)))
})

print('Number of germline variants in CNLOH regions:')
print(table(countsTable[2])["TRUE"])
print('Number of pathogenic germline variants in CNLOH regions:')
print(table(countsTable[3])["TRUE"])

print('Total number of patients analyzed:')
print(dim(countsTable)[1])
