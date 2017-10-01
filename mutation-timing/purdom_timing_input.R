' purdom_timing_input.R

Generates the appropriate timing input for use in the package cancerTiming.

Usage: purdom_timing_input.R -v VCFPATH -l LOHPATH -o OUTPUT

Options:
    -v VCFPATH      Path to input VCF file containing SNV data
    -l LOHPATH      Path to file containing LOH segs after contig and filtering steps
    -o OUTPUT       Path to output table, ready for Purdom analysis

' -> doc

library(docopt)
args <- docopt(doc)

library(VariantAnnotation)
library(tidyverse)

import_vcf <- function(vcfPath, genomeVersion='hg19') {
    snvRanges <- readVcf(vcfPath, 'hg19')

    vrObj = rowRanges(snvRanges)
    vrObj$depth = geno(snvRanges)$DP[, 2]
    vrObj$A = geno(snvRanges)$AU[, 2, 1]
    vrObj$C = geno(snvRanges)$CU[, 2, 1]
    vrObj$G = geno(snvRanges)$GU[, 2, 1]
    vrObj$T = geno(snvRanges)$TU[, 2, 1]

    vrObj <- vrObj[elementNROWS(vrObj$ALT) == 1, ]
    vrObj$ALT = unlist(CharacterList(vrObj$ALT))

    seqlevels(snvRanges) = gsub('chr', '', as.character(seqlevels(snvRanges)))
    seqlevels(snvRanges) = gsub('23', 'X', as.character(seqlevels(snvRanges)))
    seqlevels(snvRanges) = gsub('24', 'Y', as.character(seqlevels(snvRanges)))

    return(vrObj)
}

loh_path <- args[['l']]
vcf_path <- args[['v']]
output_path <- args[['o']]

vcf <- import_vcf(vcf_path); vcf <- vcf[order(vcf), ]
loh <- read_tsv(loh_path, col_types = cols(chr = col_character())) %>% 
    filter(! lohtype %in% c('ALOH', 'DLOH', 'HOMD'), 
           copynumber <= 4 & copynumber > 1, 
           ! (lohtype == 'BCNA' & copynumber == 4)) %>%
    as.data.frame() %>%
    GRanges()
loh <- loh[order(loh), ]

overlaps <- findOverlaps(loh, vcf) %>% as_tibble()

loh_overlap <- loh[overlaps[['queryHits']], ] %>% as_tibble() %>% select(lohtype, copynumber, event_id)
vcf_overlap <- vcf[overlaps[['subjectHits']], ] %>% as_tibble()
vcf_loh <- cbind(vcf_overlap, loh_overlap) %>% as_tibble()

mutant_allele_depth = vcf_overlap %>%
    select(ALT, A, C, G, T) %>% 
    mutate(i = row_number()) %>% 
    gather(base, base_depth, A, C, G, T) %>% 
    filter(ALT == base) %>% 
    arrange(i) %>% 
    .$base_depth

purdom_input <- vcf_loh %>%
    select(chr = seqnames,
           pos = start,
           segId = event_id,
           depth, REF, ALT, copynumber, lohtype) %>%
    mutate(nMutAllele = mutant_allele_depth,
           nReads = depth,
           fraction = nMutAllele / nReads,
           mutationId = row_number()
           ) %>%
    mutate(type = case_when(.$lohtype == 'NLOH' ~ 'CNLOH', 
                            .$copynumber == 3 ~ 'SingleGain', 
                            .$copynumber == 4 ~ 'DoubleGain', 
                            .$copynumber == 2 ~ 'HET',
                            TRUE ~ 'NONE'))

purdom_input %>% write_tsv(output_path)

print(sprintf('Output written to %s', output_path))
