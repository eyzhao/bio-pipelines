' mutation_catalog.R

Uses SomaticSignatures to extract mutation catalogs for a sample. Accepts either VCF or TSV.

Usage: 
    mutation_catalog.R -v VCFDATA -o OUTPUT
    mutation_catalog.R -t TSVDATA -o OUTPUT [ --multisample ]

Options:
    -v VCFDATA              Variants in VCF form
    -t TSVDATA              Variants in TSV form. Required columns: chr, pos, ref, alt (additional columns ignored)
    -o OUTPUT               Path to output RData of mutation catalog

    -m --multisample        Report on multiple samples (requires additional column called "sample").
                                If this flag is not called, then the sample column is ignored and all
                                variants are collapsed together. 
                                NOTE: multisample currently only works with TSV input.
' -> doc

mutation_catalog <- function(snv) {
    catalog <- mut.to.sigs.input(snv %>% as.data.frame,
                      sample.id = 'sample',
                      chr = 'chr',
                      pos = 'pos',
                      ref = 'ref',
                      alt = 'alt')

    catalog %>%
        as.data.frame %>%
        rownames_to_column('sample') %>%
        gather(mutation_type, count, -sample) %>%
        mutate(count = as.numeric(count)) %>%
        as_tibble %>%
        arrange(sample, mutation_type)
}

library('docopt')

args <- docopt(doc)
print(args)

library('tidyverse')
library('deconstructSigs')

vcf_path <- args[['v']]
tsv_path <- args[['t']]
multisample <- args[['multisample']]

if (! is.null(vcf_path)) {
    library(VariantAnnotation)

    vcf <- readVcf(args[['v']]) %>%
        rowRanges()
    vcf <- vcf[
        elementNROWS(vcf$REF) == 1 &
        elementNROWS(vcf$ALT) == 1
    ]
    vcf$REF <- as.character(vcf$REF)
    vcf$ALT <- unlist(CharacterList(vcf$ALT))

    snv <- as.data.frame(vcf) %>%
        dplyr::select(
            chr = seqnames,
            pos = start,
            ref = REF,
            alt = ALT) %>%
        as_tibble()
} else if (! is.null(tsv_path)) {
    snv <- read_tsv(tsv_path) %>%
        mutate(chr = paste0('chr', gsub('chr', '', chr)))

    if (multisample && ! 'sample' %in% colnames(snv)) {
        stop('In order to run multisample, the file must contain a column called sample')
    }

    if (multisample) {
        snv <- snv %>% dplyr::select(sample, chr, pos, ref, alt)
    } else {
        snv <- snv %>% dplyr::select(chr, pos, ref, alt)
    }
}

print(snv)

snv <- snv %>%
    mutate(chr = gsub('chr', '', chr)) %>%
    filter(chr %in% as.character(1:22)) %>%
    mutate(chr = paste0('chr', chr))

if (! multisample) {
    snv <- snv %>% mutate(sample = 'x')
}

catalog <- mutation_catalog(snv)

if (! multisample) {
    catalog <- catalog %>% dplyr::select(-sample)
}

catalog %>% write_tsv(args[['o']])
