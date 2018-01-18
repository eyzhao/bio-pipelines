' vcf_to_tsv.R - Converts VCFs into TSV format, with columns chr, pos, ref, alt

Usage: vcf_to_tsv.R -v VCF -o OUTPUT

Options:
    -v --vcf VCF        Path to variant call file
    -o --output OUTPUT  Path to output TSV file with columns chr, pos, ref, alt
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(VariantAnnotation)

vcf <- readVcf(args[['vcf']]) %>%
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

snv %>%
    write_tsv(args[['output']])
