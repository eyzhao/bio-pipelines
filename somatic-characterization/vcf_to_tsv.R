' vcf_to_tsv.R

Usage: vcf_to_tsv.R -v VCF -o OUTPUT

Options:
    -v --vcf VCF        Path to input VCF
    -o --output OUTPUT  Path to output TSV. Default columns will be chr, pos, ref, alt.
' -> doc

library(docopt)
args <- docopt(doc)

library(VariantAnnotation)
library(tidyverse)

vcf <- readVcf(args[['vcf']]) %>%
    rowRanges()
vcf <- vcf[
    elementNROWS(vcf$REF) == 1 &
    elementNROWS(vcf$ALT) == 1
]
vcf$REF <- as.character(vcf$REF)
vcf$ALT <- unlist(CharacterList(vcf$ALT))

as.data.frame(vcf) %>%
    dplyr::select(
        chr = seqnames,
        pos = start,
        ref = REF,
        alt = ALT) %>%
    as_tibble() %>%
    write_tsv(args[['output']])
