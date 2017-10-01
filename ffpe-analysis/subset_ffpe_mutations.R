' subset_ffpe_mutations.R

Usage: subset_ffpe_mutations.R -f FFPE -s SNV -o OUTPUT

Options:
    -f --ffpe FFPE          Path to FFPE VCF file
    -s --snv SNV            Path to non-FFPE VCF file to filter with
    -o --output OUTPUT      Path to output (.tsv file) filtered FFPE variants
' -> doc

library(docopt)
args <- docopt(doc)

library(VariantAnnotation)
library(tidyverse)

ffpe <- readVcf(args[['ffpe']])
snv <- readVcf(args[['snv']])

overlaps <- findOverlaps(
    rowRanges(ffpe),
    rowRanges(snv)
)

filtered_ffpe <- rowRanges(ffpe)[overlaps %>% queryHits %>% unique]
filtered_ffpe <- filtered_ffpe[elementNROWS(filtered_ffpe$ALT) == 1]
filtered_ffpe$ALT <- filtered_ffpe$ALT %>% CharacterList %>% unlist
filtered_ffpe$REF <- filtered_ffpe$REF %>% as.character
filtered_ffpe %>% 
    as.data.frame %>%
    dplyr::select(
        chr = seqnames,
        pos = start,
        ref = REF,
        alt = ALT
    ) %>% 
    write_tsv(args[['output']])
