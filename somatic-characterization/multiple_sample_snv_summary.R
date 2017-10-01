' multiple_sample_snv_summary.R

Usage: multiple_sample_snv_summary.R -i INPUT -o OUTPUT [ -r REGEX ]

Options:
    -i --input INPUT            Comma-separated list of VCF files with SNV mutation data to merge.
    -o --output OUTPUT          Path to output file, a TSV reporting which mutations exist in which
                                    file(s).
    -r --regex REGEX            Regex string, where groups will be captured and concatenated with
                                    underscores to form the sample IDs. If not provided, the sample
                                    IDs will be filenames.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(stringr)
library(VariantAnnotation)

files <- strsplit(args[['input']], ',')[[1]]
if (! is.null(args[['regex']])) {
    sample_id = str_match(files, args[['regex']]) %>%
        as_tibble %>% 
        dplyr::select(-V1) %>% 
        unite(sample_id, everything()) %>%
        .$sample_id
} else {
    sample_id = files
}

vranges_list <- lapply(files, function(path) {
    message(sprintf('Reading file %s', path))
    readVcf(path)
})
names(vranges_list) <- sample_id

message('Processing VRanges')
vranges_list <- plyr::llply(vranges_list, function(vcf) {
    vr = rowRanges(vcf)
    vr$depth = geno(vcf)$DP[, 2]
    vr$A = geno(vcf)$AU[, 2, 1]
    vr$C = geno(vcf)$CU[, 2, 1]
    vr$G = geno(vcf)$GU[, 2, 1]
    vr$T = geno(vcf)$TU[, 2, 1]

    vr <- vr[vr$ALT %>% elementNROWS == 1, ]
    vr$ALT <- vr$ALT %>% CharacterList %>% unlist

    seqlevels(vr) = gsub('chr', '', as.character(seqlevels(vr)))
    seqlevels(vr) = gsub('23', 'X', as.character(seqlevels(vr)))
    seqlevels(vr) = gsub('24', 'Y', as.character(seqlevels(vr)))

    mutant_allele_depth = vr %>%
        as_tibble %>%
        dplyr::select(ALT, A, C, G, T) %>%
        mutate(i = row_number()) %>%
        gather(base, base_depth, A, C, G, T) %>%
        filter(ALT == base) %>%
        arrange(i) %>%
        .$base_depth

    vr$altDepth = mutant_allele_depth

    return(vr[, c('REF', 'ALT', 'depth', 'altDepth')])
})

message('Combining samples')
combined <- vranges_list %>% plyr::ldply(function(vr) {
    vr %>% as_tibble
}) %>%
    as_tibble %>%
    dplyr::select(
        sample = .id,
        chr = seqnames,
        pos = start,
        ref = REF,
        alt = ALT,
        altDepth,
        depth
    ) %>%
    arrange(
        chr, pos, ref, alt
    )

combined %>%
    write_tsv(args[['output']])
message(sprintf('Wrote combined table to %s', args[['output']]))
