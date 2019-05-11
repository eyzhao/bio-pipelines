' mutation_catalog.R

Uses SomaticSignatures to extract mutation catalogs for a sample. Accepts either VCF or TSV.

Usage: 
    mutation_catalog.R -v VCFDATA -o OUTPUT [ --genome GENOME ]
    mutation_catalog.R -t TSVDATA -o OUTPUT [ --multisample --colnames COLNAMES --genome GENOME ]

Options:
    -v VCFDATA              Variants in VCF form
    -t TSVDATA              Variants in TSV form. Required columns: chr, pos, ref, alt (additional columns ignored)

    --colnames COLNAMES     Optionally specify the four colnames instead of using chr, pos, ref, alt.
                                Provide these names as a comma-separated list.

    -o OUTPUT               Path to output RData of mutation catalog

    -m --multisample        Report on multiple samples (requires additional column called "sample").
                                If this flag is not called, then the sample column is ignored and all
                                variants are collapsed together. 
                                NOTE: multisample currently only works with TSV input.

    --genome GENOME         String matching a BSgenome name. If empty, will default to BSgenome.Hsapiens.UCSC.hg19.
' -> doc

mutation_catalog <- function(snv, genome) {
    print(snv)
    print(genome)

    if (nrow(snv) > 0) {
        catalog <- mut.to.sigs.input(
            snv %>% as.data.frame,
            sample.id = 'sample',
            chr = 'chr',
            pos = 'pos',
            ref = 'ref',
            alt = 'alt',
            bsg = genome
        )
        print(catalog)

        out <- catalog %>%
            as.data.frame %>%
            rownames_to_column('sample') %>%
            gather(mutation_type, count, -sample) %>%
            mutate(count = as.integer(count)) %>%
            as_tibble %>%
            arrange(sample, mutation_type)
    } else {
        bases = c('A', 'C', 'G', 'T')
        changes = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
        df = crossing(five_prime = bases, base_change = changes, three_prime = bases)
        mutation_types = sprintf('%s[%s]%s', df[['five_prime']], df[['base_change']], df[['three_prime']])
        out <- tibble(
            sample = 'x',
            mutation_type = mutation_types,
            count = 0
        )
    }

    return(out)
}

library('docopt')

args <- docopt(doc)
print(args)

library('tidyverse')
library('deconstructSigs')
library(BSgenome)

if (is.null(args[['genome']])) {
    genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
} else {
    genome = getBSgenome(args[['genome']])
}

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
    snv <- read_tsv(tsv_path, col_types = cols(.default = 'c'))

    if (! is.null(args[['colnames']])) {
        old_colnames = strsplit(args[['colnames']], ',')[[1]]
        print(paste('Using column name spec:', paste(old_colnames, collapse=', ')))
        stopifnot(length(old_colnames) == 4)
        
        new_colnames = c('chr', 'pos', 'ref', 'alt')
        old_colnames = setNames(old_colnames, new_colnames)
        snv <- snv %>% rename_(.dots = old_colnames)
    }

    snv <- snv %>%
        mutate(
            pos = as.integer(pos),
            chr = paste0('chr', gsub('chr', '', chr))
        )

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
    filter(chr %in% c(as.character(1:22), 'X', 'Y')) %>%
    mutate(chr = paste0('chr', chr))

if (! multisample) {
    snv <- snv %>% mutate(sample = 'x')
}

catalog <- mutation_catalog(snv, genome)

if (! multisample) {
    catalog <- catalog %>% dplyr::select(-sample)
}

catalog %>% write_tsv(args[['o']])
