' replication_strand_bias_catalogs.R

Usage: replication_strand_bias_catalogs.R -i INPUT -o OUTPUT [ --genome GENOME ]

Options:
    -i --input INPUT            Path to mutations as TSV, with cols chr, pos, ref, alt.
    -o --output OUTPUT          Path to output catalog
    --genome GENOME             Name of the BSgenome to use. By default, uses BSgenome.Hsapiens.UCSC.hg19.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(signit)

if (is.null(args[['genome']])) {
    genome_name = 'BSgenome.Hsapiens.UCSC.hg19'
} else {
    genome_name = args[['genome']]
}

genome <- getBSgenome(genome_name)

mutations <- read_tsv(args[['input']], col_types = cols(chr = col_character()))
    
mutations %>%
  filter(
    chr %in% as.character(1:22),
    ref %in% c('A', 'C', 'G', 'T')
  ) %>%
  mutate(
    strand = case_when(
        ref %in% c('G', 'A') ~ 'leading',
        ref %in% c('C', 'T') ~ 'lagging'
    )
  ) %>%
  plyr::ddply('strand', function(z) {
    with(z, mutations_to_catalog(chr, pos, ref, alt))
  }) %>%
  write_tsv(args[['output']])

