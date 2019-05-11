' strand_bias_catalogs.R

Usage: strand_bias_catalogs.R -i INPUT -o OUTPUT -t TRANSCRIPTS [ --genome GENOME ]

Options:
    -i --input INPUT                Path to input TSV of mutations, with columns chr, pos, ref, and alt.
    -o --output OUTPUT              Path to output of strand-bias catalogs
    -t --transcripts TRANSCRIPTS    Path to transcripts from ENSEMBL biomart. This is a TSV file with
                                    at minimum the columns Chromosome Name, Transcript Start (bp),
                                    Transcript End (bp), and Strand.
    --genome GENOME                 Name of the BSgenome to use. By default, uses BSgenome.Hsapiens.UCSC.hg19.
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
transcripts <- read_tsv(args[['transcripts']])

transcripts_gr <- transcripts %>%
  filter(! is.na(`Ensembl Protein ID`)) %>%
  filter(`Chromosome Name` %in% as.character(1:22)) %>%
  select(
    chr = `Chromosome Name`,
    start = `Transcript Start (bp)`,
    end = `Transcript End (bp)`,
    strand = Strand
  ) %>%
  mutate(
    chr = paste0('chr', chr),
    strand = case_when(
      strand == 1 ~ '+',
      strand == -1 ~ '-',
      TRUE ~ 'none'
    )
  ) %>%
  filter(strand != 'none') %>%
  GRanges() %>%
  reduce() %>%
  .[strand(.) != '*']

bases <- c('A', 'C', 'G', 'T')

mutations <- read_tsv(args[['input']]) %>%
  filter(chr %in% 1:22) %>%
  mutate(
    chr = paste0('chr', chr),
    mutation_type = get_snv_mutation_type(chr, pos, ref, alt, genome),
    end = pos + 1
  ) %>%
  dplyr::rename(start = pos) %>%
  GRanges()

overlaps <- findOverlaps(
  transcripts_gr,
  mutations
)

merged <- bind_cols(
  transcripts_gr[queryHits(overlaps)] %>%
    as.data.frame %>%
    select(strand),
  mutations[subjectHits(overlaps)] %>%
    as.data.frame %>%
    select(
      chr = seqnames,
      pos = start,
      ref,
      alt,
      mutation_type
    )
) %>%
  as_tibble %>%
  mutate(strand = case_when(
    ref %in% c('G', 'A') & strand == '+' ~ 'transcribed',
    ref %in% c('C', 'T') & strand == '-' ~ 'transcribed',
    ref %in% c('G', 'A') & strand == '-' ~ 'coding',
    ref %in% c('C', 'T') & strand == '+' ~ 'coding'
  )) %>%
  mutate(chr = as.character(chr))

crossing(
    mutation_type = all_snv_mutation_types(),
    strand = c('transcribed', 'coding')
) %>%
    left_join(
        merged %>%
            group_by(strand, mutation_type) %>%
            summarise(count = n()),
        by = c('mutation_type', 'strand')
    ) %>%
    replace_na(list(count = 0)) %>%
    write_tsv(args[['output']])
