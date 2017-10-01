' run_contig_filter.R

Usage: run_contig_filter.R -i SEGFILE -o OUTPUT

Options:
    -i SEGFILE      Input file - seg file output from APOLLOH
    -o OUTPUT       Output file path - contig and filter applied to input
' -> doc

library(docopt)
args <- docopt(doc)

library(GenomicRanges)
library(hrdtools)
library(tidyverse)

run_contig_filter <- function(loh_path, genome_version = 'hg19') {
    loh_ranges <- import_ranges(loh_path, c('chr', 'start', 'end', 'width', 'nvar', 'copynumber', 'lohtype', 'major', 'minor'))
    loh_col = 'lohtype'
    cnv_col = 'copynumber'

    seqlevels(loh_ranges) = gsub('chr', '', as.character(seqlevels(loh_ranges)))
    seqlevels(loh_ranges) = gsub('23', 'X', as.character(seqlevels(loh_ranges)))
    seqlevels(loh_ranges) = gsub('24', 'Y', as.character(seqlevels(loh_ranges)))

    loh_ranges = loh_ranges[order(seqnames(loh_ranges), ranges(loh_ranges)), ]

    print('Running contig function...')
    loh_ranges <- hrdtools:::contig(loh_ranges, loh_col, cnv_col)    # Fill in seg gaps

    print('Running filter function. Printing number of events remaining after successive filtering iterations:')
    loh_ranges <- hrdtools:::filter_ranges(loh_ranges, loh_col, cnv_col)          # Filter out very short segs

    return(loh_ranges)
}

filtered <- run_contig_filter(args[['i']])
filtered %>% as.data.frame %>% rename(chr = seqnames) %>% dplyr::mutate(event_id = 1:n()) %>% write_tsv(args[['o']])
