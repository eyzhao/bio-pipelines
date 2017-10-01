' annotate_snvs_by_copystate_multiple.R

Usage: annotate_snvs_by_copystate_multiple.R -p MERGEDPATHS -o OUTPUT

Options:
    -p MERGEDPATHS      Table of paths with headers pog_id, comparison, cna_segs, cna_raw, loh_segs, snv
    -o OUTPUT           Path to output (saved as TSV)

' -> doc

library(docopt)
args <- docopt(doc)

library(VariantAnnotation)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(doParallel)

# args <- list(p = 'paths/w200-cna_loh_snv_merged.txt', o = 'rdata/snvs_with_copystate.R')

paths <- read_tsv(args[['p']])

registerDoParallel(cores = 31)

variant_depths <- paths %>%
    plyr::ddply(c('pog_id', 'comparison', 'sample_prefix', 'sample_prefix_comparison'), function(p) {
        print(p$snv)

        snv <- readVcf(p$snv[1])

        variant <- plyr::ldply(geno(snv)[5:8], function(z) {
            tibble(variant = rownames(z), tumour = z[, 'TUMOR', 1])
        }) %>%
            spread(.id, tumour) %>%
            separate(variant, into = c('chr', 'pos', 'ref', 'alt'), sep = '[:_\\/]') %>%
            rename(A = AU, C = CU, T = TU, G = GU)

        variant <- variant %>%
            mutate(total_depth = apply(variant, 1, function(row) {
                sum(as.numeric(c(row['A'], row['C'], row['T'], row['G'])))
            })) %>%
            mutate(alt_depth = as.numeric(apply(variant, 1, function(row) {
                idx <- which(colnames(variant) == row['alt']); return(row[idx]) 
            }))) %>% 
            as_tibble

        # overlap SNVs with CNV segs

        cnv <- read_tsv(p$cna_segs[1], col_names = c('chr', 'start', 'end', 'hmm')) %>%
            filter(hmm != 5) %>% # filter out high level amplifications as their exact CN values are unreliable.
            filter(chr %in% as.character(1:22)) %>% # filter out sex chromosomes as exact copy number estimates are trickier there
            GRanges()
        snv_gr <- GRanges(variant %>%
                          select(chr, start = pos, ref, alt) %>%
                          mutate(end = start) %>%
                          filter(chr %in% as.character(1:22)))
        overlaps <- findOverlaps(cnv, snv_gr)
        
        snv_cnv <- variant %>%
            as_tibble %>%
            slice(subjectHits(overlaps)) %>%
            cbind(cnv %>% as_tibble %>% slice(queryHits(overlaps)) %>% select(cnv_start = start, cnv_end = end, cnv_state = hmm)) %>%
            mutate(pos = as.numeric(pos)) %>%
            as_tibble

    }, .parallel = TRUE) 

variant_depths %>% write_tsv(args[['o']])
print(sprintf('Saved output to %s', args[['o']]))
