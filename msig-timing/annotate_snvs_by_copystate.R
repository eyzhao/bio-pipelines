' annotate_snvs_by_copystate_single.R

Usage: annotate_snvs_by_copystate_single.R -c CNVSEGS -s SNV -o OUTPUT

Options:
    -c CNVSEGS          CNV segs file
    -s SNV              VCF file containing SNV data
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

# args <- list(p = 'paths/w200-cna_loh_snv_merged.txt', o = 'rdata/snvs_with_copystate.R')

cnv_segs_path <- args[['c']]
snv_path <- args[['s']]
output_path <- args[['o']]

snv <- readVcf(snv_path)

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

cnv <- read_tsv(cnv_segs_path, col_names = c('chr', 'start', 'end', 'hmm')) %>%
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

snv_cnv %>% write_tsv(output_path)
print(sprintf('Saved output to %s', output_path))
