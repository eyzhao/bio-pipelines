' annotate_snvs_by_copystate.R

Usage: annotate_snvs_by_copystate.R -c CNVSEGS -s SNV -p PLOIDY -t TC -o OUTPUT

Options:
    -c CNVSEGS          CNV segs file
    -s SNV              VCF file containing SNV data
    -p PLOIDY           Ploidy as an integer (i.e. 2 is diploid)
    -t TC               Tumour content as a fraction from 0 to 1
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

tumour_content <- as.numeric(args[['t']])
ploidy <- as.integer(args[['p']])

sample_prefix <- gsub('(biop|arch\\d+)_.*', '\\1', args[['a']])

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
    cbind(cnv %>% as_tibble %>% slice(queryHits(overlaps)) %>% select(cnv_start = start, cnv_end = end, tumour_copy = hmm)) %>%
    mutate(
        pos = as.numeric(pos),
        normal_copy = ploidy,
        tumour_content = tumour_content
    ) %>%
    as_tibble

snv_cnv %>% write_tsv(output_path)
print(sprintf('Saved output to %s', output_path))
